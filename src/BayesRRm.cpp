/*
 * BayesRRm.cpp
 *
 *  Created on: 5 Sep 2018
 *      Author: admin
 */

#include <cstdlib>
#include "BayesRRm.h"
#include "data.hpp"
#include "distributions_boost.hpp"
#include "options.hpp"
#include "samplewriter.h"
#include <chrono>
#include <numeric>
#include <random>
#include <algorithm>

#ifdef USE_MPI
#include <mpi.h>
#endif

BayesRRm::BayesRRm(Data &data, Options &opt, const long memPageSize)
: data(data)
, opt(opt)
, bedFile(opt.bedFile + ".bed")
, memPageSize(memPageSize)
, outputFile(opt.mcmcSampleFile)
, seed(opt.seed)
, max_iterations(opt.chainLength)
, burn_in(opt.burnin)
, thinning(opt.thin)
, dist(opt.seed)
, usePreprocessedData(opt.analysisType == "PPBayes")
, showDebug(false)
{
    float* ptr =static_cast<float*>(&opt.S[0]);
    cva = (Eigen::Map<Eigen::VectorXf>(ptr, static_cast<long>(opt.S.size()))).cast<double>();
}

BayesRRm::~BayesRRm()
{
}

void BayesRRm::init(int K, unsigned int markerCount, unsigned int individualCount)
{
    // Component variables
    priorPi = VectorXd(K);      // prior probabilities for each component
    pi = VectorXd(K);           // mixture probabilities
    cVa = VectorXd(K);          // component-specific variance
    cVaI = VectorXd(K);         // inverse of the component variances
    logL = VectorXd(K);         // log likelihood of component
    muk = VectorXd (K);         // mean of k-th component marker effect size
    denom = VectorXd(K - 1);    // temporal variable for computing the inflation of the effect variance for a given non-zero component
    m0 = 0;                     // total number of markers in model
    v = VectorXd(K);            // variable storing the component assignment

    // Mean and residual variables
    mu = 0.0;       // mean or intercept
    sigmaG = 0.0;   // genetic variance
    sigmaE = 0.0;   // residuals variance

    // Linear model variables
    beta = VectorXd(markerCount);           // effect sizes
    y_tilde = VectorXd(individualCount);    // variable containing the adjusted residuals to exclude the effects of a given marker
    epsilon = VectorXd(individualCount);    // variable containing the residuals

    //phenotype vector
    y = VectorXd();
    //SNP column vecotr
    Cx = VectorXd();

    // Init the working variables
    const int km1 = K - 1;
    cVa[0] = 0;
    cVa.segment(1, km1) = cva;

    //vector with component class for each marker
    components=VectorXd(markerCount);
    components.setZero();

    //set priors for pi parameters
    priorPi[0] = 0.5;
    priorPi.segment(1, km1) = priorPi[0] * cVa.segment(1, km1).array() / cVa.segment(1, km1).sum();

    y_tilde.setZero();

    cVaI[0] = 0;
    cVaI.segment(1, km1) = cVa.segment(1, km1).cwiseInverse();
    beta.setZero();

    //sample from beta distribution
    sigmaG = dist.beta_rng(1.0, 1.0);

    pi = priorPi;

    //scale phenotype vector stored in data.y
    y = (data.y.cast<double>().array() - data.y.cast<double>().mean());
    y /= sqrt(y.squaredNorm() / (double(individualCount - 1)));

    //initialize epsilon vector as the phenotype vector
    epsilon = (y).array() - mu;
    sigmaE = epsilon.squaredNorm() / individualCount * 0.5;
    betasqn=0;
    epsilonsum=0;
    ytildesum=0;
}


#ifdef USE_MPI

struct Lineout {
    double sigmaE,sigmaG;
    int    iteration, rank;
} lineout;

void sample_error(int error, const char *string)
{
  fprintf(stderr, "Error %d in %s\n", error, string);
  MPI_Finalize();
  exit(-1);
}

inline void scaadd(double* __restrict__ vout, const double* __restrict__ vin1, const double* __restrict__ vin2, const double dMULT, const int N) {

    if (dMULT == 0.0) {
        for (int i=0; i<N; i++) {
            vout[i] = vin1[i];
        }
    } else {
        for (int i=0; i<N; i++) {
            vout[i] = vin1[i] + dMULT * vin2[i];
        }
    }

}

inline void scaadd(double* __restrict__ vout, const double* __restrict__ vin2, const double dMULT, const int N) {

    if (dMULT == 0.0) {
        for (int i=0; i<N; i++) {
            vout[i] = 0.0;
        }
    } else {
        for (int i=0; i<N; i++) {
            vout[i] = dMULT * vin2[i];
        }
    }

}

inline double dotprod(const double* __restrict__ vec1, const double* __restrict__ vec2, const int N) {

    double dp = 0.0d;

    for (int i=0; i<N; i++) {
        dp += vec1[i] * vec2[i]; 
    }

    return dp;
}


int BayesRRm::runMpiGibbs() {

    const size_t LENBUF=200;
    const size_t ALIGN=1024;

    char buff[LENBUF]; 
    int  nranks, rank, name_len, result;

    MPI_Init(NULL, NULL);
    MPI_Comm_size(MPI_COMM_WORLD, &nranks);

    MPI_Comm_rank(MPI_COMM_WORLD, &rank);

    MPI_Status stat;
    MPI_Request request;

    // Set up processing options
    // -------------------------
    if (rank == 0) {
        opt.printBanner();
        opt.printProcessingOptions();
    }
    unsigned shuf_mark = opt.shuffleMarkers;
    unsigned sync_rate = opt.MPISyncRate;


    // Initialize MC on each worker
    // ----------------------------
    //cout << "Instantiating dist with seed = " << opt.seed + rank*1000 << endl;
    //Distributions_boost dist((uint)(opt.seed + rank*1000));
    dist.reset_rng((uint)(opt.seed + rank*1000), rank);

    const unsigned int max_it = opt.chainLength;
    const unsigned int N(data.numInds);
    unsigned int Mtot(data.numSnps);
    if (rank == 0)
        printf("Full dataset includes %d markers and %d individuals.\n", Mtot, N);
    if (opt.numberMarkers > 0 && opt.numberMarkers < Mtot)
        Mtot = opt.numberMarkers;
    if (rank == 0)
        printf("Option passed to process only %d markers!\n", Mtot);
    fflush(stdout);
    MPI_Barrier(MPI_COMM_WORLD);


    // Define global marker indexing
    // -----------------------------
    int MrankS[nranks];
    int MrankL[nranks];
    int modu  = Mtot % nranks;
    int Mrank = int(Mtot / nranks);
    int checkM = 0;
    int start = 0;

    // For production, should handle all markers regardless the number of tasks
#ifdef USEALLMARKERS
    for (int i=0; i<nranks; ++i) {
        MrankL[i] = int(Mtot / nranks);
        if (modu != 0 && i < modu)
            MrankL[i] += 1;
        MrankS[i] = start;
        //printf("start %d, len %d\n", MrankS[i], MrankL[i]);
        start += MrankL[i];
        checkM += MrankL[i];
    }
    assert(checkM == Mtot);
#else
    // Accept loosing M%nranks markers but easier to sync
    for (int i=0; i<nranks; ++i) {
        MrankL[i] = int(Mtot / nranks);
        MrankS[i] = start;
        //printf("start %d, len %d\n", MrankS[i], MrankL[i]);
        start += MrankL[i];
        checkM += MrankL[i];
    }
#endif

    if (rank == 0)
        printf("checkM vs Mtot: %d vs %d. Will sacrify %d markers!\n", checkM, Mtot, Mtot-checkM);

    int M = MrankL[rank];
    printf("rank %3d will handle a block of %6d markers starting at %d\n", rank, MrankL[rank], MrankS[rank]);


    const double	    sigma0 = 0.0001;
    const double	    v0E    = 0.0001;
    const double        s02E   = 0.0001;
    const double        v0G    = 0.0001;
    const double        s02G   = 0.0001;
    const unsigned int  K      = int(cva.size()) + 1;
    const unsigned int  km1    = K - 1;
    double              dNm1   = (double)(N - 1);
    double              dN     = (double) N;
    std::vector<int>    markerI;
    VectorXd            components(M);
    double              mu;              // mean or intercept
    double              sigmaG;          // genetic variance
    double              sigmaE;          // residuals variance
    VectorXd            priorPi(K);      // prior probabilities for each component
    VectorXd            pi(K);           // mixture probabilities
    VectorXd            muk(K);          // mean of k-th component marker effect size
    VectorXd            denom(K-1);      // temporal variable for computing the inflation of the effect variance for a given non-zero component
    VectorXd            cVa(K);          // component-specific variance
    VectorXd            cVaI(K);         // inverse of the component variances
    double              num;             // storing dot product
    //int                 m0;              // total number of markes in model
    double              m0;
    VectorXd            v(K);            // variable storing the component assignment
    VectorXd            sum_v(K);        // To store the sum of v elements over all ranks
    MatrixXd            beta(M,1);       // effect sizes
    //VectorXd            sample(2*M+4+N); // varible containg a sample of all variables in the model, M marker effects, M component assigned to markers, sigmaE, sigmaG, mu, iteration number and Explained variance

    // Length of a column in bytes
    const size_t snpLenByt = (data.numInds % 4) ? data.numInds / 4 + 1 : data.numInds / 4;
    if (rank==0)
        printf("snpLenByt = %zu bytes.\n", snpLenByt);

    // Open the bed file
    // Each MPI process will preprocess a section of the BED file, store it in RAM and process it
    // Let's assume for now that each node will be in charge of 10,000 SNPs and 500,000 INDs.
    // Buffer size for the raw   (char) data: 5.10^5 x 10^4 x 1/4 x 1 byte  = 1.25 x 10^9 bytes = 1.25 GB
    // Memory size for the final (dble) data: 5.10^5 x 10^4       x 8 bytes = 40.0 x 10^9 bytes = 40.0 GB

    MPI_File    bedfh, outfh;
    MPI_Status  status;
    std::string bedfp = opt.bedFile;
    bedfp += ".bed";
    result = MPI_File_open(MPI_COMM_WORLD, bedfp.c_str(), MPI_MODE_RDONLY, MPI_INFO_NULL, &bedfh);
    if(result != MPI_SUCCESS) 
        sample_error(result, "MPI_File_open bed file");

    std::string outfp = opt.mcmcSampleFile;

    result = MPI_File_open(MPI_COMM_WORLD, outfp.c_str(), MPI_MODE_CREATE | MPI_MODE_WRONLY | MPI_MODE_EXCL, MPI_INFO_NULL, &outfh);
    if (result != MPI_SUCCESS) {
        int lc = sprintf(buff, "FATAL: MPI_File_open failed to open file %s", outfp.c_str());
        sample_error(result, buff);
    }

    //int          blockcounts[2] = {2, 2};
    //MPI_Datatype types[2]       = {MPI_INT, MPI_DOUBLE}; 
    //MPI_Aint     displs[2];
    //MPI_Datatype typeout;
    //MPI_Get_address(&lineout.sigmaE,    &displs[0]);
    //MPI_Get_address(&lineout.iteration, &displs[1]);
    //for (int i = 1; i >= 0; i--)
    //    displs[i] -= displs[0];
    //MPI_Type_create_struct(2, blockcounts, displs, types, &typeout);
    //MPI_Type_commit(&typeout);


    const size_t rawdata_n = size_t(M) * size_t(snpLenByt)    * sizeof(char);
    const size_t ppdata_n  = size_t(M) * size_t(data.numInds) * sizeof(double);
    char*   rawdata   = (char*)  malloc(rawdata_n); if (rawdata == NULL) { printf("malloc rawdata failed.\n"); exit (1); }
    double* ppdata    = (double*)malloc(ppdata_n);  if (ppdata  == NULL) { printf("malloc ppdata failed.\n");  exit (1); }
    printf("rank %d allocation %zu bytes (%.3f GB) for the raw data.\n",           rank, rawdata_n, double(rawdata_n/1E9));
    printf("rank %d allocation %zu bytes (%.3f GB) for the pre-processed data.\n", rank, ppdata_n,  double(ppdata_n/1E9));


    // Compute the offset of the section to read from the BED file
    // -----------------------------------------------------------
    MPI_Offset offset = size_t(3) + size_t(MrankS[rank]) * size_t(snpLenByt) * sizeof(char);


    // Read the BED file
    // -----------------
    const auto st1 = std::chrono::high_resolution_clock::now();

    // Check how many calls are needed (limited by the int type of the number of elements to read!)
    uint nmpiread = 1;
    if (rawdata_n >= pow(2,(sizeof(int)*8)-1)) {   
        printf("MPI_file_read_at capacity exceeded. Asking to read %zu elements vs max %12.0f\n", 
               rawdata_n, pow(2,(sizeof(int)*8)-1));
        nmpiread = ceil(double(rawdata_n) / double(pow(2,(sizeof(int)*8)-1)));       
    }
    assert(nmpiread >= 1);
    cout << "Will need " << nmpiread << " calls to MPI_file_read_at to load all the data." << endl;

    if (nmpiread == 1) {
        result = MPI_File_read_at(bedfh, offset, rawdata, rawdata_n, MPI_CHAR, &status);
        if(result != MPI_SUCCESS) 
            sample_error(result, "MPI_File_read_at");
    } else {
        cout << "rawdata_n = " << rawdata_n << endl;
        size_t chunk    = size_t(double(rawdata_n)/double(nmpiread));
        size_t checksum = 0;
        for (int i=0; i<nmpiread; ++i) {
            size_t chunk_ = chunk;        
            if (i==nmpiread-1)
                chunk_ = rawdata_n - (i * chunk);
            checksum += chunk_;
            printf("rank %03d: chunk %02d: read at %zu a chunk of %zu.\n", rank, i, i*chunk*sizeof(char), chunk_);
            result = MPI_File_read_at(bedfh, offset + size_t(i)*chunk*sizeof(char), &rawdata[size_t(i)*chunk], chunk_, MPI_CHAR, &status);
            if(result != MPI_SUCCESS) 
                sample_error(result, "MPI_File_read_at");
        }
        if (checksum != rawdata_n) {
            cout << "FATAL!! checksum not equal to rawdata_n: " << checksum << " vs " << rawdata_n << endl; 
        }
    }
    const auto et1 = std::chrono::high_resolution_clock::now();
    const auto dt1 = et1 - st1;
    const auto du1 = std::chrono::duration_cast<std::chrono::milliseconds>(dt1).count();
    std::cout << "rank " << rank << ", time to read the BED file: " << du1 / double(1000.0) << " s." << std::endl;

    result = MPI_File_close(&bedfh);
    if(result != MPI_SUCCESS) 
        sample_error(result, "MPI_File_close");


    // Pre-process the data
    // --------------------
    const auto st2 = std::chrono::high_resolution_clock::now();
    data.preprocess_data(rawdata, M, snpLenByt, ppdata, rank);
    const auto et2 = std::chrono::high_resolution_clock::now();
    const auto dt2 = et2 - st2;
    const auto du2 = std::chrono::duration_cast<std::chrono::milliseconds>(dt2).count();
    std::cout << "rank " << rank << ", time to preprocess the data: " << du2 / double(1000.0) << " s." << std::endl;
    

    for (int i=0; i<M; ++i)
        markerI.push_back(i);
    //std::iota(markerI.begin(), markerI.end(), 0);

    //data.ZPZdiag.resize(M);


    // Processing part
    // ---------------
    const auto st3 = std::chrono::high_resolution_clock::now();

    double *y, *y_tilde, *epsilon, *tmpEps, *deltaEps, *dEpsSum, *deltaSum, *Cx;

    const size_t NDB = size_t(N) * sizeof(double);

    y        = (double*)malloc(NDB); if (y        == NULL) { printf("malloc y failed.\n");        exit (1); }
    y_tilde  = (double*)malloc(NDB); if (y_tilde  == NULL) { printf("malloc y_tilde failed.\n");  exit (1); }
    epsilon  = (double*)malloc(NDB); if (epsilon  == NULL) { printf("malloc epsilon failed.\n");  exit (1); }
    tmpEps   = (double*)malloc(NDB); if (tmpEps   == NULL) { printf("malloc tmpEps failed.\n");   exit (1); }
    deltaEps = (double*)malloc(NDB); if (deltaEps == NULL) { printf("malloc deltaEps failed.\n"); exit (1); }
    dEpsSum  = (double*)malloc(NDB); if (dEpsSum  == NULL) { printf("malloc dEpsSum failed.\n");  exit (1); }
    deltaSum = (double*)malloc(NDB); if (deltaSum == NULL) { printf("malloc deltaSum failed.\n"); exit (1); }

    priorPi[0] = 0.5;
    cVa[0]     = 0.0;
    cVaI[0]    = 0.0;
    muk[0]     = 0.0;
    mu         = 0.0;

    for (int i=0; i<N; ++i) {
        y_tilde[i] = 0.0;
        dEpsSum[i] = 0.0;
    }

    cVa.segment(1,km1)     = cva;
    cVaI.segment(1,km1)    = cVa.segment(1,km1).cwiseInverse();
    priorPi.segment(1,km1) = priorPi[0] * cVa.segment(1,km1).array() / cVa.segment(1,km1).sum();
    sigmaG                 = dist.beta_rng(1.0, 1.0);
    pi                     = priorPi;
    beta.setZero();
    components.setZero();

    double y_mean = 0.0;
    for (int i=0; i<N; ++i) {
        y[i]    = (double)data.y(i);
        y_mean += y[i];
    }

    y_mean /= N;
    //printf("rank %d: y_mean = %20.15f\n", rank, y_mean);

    for (int i=0; i<N; ++i) {
        y[i] -= y_mean;
    }

    double y_sqn = 0.0d;
    for (int i=0; i<N; ++i) {
        y_sqn += y[i] * y[i];
    }
    y_sqn = sqrt(y_sqn / dNm1);

    sigmaE = 0.0d;
    for (int i=0; i<N; ++i) {
        y[i]       /= y_sqn;
        //epsilon[i]  = y[i] - mu;
        epsilon[i]  = y[i];
        sigmaE     += epsilon[i] * epsilon[i];
    }
    sigmaE = sigmaE / dN * 0.5d;
    //printf("sigmaE = %20.10f with epsilon = y-mu %22.15f\n", sigmaE, mu);


    double   sum_beta_squaredNorm, sum_eps, mean_eps;
    double   sigE_G, sigG_E, i_2sigE;
    double   bet, betaOld, deltaBeta, beta_squaredNorm, p, acum, e_sqn;
    size_t   markoff;
    int      marker, left;
    VectorXd logL(K);
    

    // Main iteration loop
    // -------------------
    for (int iteration=0; iteration < max_it; iteration++) {
        
        sum_eps = 0.0;
        for (int i=0; i<N; ++i) {
            epsilon[i] += mu;
            sum_eps    += epsilon[i];
        }

        mean_eps = sum_eps/dN;

        // update mu
        mu = dist.norm_rng(mean_eps, sigmaE/dN);

        // We substract again now epsilon =Y-mu-X*beta
        for (int i=0; i<N; ++i)
            epsilon[i] -= mu;

        //EO: watch out, std::shuffle is not portable, so do no expect identical
        //    results between Intel and GCC when shuffling the markers is on!!
        //------------------------------------------------------------------------
        //for (auto i = 10; i; --i)
        //    std::cout << i << " b> " << dist.rng() << std::endl;
        //EO: shuffle or not the markers (only tests)
        if (shuf_mark) {
            std::shuffle(markerI.begin(), markerI.end(), dist.rng);
            //std::random_shuffle(markerI.begin(), markerI.end());
        }
        //for (auto i = 10; i; --i)
        //    std::cout << i << " a> " << dist.rng() << std::endl;

        
        m0 = 0.0d;
        v.setZero();

        sigE_G  = sigmaE / sigmaG;
        sigG_E  = sigmaG / sigmaE;
        i_2sigE = 1.0 / (2.0 * sigmaE);

        for (int i=0; i<N; ++i)
            tmpEps[i] = epsilon[i];

        // Loop over (shuffled) markers
        for (int j = 0; j < M; j++) {

            marker  = markerI[j];
            markoff = size_t(marker) * size_t(N);
            Cx      = &ppdata[markoff];
            //printf("%d/%d/%d: Cx[0] = %20.15f / %20.15f\n", iteration, rank, marker, Cx[0], ppdata[markoff]);

            bet =  beta(marker,0);

            scaadd(y_tilde, epsilon, Cx, bet, N);

            //we compute the denominator in the variance expression to save computations
            //denom = dNm1 + sigE_G * cVaI.segment(1, km1).array();
            for (int i=1; i<=km1; ++i)
                denom(i-1) = dNm1 + sigE_G * cVaI(i);

            //we compute the dot product to save computations
            num = dotprod(Cx, y_tilde, N);

            //muk for the other components is computed according to equations
            muk.segment(1, km1) = num / denom.array();           

            //first component probabilities remain unchanged
            logL = pi.array().log();

            // Update the log likelihood for each component
            logL.segment(1,km1) = logL.segment(1, km1).array()
                - 0.5d * (sigG_E * dNm1 * cVa.segment(1,km1).array() + 1.0d).array().log() 
                + muk.segment(1,km1).array() * num * i_2sigE;

            // I use beta(1,1) because I cant be bothered in using the std::random or create my own uniform distribution, I will change it later
            p = dist.beta_rng(1.0, 1.0);

            acum = 0.d;
            if(((logL.segment(1,km1).array()-logL[0]).abs().array() > 700 ).any() ){
                acum = 0.0d;
            } else{
                acum = 1.0d / ((logL.array()-logL[0]).exp().sum());
            }

            //EO: K -> K-1 by Daniel on 20190219!            
            //-----------------------------------
            for (int k=0; k<K-1; k++) {                
                if (p <= acum) {
                    if (k==0) {
                        beta(marker,0) = 0.0d;
                    } else {
                        beta(marker,0) = dist.norm_rng(muk[k],sigmaE/denom[k-1]);
                        //printf("@@@ beta update %d/%d/%d muk[%d] = %10.8f with p=%10.8f <= acum=%10.8f\n", iteration, rank, marker, k, muk[k], p, acum);
                    }
                    v[k] += 1.0d;
                    components[marker] = k;
                    break;
                } else {
                    //if too big or too small
                    if (((logL.segment(1,km1).array()-logL[k+1]).abs().array() > 700.0d ).any() ){
                        acum += 0.0d;
                    } else{
                        acum += 1.0d / ((logL.array()-logL[k+1]).exp().sum());
                    }
                }
            }

            //epsilon = y_tilde - Cx * beta(marker,0);
            //now epsilon contains Y-mu - X*beta+ X.col(marker)*beta(marker)_old- X.col(marker)*beta(marker)_new

            betaOld   = bet;
            bet       = beta(marker,0);
            deltaBeta = betaOld - bet;
            //printf("%d/%d/%d: deltaBeta = %20.15f = %10.7f - %10.7f\n", iteration, rank, marker, deltaBeta, betaOld, bet);
            //fflush(stdout);
            //MPI_Barrier(MPI_COMM_WORLD);

            // Compute delta epsilon
            scaadd(deltaEps, Cx, deltaBeta, N);

            // Update local sum of delta epsilon
            for (int i=0; i<N; ++i)
                dEpsSum[i] += deltaEps[i];       

            if (j%sync_rate == 0 || j-1 == M) {

                // Synchronize the deltaEps
                //MPI_Allreduce(&deltaEps[0], &deltaSum[0], N, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
                MPI_Allreduce(&dEpsSum[0], &deltaSum[0], N, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);

                for (int i=0; i<N; ++i)
                    epsilon[i] = tmpEps[i] + deltaSum[i];

                // Reset local sum of delta epsilon
                for (int i=0; i<N; ++i)
                    dEpsSum[i] = 0.0;

                // Store epsilon state at last synchronization
                for (int i=0; i<N; ++i)
                    tmpEps[i] = epsilon[i];
            }
        }

        beta_squaredNorm = beta.squaredNorm();

        // Transfer global to local
        // ------------------------
        MPI_Allreduce(&beta_squaredNorm, &sum_beta_squaredNorm, 1,        MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
        MPI_Allreduce(v.data(),          sum_v.data(),          v.size(), MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
        v                = sum_v;
        m0               = double(Mtot) - v[0];
        beta_squaredNorm = sum_beta_squaredNorm;

        sigmaG  = dist.inv_scaled_chisq_rng(v0G+m0, (beta_squaredNorm * m0 + v0G*s02G) /(v0G+m0));

        printf("it %3d, rank %d: sigmaG(%15.10f, %15.10f) = %15.10f, betasq=%15.10f, m0=%10.1f\n", iteration, rank, v0G+m0,(beta_squaredNorm * m0 + v0G*s02G) /(v0G+m0), sigmaG, beta_squaredNorm, m0);
        fflush(stdout);

        e_sqn = 0.0d;
        for (int i=0; i<N; ++i) {
            e_sqn += epsilon[i] * epsilon[i];
        }
 
        sigmaE  = dist.inv_scaled_chisq_rng(v0E+N,(e_sqn + v0E*s02E) /(v0E+N));
        
        //cout<< "inv scaled parameters "<< v0G+m0 << "__"<< (beta.squaredNorm()*m0+v0G*s02G)/(v0G+m0) << endl;
        //printf("inv scaled parameters %20.15f __ %20.15f\n", v0G+m0, (beta.squaredNorm()*m0+v0G*s02G)/(v0G+m0));
        //sigmaE = dist.inv_scaled_chisq_rng(v0E+N,((epsilon).squaredNorm()+v0E*s02E)/(v0E+N));
        //printf("sigmaG = %20.15f, sigmaE = %20.15f, e_sqn = %20.15f\n", sigmaG, sigmaE, e_sqn);
        //printf("it %6d, rank %3d: epsilon[0] = %15.10f, y[0] = %15.10f, m0=%10.1f,  sigE=%15.10f,  sigG=%15.10f [%6d / %6d]\n", iteration, rank, epsilon[0], y[0], m0, sigmaE, sigmaG, markerI[0], markerI[M-1]);

        pi = dist.dirichilet_rng(v.array() + 1.0);

        // Write to output file
        //Lineout lineout;
        //lineout.sigmaE    = sigmaE;
        //lineout.sigmaG    = sigmaG;
        //lineout.iteration = iteration;
        //lineout.rank      = rank;
        //offset = size_t(iteration) * size_t(nranks) + size_t(rank) * sizeof(lineout);
        //result = MPI_File_write_at_all(outfh, offset, &lineout, 1, typeout, &status);

        left = snprintf(buff, LENBUF, "%3d, %6d, %15.10f, %15.10f, %15.10f\n", rank, iteration, sigmaE, sigmaG, sigmaG/(sigmaE+sigmaG));
        //printf("letf = %d\n", left);
        offset = (size_t(iteration) * size_t(nranks) + size_t(rank)) * strlen(buff);
        result = MPI_File_write_at_all(outfh, offset, &buff, strlen(buff), MPI_CHAR, &status);
        if (result != MPI_SUCCESS) 
            sample_error(result, "MPI_File_write_at_all");
    }

    MPI_File_close(&outfh);
    printf("rank %d finished writing to csv output file %-100s\n", rank, outfp.c_str());
    
    //MPI_Type_free(&typeout);

    free(y);
    free(y_tilde);
    free(epsilon);
    free(tmpEps);
    free(deltaEps);
    free(dEpsSum);
    free(deltaSum);
    free(rawdata);
    free(ppdata);

    // Finalize the MPI environment. No more MPI calls can be made after this
    MPI_Finalize();

    const auto et3 = std::chrono::high_resolution_clock::now();
    const auto dt3 = et3 - st3;
    const auto du3 = std::chrono::duration_cast<std::chrono::milliseconds>(dt3).count();
    std::cout << "rank " << rank << ", time to process the data: " << du3 / double(1000.0) << " s." << std::endl;

    return 0;
}

#endif


int BayesRRm::runGibbs()
{
    const unsigned int M(data.numSnps);
    const unsigned int N(data.numInds);
    const double NM1 = double(N - 1);
    const int K(int(cva.size()) + 1);
    const int km1 = K - 1;

    //initialize variables with init member function
    init(K, M, N);

    //specify how to write samples
    SampleWriter writer;
    writer.setFileName(outputFile);
    writer.setMarkerCount(M);
    writer.setIndividualCount(N);
    writer.open();

    // Sampler variables
    VectorXd sample(2*M+4+N); // variable containg a sample of all variables in the model, M marker effects, M component assigned to markers, sigmaE, sigmaG, mu, iteration number and Explained variance
    std::vector<unsigned int> markerI(M);
    std::iota(markerI.begin(), markerI.end(), 0);

    std::cout << "Running Gibbs sampling" << endl;
    const auto t1 = std::chrono::high_resolution_clock::now();


    for (unsigned int iteration = 0; iteration < max_iterations; iteration++) {
        // Output progress
        const auto iterStart = std::chrono::high_resolution_clock::now();
        if (iteration > 0 && iteration % unsigned(std::ceil(max_iterations / 10)) == 0)
            std::cout << "iteration: " << iteration << std::endl;

        epsilon = epsilon.array() + mu;//  we substract previous value
        mu = dist.norm_rng(epsilon.sum() / (double)N, sigmaE / (double)N); //update mu
        epsilon = epsilon.array() - mu;// we substract again now epsilon =Y-mu-X*beta

        std::random_shuffle(markerI.begin(), markerI.end());

        m0 = 0;
        v.setZero();

        // This for should not be parallelized, resulting chain would not be ergodic, still, some times it may converge to the correct solution
        for (unsigned int j = 0; j < M; j++) {

            double acum = 0.0;
            const auto marker = markerI[j];
            double beta_old=beta(marker);

            //read data for column with member function getSnpData
            Cx = getSnpData(marker);

            // residual update only if marker was previously included in model
            if(components(marker)!=0){
                y_tilde=epsilon+beta_old*Cx;
            }
            else{
                y_tilde=epsilon;
            }
            // muk for the zeroth component=0
            muk[0] = 0.0;

            // We compute the denominator in the variance expression to save computations
            const double sigmaEOverSigmaG = sigmaE / sigmaG;
            denom = NM1 + sigmaEOverSigmaG * cVaI.segment(1, km1).array();

            // We compute the dot product to save computations
            const double num = (Cx.dot(y_tilde));

            // muk for the other components is computed according to equations
            muk.segment(1, km1) = num / denom.array();

            // Update the log likelihood for each variance component
            const double logLScale = sigmaG / sigmaE * NM1;
            logL = pi.array().log(); // First component probabilities remain unchanged
            logL.segment(1, km1) = logL.segment(1, km1).array()
                    		        - 0.5 * ((logLScale * cVa.segment(1, km1).array() + 1).array().log())
                    		        + 0.5 * (muk.segment(1, km1).array() * num) / sigmaE;

            double p(dist.unif_rng());
            if (((logL.segment(1, km1).array() - logL[0]).abs().array() > 700).any()) { // this quantity will be exponentiated and go to huge values and we set to infinity, thus the reciprocal is set to zero
                acum = 0;
            } else {
                acum = 1.0 / ((logL.array() - logL[0]).exp().sum());
            }

            for (int k = 0; k < K; k++) {
                if (p <= acum) {
                    //if zeroth component
                    if (k == 0) {
                        beta(marker) = 0;
                    } else {
                        beta(marker) = dist.norm_rng(muk[k], sigmaE/denom[k-1]);
                    }
                    v[k] += 1.0;
                    components[marker] = k;
                    break;
                } else {
                    if (((logL.segment(1, km1).array() - logL[k+1]).abs().array() > 700).any()) {//again controlling for large values for the exponent
                        acum += 0;
                    } else {
                        acum += 1.0 / ((logL.array() - logL[k+1]).exp().sum());
                    }
                }
            }
            betasqn+=beta(marker)*beta(marker)-beta_old*beta_old;
            // residual update only if updated marker is included in model
            // Now epsilon contains Y-mu - X*beta + X.col(marker) * beta(marker)_old - X.col(marker) * beta(marker)_new
            if(components(marker)!=0){
                epsilon=y_tilde-beta(marker)*Cx;
            }
            else{
                epsilon=y_tilde;
            }
        }
        //set no. of markers included in the model
        m0 = int(M) - int(v[0]);
        //sample sigmaG from inverse gamma
        sigmaG = dist.inv_scaled_chisq_rng(v0G + m0, (betasqn * m0 + v0G * s02G) / (v0G + m0));

        const double epsilonSqNorm=epsilon.squaredNorm();
        //sample residual variance sigmaE from inverse gamma
        sigmaE = dist.inv_scaled_chisq_rng(v0E + N, (epsilonSqNorm + v0E * s02E) / (v0E + N));
        //sample hyperparameter pi from dirichlet
        pi = dist.dirichilet_rng(v.array() + 1.0);

        if (showDebug)
            printDebugInfo();

        //write samples
        if (iteration >= burn_in && iteration % thinning == 0) {
            sample << iteration, mu, beta, sigmaE, sigmaG, components, epsilon;
            writer.write(sample);
        }

        //output time taken for each iteration
        const auto endTime = std::chrono::high_resolution_clock::now();
        const auto dif = endTime - iterStart;
        const auto iterationDuration = std::chrono::duration_cast<std::chrono::milliseconds>(dif).count();
        std::cout << iterationDuration / double(1000.0) << "s" << std::endl;

        //end of iteration
    }
    //show info on total time spent
    const auto t2 = std::chrono::high_resolution_clock::now();
    const auto duration = std::chrono::duration_cast<std::chrono::seconds>(t2 - t1).count();
    std::cout << "duration: " << duration << "s" << std::endl;

    return 0;
}

VectorXd BayesRRm::getSnpData(unsigned int marker) const
{
    if (!usePreprocessedData) {
        //read column from RAM loaded genotype matrix.
        return data.Z.col(marker).cast<double>();
    } else {
        //read column from preprocessed and memory mapped genotype matrix file.
        return data.mappedZ.col(marker).cast<double>();
    }
}

void BayesRRm::printDebugInfo() const
{
    const unsigned int N(data.numInds);
    cout << "inv scaled parameters " << v0G + m0 << "__" << (beta.squaredNorm() * m0 + v0G * s02G) / (v0G + m0);
    cout << "num components: " << opt.S.size();
    cout << "\nMixture components: " << cva[0] << " " << cva[1] << " " << cva[2] << "\n";
    cout << "sigmaG: " << sigmaG << "\n";
    cout << "y mean: " << y.mean() << "\n";
    cout << "y sd: " << sqrt(y.squaredNorm() / (double(N - 1))) << "\n";
    // cout << "x mean " << Cx.mean() << "\n";
    //   cout << "x sd " << sqrt(Cx.squaredNorm() / (double(N - 1))) << "\n";
}
