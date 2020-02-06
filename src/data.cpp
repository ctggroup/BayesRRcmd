#include "data.hpp"
#include <Eigen/Eigen>
#include <sys/mman.h>
#include <sys/stat.h>
#include <fcntl.h>
#include <iterator>
#include <filesystem>
#include "compression.h"
#include "markersubset.h"
#include "options.hpp"

namespace fs = std::filesystem;

#define handle_error(msg)                               \
    do { perror(msg); exit(EXIT_FAILURE); } while (0)

Data::Data()
    : ppBedFd(-1)
    , ppBedMap(nullptr)
    , mappedZ(nullptr, 1, 1)
    , ppbedIndex()
    , X(nullptr, 0, 0)
{
}

Data::~Data()
{
    unmapPreprocessedBedFile();
}

void Data::mapPreprocessBedFile(const Options &options)
{
    // Read the MarkerSubset this preprocessed file contains
    MarkerSubset subset = {0, 0};
    const auto ppSubsetFile = ppSubsetFileForType(options);
    if (fs::exists(ppSubsetFile)) {
        ifstream subsetStream(ppSubsetFile, std::ifstream::binary);
        if (!subsetStream)
            throw("Error: Failed to open preprocessed subset file");

        subsetStream >> subset;
        if (subsetStream.fail())
            throw("Error: Failed to read preprocessed subset file");
    }
    setMarkerSubset(subset);

    // Load the index to the compressed preprocessed bed file
    ppbedIndex.resize(numSnps);
    const auto indexFile = ppIndexFileForType(options);
    ifstream indexStream(indexFile, std::ifstream::binary);
    if (!indexStream)
        throw("Error: Failed to open compressed preprocessed bed file index");

    indexStream.read(reinterpret_cast<char *>(ppbedIndex.data()),
                     numSnps * 3 * sizeof(unsigned long));

    // Open and mmap the preprocessed bed file
    const auto preprocessedBedFile = ppFileForType(options);
    ppBedFd = open(preprocessedBedFile.c_str(), O_RDONLY);
    if (ppBedFd == -1)
        throw("Error: Failed to open preprocessed bed file [" + preprocessedBedFile + "]");

    const size_t ppBedSize = preprocessedDataMapSize();
    ppBedMap = reinterpret_cast<double *>(mmap(nullptr, ppBedSize, PROT_READ, MAP_SHARED, ppBedFd, 0));
    if (ppBedMap == MAP_FAILED)
        throw("Error: Failed to mmap preprocessed bed file");
}

void Data::unmapPreprocessedBedFile()
{
    if (ppBedMap == nullptr)
        return;

    // Unmap the data from the Eigen accessors
    new (&mappedZ) Map<MatrixXd>(nullptr, 1, 1);

    // Unmap the preprocessed data
    const size_t ppBedSize = preprocessedDataMapSize();
    munmap(ppBedMap, ppBedSize);
    ppBedMap = nullptr;

    close(ppBedFd);

    ppbedIndex.clear();
}

void Data::readFamFile(const string &famFile){
    // ignore phenotype column
    ifstream in(famFile.c_str());
    if (!in) throw ("Error: can not open the file [" + famFile + "] to read.");

    cout << "Reading PLINK FAM file from [" + famFile + "]." << endl;

    indInfoVec.clear();
    indInfoMap.clear();
    string fid, pid, dad, mom, sex, phen;
    unsigned idx = 0;
    while (in >> fid >> pid >> dad >> mom >> sex >> phen) {
        auto ind = std::make_shared<IndInfo>(idx++, fid, pid, dad, mom, atoi(sex.c_str()));
        indInfoVec.push_back(ind);
        if (indInfoMap.insert(pair<string, IndInfoPtr>(ind->catID, ind)).second == false) {
            throw ("Error: Duplicate individual ID found: \"" + fid + "\t" + pid + "\".");
        }
    }
    in.close();
    totalIndividuals = static_cast<unsigned int>(indInfoVec.size());
    activeIndividuals = totalIndividuals;

    cout << totalIndividuals << " individuals to be included from [" + famFile + "]." << endl;
}

void Data::readBimFile(const string &bimFile) {
    // Read bim file: recombination rate is defined between SNP i and SNP i-1
    ifstream in(bimFile.c_str());
    if (!in) throw ("Error: can not open the file [" + bimFile + "] to read.");

    cout << "Reading PLINK BIM file from [" + bimFile + "]." << endl;
    snpInfoVec.clear();
    snpInfoMap.clear();
    string id, allele1, allele2;
    unsigned chr, physPos;
    float genPos;
    unsigned idx = 0;
    while (in >> chr >> id >> genPos >> physPos >> allele1 >> allele2) {
        auto snp = std::make_shared<SnpInfo>(idx++, id, allele1, allele2, chr, genPos, physPos);
        snpInfoVec.push_back(snp);
        if (snpInfoMap.insert(pair<string, SnpInfoPtr>(id, snp)).second == false) {
            throw ("Error: Duplicate SNP ID found: \"" + id + "\".");
        }
    }
    in.close();
    numSnps = (unsigned) snpInfoVec.size();
    G = vector<int>(numSnps, 0);

    cout << numSnps << " SNPs to be included from [" + bimFile + "]." << endl;
}


void Data::readBedFile_noMPI(const string &bedFile){
    unsigned i = 0, j = 0, k = 0;

    if (numSnps == 0) throw ("Error: No SNP is retained for analysis.");
    if (activeIndividuals == 0) throw ("Error: No individual is retained for analysis.");

    Z.resize(activeIndividuals, numSnps);
    ZPZdiag.resize(numSnps);
    snp2pq.resize(numSnps);

    // Read bed file
    char ch[1];
    bitset<8> b;
    unsigned allele1=0, allele2=0;
    ifstream BIT(bedFile.c_str(), ios::binary);
    if (!BIT) throw ("Error: can not open the file [" + bedFile + "] to read.");
    cout << "Reading PLINK BED file from [" + bedFile + "] in SNP-major format ..." << endl;
    for (i = 0; i < 3; i++) BIT.read(ch, 1); // skip the first three bytes
    SnpInfoPtr snpInfo = nullptr;
    unsigned snp = 0, ind = 0;
    unsigned nmiss = 0;
    float mean = 0.0;

    for (j = 0, snp = 0; j < numSnps; j++) { // Read genotype in SNP-major mode, 00: homozygote AA; 11: homozygote BB; 10: hetezygote; 01: missing
        snpInfo = snpInfoVec[j];
        mean = 0.0;
        nmiss = 0;
        if (!snpInfo->included) {
            for (i = 0; i < totalIndividuals; i += 4) BIT.read(ch, 1);
            continue;
        }
        for (i = 0, ind = 0; i < totalIndividuals && ind < activeIndividuals;) {
            BIT.read(ch, 1);
            if (!BIT) throw ("Error: problem with the BED file ... has the FAM/BIM file been changed?");
            b = ch[0];
            k = 0;
            while (k < 7 && i < totalIndividuals && ind < activeIndividuals) {
                if (!indInfoVec[i]->kept) k += 2;
                else {
                    allele1 = (!b[k++]);
                    allele2 = (!b[k++]);
                    if (allele1 == 0 && allele2 == 1) {  // missing genotype
                        Z(ind++, snp) = -9;
                        ++nmiss;
                    } else {
                        mean += Z(ind++, snp) = allele1 + allele2;
                    }
                }
                i++;
            }
        }
        // fill missing values with the mean
        float sum = mean;
        mean /= float(activeIndividuals-nmiss);


        if (nmiss) {
            for (i=0; i<activeIndividuals; ++i) {
                if (Z(i,snp) == -9) Z(i,snp) = mean;
            }
        }

        // compute allele frequency
        snpInfo->af = 0.5f*mean;
        snp2pq[snp] = 2.0f*snpInfo->af*(1.0f-snpInfo->af);

        // standardize genotypes


        //	Z.col(j).array() -= mean;

        //	float sqn = Z.col(j).squaredNorm();
        float sqn = Z.col(j).squaredNorm() - activeIndividuals * mean * mean;
        Z.col(j).array() -= mean;
        float std_ = 1.f / (sqrt(sqn / float(activeIndividuals)));

        Z.col(j).array() *= std_;

        ZPZdiag[j] = sqn;

        if (++snp == numSnps) break;
    }
    BIT.clear();
    BIT.close();
    // standardize genotypes
    for (i=0; i<numSnps; ++i) {
        Z.col(i).array() -= Z.col(i).mean();
        ZPZdiag[i] = Z.col(i).squaredNorm();
    }
    cout << "Genotype data for " << activeIndividuals << " individuals and " << numSnps << " SNPs are included from [" + bedFile + "]." << endl;
}

void Data::readBedFile_noMPI_unstandardised(const string &bedFile){
    unsigned i = 0, j = 0, k = 0;

    if (numSnps == 0) throw ("Error: No SNP is retained for analysis.");
    if (activeIndividuals == 0) throw ("Error: No individual is retained for analysis.");

    Zones.resize(numSnps);
    Ztwos.resize(numSnps);
    means.resize(numSnps);
    sds.resize(numSnps);
    mean_sd_ratio.resize(numSnps);

    //	Z.resize(numInds, numSnps);
    //	ZPZdiag.resize(numSnps);
    snp2pq.resize(numSnps);

    // Read bed file
    char ch[1];
    bitset<8> b;
    unsigned allele1=0, allele2=0;
    ifstream BIT(bedFile.c_str(), ios::binary);
    if (!BIT) throw ("Error: can not open the file [" + bedFile + "] to read.");
    cout << "Reading PLINK BED file from [" + bedFile + "] in SNP-major format ..." << endl;
    for (i = 0; i < 3; i++) BIT.read(ch, 1); // skip the first three bytes
    SnpInfoPtr snpInfo = nullptr;
    unsigned snp = 0, ind = 0;
    unsigned nmiss = 0;
    float mean = 0.0;
    float sqn = 0.0;

    for (j = 0, snp = 0; j < numSnps; j++) { // Read genotype in SNP-major mode, 00: homozygote AA; 11: homozygote BB; 10: hetezygote; 01: missing
        snpInfo = snpInfoVec[j];
        mean = 0.0;
        sqn = 0.0;
        nmiss = 0;
        if (!snpInfo->included) {
            for (i = 0; i < totalIndividuals; i += 4) BIT.read(ch, 1);
            continue;
        }
        for (i = 0, ind = 0; i < totalIndividuals && ind < activeIndividuals;) {
            BIT.read(ch, 1);
            if (!BIT) throw ("Error: problem with the BED file ... has the FAM/BIM file been changed?");
            b = ch[0];
            k = 0;
            while (k < 7 && i < totalIndividuals && ind < activeIndividuals) {
                if (!indInfoVec[i]->kept) k += 2;
                else {
                    allele1 = (!b[k++]);
                    allele2 = (!b[k++]);
                    //Assume no missing genotypes
                    //if (allele1 == 0 && allele2 == 1) {  // missing genotype
                    //	Z(ind++, snp) = -9;
                    //++nmiss;
                    //	} else {
                    int all_sum = allele1 + allele2;
                    if(all_sum == 1 ){
                        Zones[j].push_back(ind++); //Save the index of the individual to the vector of ones
                    }else if(all_sum == 2){
                        Ztwos[j].push_back(ind++);
                    }else{
                        ind++;
                    }

                    if (allele1 == 0 && allele2 == 1) {  // missing genotype
                        cout << "MISSING " << endl;
                    }

                    mean += all_sum;
                    sqn += all_sum*all_sum;
                    //	}
                }
                i++;
            }
        }
        // fill missing values with the mean
        mean /= float(activeIndividuals-nmiss);

        //Assume non-missingness

        /*	if (nmiss) {
            for (i=0; i<numInds; ++i) {
                if (Z(i,snp) == -9) Z(i,snp) = mean;
            }
        }
         */
        // compute allele frequency
        snpInfo->af = 0.5f*mean;
        snp2pq[snp] = 2.0f*snpInfo->af*(1.0f-snpInfo->af);

        // standardize genotypes

        sqn -= activeIndividuals * mean * mean;
        means(j) = mean;
        sds(j) = sqrt(sqn / float(activeIndividuals));

        mean_sd_ratio(j) = means(j)/sds(j);

        if (++snp == numSnps) break;
    }
    BIT.clear();
    BIT.close();

    cout << "Genotype data for " << activeIndividuals << " individuals and " << numSnps << " SNPs are included from [" + bedFile + "]." << endl;
}


void Data::readPhenotypeFile(const string &phenFile) {
    // NA: missing phenotype
    ifstream in(phenFile.c_str());
    if (!in) throw ("Error: can not open the phenotype file [" + phenFile + "] to read.");

    cout << "Reading phenotypes from [" + phenFile + "]." << endl;
    map<string, IndInfoPtr>::iterator it, end=indInfoMap.end();
    IndInfoPtr ind = nullptr;
    Gadget::Tokenizer colData;
    string inputStr;
    string sep(" \t");
    string id;
    unsigned line=0;
    //correct loop to go through numInds
    y = VectorXf::Zero(activeIndividuals);

    while (getline(in,inputStr)) {
        colData.getTokens(inputStr, sep);
        id = colData[0] + ":" + colData[1];
        it = indInfoMap.find(id);
        // First one corresponded to mphen variable (1+1)
        if (it != end && colData[1+1] != "NA") {
            ind = it->second;
            ind->phenotype = atof(colData[1+1].c_str());
            //fill in phenotype y vector
            if (ind->kept) {
                y[line] = ind->phenotype;
                ++line;
            }
        }
    }

    in.close();
}

//Copied from
//https://gist.github.com/infusion/43bd2aa421790d5b4582
void Data::readFixedEffects(const string &filename, unsigned int cols) {

    // Clear the existing map;
    numFixedEffects = 0;
    new (&X) Map<MatrixXf>(nullptr, 0, 0);
    m_fixedEffectsData = std::make_shared<std::vector<float>>();
    m_removedIndividualIndexes = std::make_shared<std::vector<int>>();

    std::ifstream in(filename);
    if (!in.is_open()) {
        cout << "Failed to open fixed effects file: " << filename << std::endl;
        return;
    }

    cout << "Reading " << cols << " fixed effects from: " << filename << endl;

    m_fixedEffectsData->reserve(totalIndividuals * static_cast<unsigned int>(cols));
    auto rowStartItr = m_fixedEffectsData->end();

    std::string line;
    unsigned int row = 0;
    int individual = 0;

    auto isNA = [](const char* str) {
        static const char* na = "NA";
        return strncmp(str, na, 2) == 0;
    };

    while (std::getline(in, line) && row < totalIndividuals) {

        auto emplaceItr = rowStartItr;
        const char *ptr = line.c_str();
        size_t len = line.length();
        bool naFound = false;

        char *start = const_cast<char *>(ptr);
        for (size_t i = 0; i < len; i++) {
            if (ptr[i] == ',') {
                if (isNA(start)) {
                    naFound = true;
                    break;
                }
                m_fixedEffectsData->emplace(emplaceItr++, atof(start));
                start = const_cast<char *>(ptr) + i + 1;
            }
        }

        if (naFound || isNA(start)) {
            m_fixedEffectsData->erase(rowStartItr, m_fixedEffectsData->end());
            m_removedIndividualIndexes->emplace_back(individual++);
            continue;
        }

        m_fixedEffectsData->emplace(emplaceItr++, atof(start));
        rowStartItr = m_fixedEffectsData->end();
        ++row;
        ++individual;
    }

    if (row >= totalIndividuals)
        cerr << "Warning: fixed effects data contains more individuals than the data set!" << endl;

    in.close();

    new (&X) Map<MatrixXf>(m_fixedEffectsData->data(), row, cols);
    numFixedEffects = cols;

    cout << "Read " << numFixedEffects << " fixed effects and found "
         << m_removedIndividualIndexes->size() << " NA individuals." << endl;
}

void Data::updateNaIndividuals()
{
    if (!m_removedIndividualIndexes || m_removedIndividualIndexes->empty())
        return;

    cout << "Removing " << m_removedIndividualIndexes->size() << " NA individuals from the data set." << endl;

    std::for_each(indInfoVec.begin(), indInfoVec.end(), [&naIndices = m_removedIndividualIndexes](const IndInfoPtr& individual){
        const auto findItr = std::find(naIndices->cbegin(), naIndices->cend(), individual->index);
        individual->kept = findItr == naIndices->cend();
    });

    activeIndividuals -= m_removedIndividualIndexes->size();
}


void Data::readFailureFile(const string &failureFile){
    ifstream input(failureFile);
    vector<double> tmp;
    tmp.reserve(totalIndividuals);

    int col;
    if(!input.is_open()){
        cout << "Error opening the file" << endl;
        return;
    }

    unsigned int index = 0;
    while(true){
        input >> col ;
        if(input.eof()) break;

        if (!indInfoVec[index++]->kept) {
            continue;
        }

        tmp.push_back(col);
    }
    input.close();
    fail = Eigen::VectorXd::Map(tmp.data(), tmp.size());
}

size_t Data::preprocessedDataMapSize() const
{
    // Calculate the expected file sizes - return size_t so that we don't
    // overflow the unsigned int's that we would otherwise get as intermediate
    // variables!
    const auto lastIndex = ppbedIndex[m_markerSubset.last()];
    return lastIndex.pos + lastIndex.compressedSize;
}

void Data::readGroupFile(const string& groupFile) {
    G.clear();
    numGroups = 1;

    ifstream input(groupFile);
    if(!input.is_open()){
        cout<<"Error opening the file"<< endl;
        return;
    }

    // Read the groups
    G.reserve(numSnps);

    string col1;
    int col2;
    while(true){
        input >> col1 >> col2;
        if(input.eof()) break;
        G.emplace_back(col2);
    }

    // Calculate the number of groups
    auto temp{G};
    std::sort(temp.begin(), temp.end());
    auto last = std::unique(temp.begin(), temp.end());
    temp.erase(last, temp.end());
    numGroups = temp.size();

    cout << "Groups read from file: " << numGroups << endl;
}

void Data::preprocessCSVFile(const string&csvFile,const string &preprocessedCSVFile, const string &preprocessedCSVIndexFile, bool compress)
{
    cout << "Preprocessing csv file:" << csvFile << ", Compress data =" << (compress ? "yes" : "no") << endl;

    VectorXd snpData = VectorXd::Zero(activeIndividuals);

    std::ifstream indata;
    indata.open(csvFile);
    std::string line;
    uint rows = 0;
    uint cols =0;
    unsigned int individual = 0;
    ofstream ppCSVOutput(preprocessedCSVFile.c_str(), ios::binary);
    if (!ppCSVOutput) {
        cerr << "Error: Unable to open the preprocessed bed file [" + preprocessedCSVFile + "] for writing." << endl;
        throw("Error: Unable to open the preprocessed bed file [" + preprocessedCSVFile + "] for writing.");
    }
    ofstream ppCSVIndexOutput(preprocessedCSVIndexFile.c_str(), ios::binary);
    if (!ppCSVIndexOutput) {
        cerr << "Error: Unable to open the preprocessed bed index file [" + preprocessedCSVIndexFile + "] for writing." << endl;
        throw("Error: Unable to open the preprocessed bed index file [" + preprocessedCSVIndexFile + "] for writing.");
    }

    // How much space do we need to compress the data (if requested)
    const auto maxCompressedOutputSize = compress ? maxCompressedDataSize<double>(activeIndividuals) : 0;
    unsigned char *compressedBuffer = nullptr;
    unsigned long pos = 0;
    if (compress)
        compressedBuffer = new unsigned char[maxCompressedOutputSize];
    while (std::getline(indata, line))
    {
        std::stringstream lineStream(line);
        std::string cell;
        cols=0;
        individual = 0;


        while (std::getline(lineStream, cell, ','))
        {
            if (cell.empty()) {
                cerr << "Error, there are missing values in the file" << endl;
                throw("Error, there are missing values in the file");
            }

            if (!indInfoVec[individual++]->kept) {
                continue;
            }

            snpData[cols++]=std::stod(cell);
        }

        if (!compress)
        {
            writeUncompressedDataWithIndex(reinterpret_cast<unsigned char *>(&snpData[0]), activeIndividuals * sizeof(double), ppCSVOutput, ppCSVIndexOutput, pos);
        }
        else
        {
            compressAndWriteWithIndex(snpData, ppCSVOutput, ppCSVIndexOutput, pos, compressedBuffer, maxCompressedOutputSize);
        }
        ++rows;
    }
    if (compress)
        delete[] compressedBuffer;
    indata.clear();
    indata.close();

    cout << "csv file for" << activeIndividuals << " individuals and " << numSnps << " Variables are included from [" +  csvFile + "]." << endl;
}

bool Data::setMarkerSubset(const MarkerSubset &subset)
{
    if (!subset.isValid(numSnps))
        return false;

    if (subset.first() == 0 && subset.size() == 0) {
        m_markerSubset = MarkerSubset(0, numSnps);
    } else {
        m_markerSubset = subset;
    }

    return true;
}

bool Data::validMarkerSubset() const
{
    return m_markerSubset.isValid(numSnps);
}

std::vector<unsigned int> Data::getMarkerIndexList() const
{
    return m_markerSubset.toMarkerIndexList(numSnps);
}

//We asume the csv is well formed and individuals are columns and  markers are rows
void Data::readCSVFile( const string &csvFile)
{
    std::ifstream indata;
    indata.open(csvFile);
    if (!indata)
        throw("Error: Unable to open the CSV data file [" + csvFile + "] for reading.");
    std::string line;
    uint rows = 0;
    uint cols = 0;
    while (std::getline(indata, line))
    {
        if(rows == 0){
            std::stringstream lineStream(line);
            std::string cell;

            while (std::getline(lineStream, cell, ','))
            {
                indInfoVec.emplace_back(std::make_shared<IndInfo>(cols++));
            }
        }
        ++rows;
    }
    totalIndividuals = cols;
    activeIndividuals = totalIndividuals;
    numSnps = rows;
    G = vector<int>(numSnps, 0);
    indata.clear();
    indata.close();
    cout << totalIndividuals << " individuals to be included from [" + csvFile + "]." << endl;
    cout << numSnps << " markers to be included from [" + csvFile + "]." << endl;
}

void Data::readCSVPhenFile( const string &csvFile)
{
    std::ifstream indata;
    indata.open(csvFile);
    if (!indata)
        throw("Error: Unable to open the CSV phenotype file [" + csvFile + "] for reading.");
    std::string line;
    uint cols = 0;
    y = VectorXf::Zero(activeIndividuals);
    std::getline(indata, line);

    std::stringstream lineStream(line);
    std::string cell;

    unsigned int individual = 0;

    while (std::getline(lineStream, cell, ',') && cols<activeIndividuals)
    {
        if (cell.empty())
            throw("Error, there are missing values in the file");

        if (!indInfoVec[individual++]->kept) {
            continue;
        }

        y[cols++]= std::stof(cell);
    }

    indata.clear();
    indata.close();

    cout << cols << " phenotype measures to be included from [" + csvFile + "]." << endl;
}
