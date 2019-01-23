/*
 * Brrcmdtest.cpp
 *
 *  Created on: Jan 15, 2019
 *      Author: daniel
 */

#include<gtest/gtest.h>
#include "data.hpp"
#include<Eigen/Eigen>
#include<distributions_boost.hpp>
#include "gadgets.hpp"
#include "options.hpp"
#include "BayesRRm.h"
#include "BayesRRmz.h"
#include "limitsequencegraph.cpp"

TEST(Brrcmdtest,testbedfilereading){

	string bedFile="test/data/uk10k_chr1_1mb";
	string phenotypeFile ="test/data/test.phen";
	bool success=1;

	Data data;
	data.readFamFile("test/data/uk10k_chr1_1mb.fam");
    data.readPhenotypeFile(phenotypeFile, 1);
	data.readBimFile(bedFile + ".bim");

	ASSERT_EQ(data.numSnps,6717);
	ASSERT_EQ(data.numInds,3642);
	ASSERT_TRUE(success);
}

TEST(Brrcmdtest,testBrrm_preprocess){
	string bedFile="test/data/uk10k_chr1_1mb";
		string phenotypeFile ="test/data/test.phen";
		bool success=1;
        std::vector<float> S={0.001,0.01};
		Data data;
		data.readFamFile("test/data/uk10k_chr1_1mb.fam");
	    data.readPhenotypeFile(phenotypeFile, 1);
		data.readBimFile(bedFile + ".bim");
		data.numKeptInds=data.numInds;//TODO get rid of these variables
		data.numIncdSnps=data.numSnps;
		Options opt;
		opt.bedFile=bedFile;
		opt.phenotypeFile=phenotypeFile;
		opt.chainLength=10;
		opt.burnin=5;
		opt.mcmcSampleFile=("test/data/testBrrm.csv");
		opt.thin=1;

	    data.preprocessBedFile("test/data/uk10k_chr1_1mb.bed",
		                                   "test/data/uk10k_chr1_1mb.ppbed",
		                                   "test/data/uk10k_chr1_1mb.ppbedindex",
		                                   0);
	    data.mapPreprocessBedFile("test/data/uk10k_chr1_1mb.ppbed");

	    const Eigen::VectorXd &Cx0 = data.mappedZ.col(0);
        const Eigen::VectorXd &Cxend= data.mappedZ.col(data.numIncdSnps-1);

	    data.unmapPreprocessedBedFile();
	    ASSERT_LT((double)Cx0.mean(),1e-8)<<"BED files correctly preprocessed and mapped";
	    ASSERT_LT((double)Cxend.mean(),1e-8);
	    ASSERT_LT(abs((double)Cx0.squaredNorm()-(data.numKeptInds-1)),1e-8);
	    ASSERT_LT(abs((double)Cxend.squaredNorm()-(data.numKeptInds-1)),1e-8);


}

TEST(Brrcmdtest,testBrrm_compress){
	string bedFile="test/data/uk10k_chr1_1mb";
		string phenotypeFile ="test/data/test.phen";
		bool success=1;
        std::vector<float> S={0.001,0.01};
		Data data;
		data.readFamFile("test/data/uk10k_chr1_1mb.fam");
	    data.readPhenotypeFile(phenotypeFile, 1);
		data.readBimFile(bedFile + ".bim");
		data.numKeptInds=data.numInds;//TODO get rid of these variables
		data.numIncdSnps=data.numSnps;
		Options opt;
		opt.bedFile=bedFile;
		opt.phenotypeFile=phenotypeFile;
		opt.chainLength=10;
		opt.burnin=5;
		opt.mcmcSampleFile="test/data/testBrrm.csv";
		opt.thin=1;
		opt.S=S;
		opt.seed=1;
		data.preprocessBedFile("test/data/uk10k_chr1_1mb.bed",
				                                   "test/data/uk10k_chr1_1mb.ppbed",
				                                   "test/data/uk10k_chr1_1mb.ppbedindex",
				                                   1);
		data.mapCompressedPreprocessBedFile("test/data/uk10k_chr1_1mb.ppbed",
				"test/data/uk10k_chr1_1mb.ppbedindex");
        data.unmapCompressedPreprocessedBedFile();

	    ASSERT_TRUE(success)<<"BED files correctly preprocessed compressed and mapped";

}


TEST(Brrcmdtest,testBrrm_sparse){
	string bedFile="test/data/uk10k_chr1_1mb";
			string phenotypeFile ="test/data/test.phen";
			bool success=1;
	        std::vector<float> S={0.001,0.01};
			Data data;
			data.readFamFile("test/data/uk10k_chr1_1mb.fam");
		    data.readPhenotypeFile(phenotypeFile, 1);
			data.readBimFile(bedFile + ".bim");
			data.numKeptInds=data.numInds;//TODO get rid of these variables
			data.numIncdSnps=data.numSnps;
			data.scanBedFile("test/data/uk10k_chr1_1mb.bed");

	        ASSERT_TRUE(success);

}
