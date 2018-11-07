#!/bin/bash

# E. Orliac
# 04 Sep 2018
#
# Run reference solution

dataset=uk10k_chr1_1mb

out=$dataset.mmap

echo Running mmap solution on dataset $dataset. Output in $out

./brr --bayes bayesMmap --bfile ../test/data/$dataset --group /Users/mpatxotb/Desktop/read.groupsCpp.C/groupfile.txt --pheno ../test/data/test.phen --chain-length 1000 --burn-in 500 --thin 2 --mcmc-samples ./bayesOutput_bayesm.csv --S 0.01,0.001,0.0001 


#./brr --bayes gbayes --bfile ../test/data/$dataset --group /Users/mpatxotb/Desktop/read.groupsCpp.C/groupfile.txt --mS '0.0001,0.001,0.01;0.0002,0.002,0.02;0.0005,0.005,0.05' --pheno ../test/data/test.phen --chain-length 1000 --burn-in 500 --thin 2 --mcmc-samples ./bayesOutput_bayesg.csv


out=$dataset.pp
echo Running preprocessing solution on dataset $dataset. Output in $out

#./brr --preprocess  --bfile ../test/data/$dataset --pheno ../test/data/test.phen > $out
#./brr --ppbayes  bayes --bfile ../test/data/$dataset --pheno ../test/data/test.phen>>$out
echo done
