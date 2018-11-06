#!/bin/bash

# E. Orliac
# 04 Sep 2018
#
# Run reference solution

dataset=uk10k_chr1_1mb

out=$dataset.mmap

echo Running mmap solution on dataset $dataset. Output in $out

./brr --bayes gbayes --bfile ../test/data/$dataset --group /Users/mpatxotb/Desktop/read.groupsCpp.C/groupfile.txt --mS '0.0001,0.001,0.01;0.0001,0.001,0.01;0.0001,0.001,0.01' --pheno ../test/data/test.phen --chain-length 10 --burn-in 5 --thin 2 --mcmc-samples ./bayesOutput.csv --S 0.01,0.001,0.0001 

out=$dataset.pp
echo Running preprocessing solution on dataset $dataset. Output in $out

#./brr --preprocess  --bfile ../test/data/$dataset --pheno ../test/data/test.phen > $out
#./brr --ppbayes  bayes --bfile ../test/data/$dataset --pheno ../test/data/test.phen>>$out
echo done
