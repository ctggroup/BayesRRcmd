# script to check bayesG output
# 24.05.19 Marion Patxot

#libraries
library(data.table)
library(dplyr)

#read output file
out<-fread('output.csv',header=TRUE)
out %>% select('sigmaG[1]','sigmaG[2]')

#read true variance
var<-fread('true.var.txt',header=F)
