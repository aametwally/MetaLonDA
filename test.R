library(ggplot2)
library(gss)
library(plyr)
library(plyr)
library(parallel)
library(doParallel)
library(metagenomeSeq)
library(DESeq2)
library(edgeR)


setwd("/Users/ahmedmetwally/Box Sync/metalonda/Code/metalonda_github/MetaLonDA/")
source("R/Metalonda.R")
source("R/CurveFitting.R")
source("R/Permutation.R")
source("R/Visualization.R")
source("R/Normalization.R")
source("R/metalonda_test_data.R")




data(metalonda_test_data)
n.sample = 5
n.timepoints = 10
n.group = 2
Group = factor(c(rep(0, n.sample*n.timepoints), rep(1,n.sample*n.timepoints)))
Time = rep(rep(1:n.timepoints, times = n.sample), 2)
ID = factor(rep(1:(2*n.sample), each = n.timepoints))
points = seq(1, 10, length.out = 10)
output.nbinomial = metalonda(Count = metalonda_test_data[1,], Time = Time, Group = Group,
ID = ID, fit.method =  "nbinomial", n.perm = 10, points = points,
text = rownames(metalonda_test_data)[1], parall = FALSE, pvalue.threshold = 0.05, adjust.method = "BH")



output.nbinomial = metalondaAll(Count = metalonda_test_data, Time = Time, Group = Group,
ID = ID, n.perm = 9, fit.method =  "nbinomial", num.intervals = 100, 
parall = FALSE, pvalue.threshold = 0.05, adjust.method = "BH", time.unit = "hours", norm.method = "none",
prefix = "Test")
