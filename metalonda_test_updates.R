setwd("~/Dropbox/metalonda_work/updateMetalonda/MetaLonDA")
detach("package:MetaLonDA", unload=TRUE)
source("R/CurveFitting.R")
source("R/Visualization.R")
source("R/Permutation.R")
source("R/Normalization.R")
source("R/Metalonda.R")
source("R/metalonda_test_data.R")



library(ggplot2)
library(gss)
library(plyr)
library(caTools)
library(parallel)
library(doParallel)
library(zoo)



data(metalonda_test_data)
View(metalonda_test_data[,1:20])


## Create Group, Time, and ID annotation vectors
n.group = 2
n.sample = 5 
n.timepoints = 10
Group = factor(c(rep("A", n.sample*n.timepoints), rep("B",n.sample*n.timepoints)))
Time = rep(rep(1:n.timepoints, times = n.sample), 2)
ID = factor(rep(1:(2*n.sample), each = n.timepoints))

## Define the prediction timeponits 
points = seq(1, 10, length.out = 100)




Count= metalonda_test_data[5,]
Time = Time
Group = Group
ID = ID
n.perm = 5
fit.method = "nbinomial"
points = points
text = "testfeature"
parall = FALSE
pvalue.threshold = 0.05
adjust.method = "BH"
time.unit = "pppp"
ylabel = "Normalized Count"
col = c("blue", "firebrick")


output.metalonda.f5 = metalonda(Count = metalonda_test_data[5,], Time = Time, Group = Group,
                                ID = ID, n.perm = 5, fit.method = "nbinomial", points = points,
                                text = rownames(metalonda_test_data)[5], parall = FALSE, pvalue.threshold = 0.05,     
                                adjust.method = "BH", col = c("black", "green"))


output.metalonda.all = metalondaAll(Count = metalonda_test_data, Time = Time, Group = Group,
                                    ID = ID, n.perm = 100, fit.method = "nbinomial", num.intervals = 100, 
                                    parall = FALSE, pvalue.threshold = 0.05, adjust.method = "BH", time.unit = "hours", 
                                    norm.method = "none", prefix = "Test", ylabel = "Read Counts", col = c("black", "green"))

