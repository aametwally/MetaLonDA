setwd("~/Dropbox/metalonda_work/updateMetalonda/MetaLonDA")
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


output.metalonda.f5 = metalonda(Count = metalonda_test_data[5,], Time = Time, Group = Group,
                                ID = ID, n.perm = 5, fit.method = "nbinomial", points = points,
                                text = rownames(metalonda_test_data)[5], parall = FALSE, pvalue.threshold = 0.05,     
                                adjust.method = "BH")


output.metalonda.all = metalondaAll(Count = metalonda_test_data, Time = Time, Group = Group,
                                    ID = ID, n.perm = 5, fit.method = "nbinomial", num.intervals = 100, 
                                    parall = FALSE, pvalue.threshold = 0.05, adjust.method = "BH", time.unit = "hours", 
                                    norm.method = "none", prefix = "Test")








# data <- read.csv("DataTable_Ahmed.csv", as.is = TRUE, check.names = FALSE)
# 
# ### Remove samples with NAs
# data = data[ , apply(data, 2, function(x) !any(is.na(x)))]
# 
#  
# # Extract Time, Group, and ID vector
# time.vec = as.integer(colnames(data)[-1])
# group.vec = factor(as.character(data[1,])[-1]) 
# id.vec = factor(as.character(data[2,])[-1])
# 
# ### remove last token from subject IDs
# id.vec = gsub("R.F.","RF_", id.vec)
# id.vec = gsub("\\..*","", id.vec)
# length(id.vec)
# length(unique(id.vec)) ## Seems that after the filterationw e end up with only 15 subjects!!!!
# 
# 
# 
# ## Extract countTable
# countTable = data[-c(1,2),]
# colnames(countTable) = as.character(data[3,])
# countTable = countTable[-1,]
# rownames(countTable) = countTable[,1]
# countTable = countTable[,-1]
# 
# 
# ### convert CountTable to integer
# countTable2 <- data.frame(sapply(countTable, function(x) as.numeric(as.character(x))))
# rownames(countTable2) = rownames(countTable)
# 
# # Check dimentions
# length(time.vec)
# length(group.vec)
# length(id.vec)
# dim(countTable)
# 
# 
# output.metalonda.all = metalondaAll(Count = countTable2, Time = time.vec, Group = group.vec, ID = id.vec, n.perm = 100, 
#                                     fit.method = "nbinomial", num.intervals = 100, parall = FALSE, 
#                                     pvalue.threshold = 0.05, adjust.method = "none", 
#                                     time.unit = "Days???", norm.method = "none", prefix = "RF_BH_Con_L6_timeline2")
