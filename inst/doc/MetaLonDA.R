## ------------------------------------------------------------------------
library(MetaLonDA)

## Load read counts of 8 features from 100 samples. Samples are from 2 groups, 5 subjects per group, and 10 time points per subject.
data(metalonda_test_data)

## ------------------------------------------------------------------------
## Create Group, Time, and ID annotation vectors
n.group = 2
n.sample = 5 
n.timepoints = 10
Group = factor(c(rep("A", n.sample*n.timepoints), rep("B",n.sample*n.timepoints)))
Time = rep(rep(1:n.timepoints, times = n.sample), 2)
ID = factor(rep(1:(2*n.sample), each = n.timepoints))

## Define the prediction timeponits 
points = seq(1, 10, length.out = 100)
output.metalonda.f5 = metalonda(Count = metalonda_test_data[5,], Time = Time, Group = Group,
                                ID = ID, n.perm = 20, fit.method = "nbinomial", points = points,
                                text = rownames(metalonda_test_data)[5], parall = FALSE, pvalue.threshold = 0.05,     
                                adjust.method = "BH", time.unit = "hours", ylabel = "Read Counts", col = c("chartreuse",
                                                                                                           "blue4"))

## ------------------------------------------------------------------------
## Identify significant time intervals for all features: 
output.metalonda.all = metalondaAll(Count = metalonda_test_data, Time = Time, Group = Group,
                                    ID = ID, n.perm = 20, fit.method = "nbinomial", num.intervals = 100, 
                                    parall = FALSE, pvalue.threshold = 0.05, adjust.method = "BH", time.unit = "hours", 
                                    norm.method = "none", prefix = "Test", ylabel = "Read Counts", col = c("chartreuse",
                                                                                                           "blue4"))

