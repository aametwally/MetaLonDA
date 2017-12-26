# What is MetaLonDA?

[![Build Status](https://travis-ci.org/aametwally/MetaLonDA.svg?branch=master)](https://travis-ci.org/aametwally/MetaLonDA)


MetaLonDA (METAgenomic LONgitudinal Differential Abundance method) is a method that identifies the significant time intervals of microbial features in longitudianl studies. MetaLonDA has the ability to handle the inconsistencies and common challenges associated with human studies, such as variable sample collection times and uneven number of time points along the subjects’ longitudinal study. The method employs a negative binomial distribution in conjunction with a semi-parametric SS-ANOVA to model the count reads. Then, it performs the significance testing based on unit time intervals using permutation testing procedure.



## Publication:
* Ahmed A. Metwally, Patricia W. Finn, Yang Dai, and David L. Perkins. "Detection of Differential Abundance Intervals in Longitudinal Metagenomic Data Using Negative Binomial Smoothing Spline ANOVA." ACM BCB (2017) 




# Getting Started
This section details steps for installing and running MetaLonDA. If you experience difficulty installing or running the software, please contact (Ahmed Metwally: ametwa2@uic.edu).

## Prerequisites

* R(>= 3.2.0)


## Installing and Running

#### Installation:
```
install.packages("MetaLonDA")
```


#### Example:
```
library(MetaLonDA)


## Load read counts of 9 features from 100 samples. Samples are from 2 groups, 5 subjects per group, and 10 time points per subject.
data(metalonda_test_data)

## Create Group, Time, and ID annotation vectors
n.sample = 5
n.timepoints = 10
n.group = 2
Group = factor(c(rep(0, n.sample*n.timepoints), rep(1,n.sample*n.timepoints)))
Time = rep(rep(1:n.timepoints, times = n.sample), 2)
ID = factor(rep(1:(2*n.sample), each = n.timepoints))

## Define the prediction timeponits 
points = seq(1, 10, length.out = 100)

## Identify significant time intervals of the first feature: 
output.metalonda.f1 = metalonda(Count = metalonda_test_data[1,], Time = Time, Group = Group,
  ID = ID, fit.method = "nbinomial", n.perm = 20, points = points,
  text = rownames(metalonda_test_data)[1], parall = FALSE, pvalue.threshold = 0.05, adjust.method = "BH")

## Identify significant time intervals for all features: 
output.metalonda.all = metalondaAll(Count = metalonda_test_data, Time = Time, Group = Group,
  ID = ID, n.perm = 20, fit.method = "nbinomial", num.intervals = 100, 
  parall = FALSE, pvalue.threshold = 0.05, adjust.method = "BH", time.unit = "hours", norm.method = "none",
  prefix = "Test")
  
```


### Bugs and Suggestions
MetaLonDA is under active research development. Please report any bugs/suggestions to Ahmed Metwally (ametwa2@uic.edu).
