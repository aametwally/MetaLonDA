# What is MetaLonDA?

[![Build Status](https://travis-ci.org/aametwally/MetaLonDA.svg?branch=master)](https://travis-ci.org/aametwally/MetaLonDA)


MetaLonDA (METAgenomic LONgitudinal Differential Abundance method) is a method that identifies the significant time intervals of microbial features in longitudianl studies. MetaLonDA has the ability to handle the inconsistencies and common challenges associated with human studies, such as variable sample collection times and uneven number of time points along the subjects’ longitudinal study. The method employs a negative binomial distribution in conjunction with a semi-parametric SS-ANOVA to model the count reads. Then, it performs the significance testing based on unit time intervals using permutation testing procedure.



# Getting Started
This section details steps for installing and running MetaLonDA. If you experience difficulty installing or running the software, please contact (Ahmed Metwally: ametwa2@uic.edu).

## Prerequisites

* R(>= 3.2.0)


## Installing, Testing, and Running

#### Installation:
```
install.packages("MetaLonDA")
```


#### Example:
```
library(MetaLonDA)
data(metalonda_test_data)

n.sample = 5
n.timepoints = 10
n.group = 2
Group = factor(c(rep(0, n.sample*n.timepoints), rep(1,n.sample*n.timepoints)))
Time = rep(rep(1:n.timepoints, times = n.sample), 2)
ID = factor(rep(1:(2*n.sample), each = n.timepoints))
points = seq(1, 10, length.out = 10)

## Test the first feature 
output.nbinomial = metalonda(Count = metalonda_test_data[1,], Time = Time, Group = Group,
  ID = ID, fit.method =  "nbinomial", n.perm = 10, points = points,
  text = rownames(metalonda_test_data)[1], parall = FALSE, pvalue.threshold = 0.05, adjust.method = "BH")

## Test all features
output.nbinomial = metalondaAll(Count = metalonda_test_data, Time = Time, Group = Group,
  ID = ID, n.perm = 9, fit.method =  "nbinomial", num.intervals = 100, 
  parall = FALSE, pvalue.threshold = 0.05, adjust.method = "BH", time.unit = "hours", norm.method = "none",
  prefix = "Test")
  
```


### Bugs and Suggestions
MetaLonDA is under active research development. Please report any bugs/suggestions to Ahmed Metwally (ametwa2@uic.edu).
