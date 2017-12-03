# What is MetaLonDA?
MetaLonDA METAgenomic LONgitudinal Differential Abundance method) is a method that identify the significant time intervals of microbial features in longitudianl studies. MetaLonDA has the ability to handle the inconsistencies and common challenges associated with human studies, such as variable sample collection times and uneven number of time points along the subjects’ longitudinal study. The method employs a negative binomial distribution in conjunction with  a semi-parametric SS-ANOVA to model the count reads. MetaLonDA performs the significance testing based on unit time intervals.


### Publication:
Ahmed A. Metwally, Patricia W. Finn, Yang Dai, and David L. Perkins. "Detection of Differential Abundance Intervals in Longitudinal Metagenomic Data Using Negative Binomial Smoothing Spline ANOVA." ACM BCB (2017)


# Getting Started
This section details steps for installing and running MetaLonDA. If you experience difficulty installing or running the software, please contact (Ahmed Metwally: ametwa2@uic.edu).

## Prerequisites

* R(>= 3.1.2)


## Installing, Testing, and Running

#### Installation:
```
install.packages("MetaLonDA")
library(MetaLonDA)
data(metalonda_test_data)
```


#### Example:
```
n.sample = 5 # sample size;
n.timepoints = 10 # time point;
n.group= 2 # number of group;
Group = factor(c(rep(0,n.sample*n.timepoints), rep(1,n.sample*n.timepoints)))
Time = rep(rep(1:n.timepoints, times = n.sample), 2)
ID = factor(rep(1:(2*n.sample), each = n.timepoints))
points = seq(1, 10, length.out = 10)


output_all_nbinomial = metalondaAll(data = metalonda_test_data, Time = Time, Group = Group, ID = ID, 
                                     log = FALSE, fit.method = "nbinomial", n.perm = 10, points = points, pvalue_threshold=0.05)

output_1_nbinomial = metalonda(Count = metalonda_test_data[1,], Time = Time, Group = Group, ID = ID, log = log,
           fit.method =  "nbinomial", n.perm = 10, points = points,
           text=rownames(metalonda_test_data)[1], parall = FALSE, pvalue_threshold=0.05, adjust.method="BH")

```


### Bugs and Suggestions
MetaLonDA is under active research development. Please report any bugs/suggestions to Ahmed Metwally (ametwa2@uic.edu).
