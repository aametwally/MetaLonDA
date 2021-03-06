% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/Visualization.R
\name{visualizeTimeIntervals}
\alias{visualizeTimeIntervals}
\title{Visualize all significant time intervals for all tested features}
\usage{
visualizeTimeIntervals(
  interval.details,
  prefix = "Test",
  unit = "days",
  col = c("blue", "firebrick")
)
}
\arguments{
\item{interval.details}{Dataframe has infomation about significant interval (feature name, start, end, dominant, p-value)}

\item{prefix}{prefix for the output figure}

\item{unit}{time unit used in the Time vector (hours, days, weeks, months, etc.)}

\item{col}{two color to be used for the two groups (eg., c("red", "blue")).}
}
\description{
Visualize all significant time intervals for all tested features
}
\references{
Ahmed Metwally (ametwall@stanford.edu)
}
