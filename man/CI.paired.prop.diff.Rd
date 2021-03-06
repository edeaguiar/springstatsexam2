% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/ConfidenceIntervals.R
\name{CI.paired.prop.diff}
\alias{CI.paired.prop.diff}
\title{Margin of error for the difference between two population (paired) proportions, large}
\usage{
CI.paired.prop.diff(CI, S1, S2, n1, n2)
}
\arguments{
\item{CI}{Confidence interval, expressed as a decimal}

\item{S1}{Number of subjects 1 (numerator of proportion)}

\item{S2}{Number of subjects 2 (numerator of proportion)}

\item{n1}{Sample size 1}

\item{n2}{Sample size 2}
}
\value{

}
\description{
Margin of error for the difference between two population (paired) proportions, large
}
\examples{
To compute p1 and p2 from the S and n for calculating the lower and upper bounds: p1=(s1+1)/(n1+2) and p2=(s2+1)/(n2+2). Subtract p1 (the bigger proportion) from the smaller p2 and then subtract/add the margin of error
}
