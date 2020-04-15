#' Calculates basic test statistic, variance is known, n is large
#'
#' @param xbar Sample mean
#' @param mu Hypothesized sample mean
#' @param sd Standard deviation
#' @param n Sample size
#'
#' @return The test statistic
#' @export
#'
#' @examples
heart <- function(xbar, mu, sd, n)
{
  (xbar-mu)/(sd/sqrt(n))
}

#' Calculates test statistic for paired data, variance is unknown, n is large
#'
#' @param p1 How many out of the sample size for first data
#' @param s1 Sample size for first data
#' @param p2 How many out of the sample size for first data
#' @param s2 Sample size for second data
#'
#' @return
#' @export
#'
#' @examples
ptstatnissmall<-function(p1,s1,p2,s2)
{
  ((p2/s2)-(p1/s1))/(sqrt(((p2+p1)/(s2+s1))*(1-((p2+p1)/(s2+s1)))*(1/s1+1/s2)))
}
#' Calculates test statistic for paired data, variances unknown and not assumed equal
#'
#' @param mu1 Average for first data
#' @param mu2 Average for second data
#' @param sd1 Standard Deviation of the first data
#' @param sd2 Standard Deviation of the second data
#' @param n1 First sample size
#' @param n2 Second sample size
#'
#' @return
#' @export
#'
#' @examples
ptstatvarunequal<-function(mu1,mu2,sd1,sd2,n1,n2)
{
  (mu1-mu2)/sqrt((sd1^2/n1)+(sd2^2/n2))
}
#' Calculates test statistic, variance is unknown, N is large
#'
#' @param p1 The numerator of the actual proportion
#' @param n The sample size (denominator of the actual proportion)
#' @param p0 The wanted proportion as a decimal
#'
#' @return
#' @export
#'
#' @examples
tstatnislarge<-function(p1,n,p0)
{
  ((p1/n)-p0)/sqrt(p0*(1-p0)/n)
}

#' Standard deviation p hat for significance level
#'
#' @param p proportion to be tested, a decimal
#' @param n Sample size
#'
#' @return
#' @export
#'
#' @examples
sdphatsiglvl<-function(p,n)
{
  sqrt(p*(1-p)/n)
}

#' Calculates the value of p hat for calculating significance level
#'
#' @param H0 null hypothesis threshold
#' @param n sample size
#'
#' @return
#' @export
#'
#' @examples
phatsiglevel<-function(H0,n)
{
  H0/n
}

#' Determines the significance level
#'
#' @param p proportion as a decimal
#' @param H0 null hypothesis threshold
#' @param n Sample size
#'
#' @return
#' @export
#'
#' @examples
detsiglevel<-function(p,H0,n)
{
  pnorm(phatsiglevel(H0,n),p,sdphatsiglvl(p,n))
}

#' Computes P-value when the alternate hypothesis =/ the null
#'
#' @param xbar Sample mean
#' @param mu Hypothesized sample mean
#' @param sd Standard Deviation
#' @param n Sample Size
#'
#' @return
#' @export
#'
#' @examples
altnotnull<-function(xbar,mu,sd,n)
{
  2*pnorm(abs(heart(xbar, mu, sd, n)),lower.tail = F)
}

#' Computes P-value when the alternate hypothesis is less than the null
#'
#' @param xbar Sample mean
#' @param mu Hypothesized sample mean
#' @param sd Standard Deviation
#' @param n Sample Size
#'
#' @return
#' @export
#'
#' @examples
altlessthannull<-function(xbar,mu,sd,n)
{
  pnorm(heart(xbar, mu, sd, n))
}

#' Calculates P-Value when the alternate hypothesis is greater than the null
#'
#' @param xbar Sample mean
#' @param mu Hypothesized sample mean
#' @param sd Standard Deviation
#' @param n Sample Size
#'
#' @return
#' @export
#'
#' @examples
altgreaterthannull<-function(xbar, mu, sd, n)
{
  pnorm(heart(xbar, mu, sd, n),lower.tail = F)
}

#' T for computing P value of a vector with a hypothesized mean
#'
#' @param x Vector
#' @param mu Hypothesizedsample mean
#'
#' @return
#' @export
#'
#' @examples
unknownpopvarsmallsamplesizetvalue<-function(x,mu)
{
  (mean(x)-mu)/(sd(x)/sqrt(length(x)))
}

#' Unnecessary degrees of freedom function
#'
#' @param x Vector
#'
#' @return
#' @export
#'
#' @examples
degfree<-function(x)
{
  length(x)-1
}

#' Calculates P value when the alternate hypothesis is greater than null, sample size is small, population variance is unknown
#'
#' @param x Vector
#' @param mu Hypothesized sample mean
#'
#' @return
#' @export
#'
#' @examples
gunknown<-function(x,mu)
{
  pt(unknownpopvarsmallsamplesizetvalue(x,mu),degfree(x),lower.tail = F)
}

#' Calculates P value when the alternate hypothesis is less than null, sample size is small, population variance is unknown
#'
#' @param x Vector
#' @param mu Hypothesized sample mean
#'
#' @return
#' @export
#'
#' @examples
lunknown<-function(x,mu)
{
  pt(unknownpopvarsmallsamplesizetvalue(x,mu),degfree(x))
}

#' Calculates P value when the alternate hypothesis =/ null, sample size is small, population variance is unknown
#'
#' @param x Vector
#' @param mu Hypothesized sample mean
#'
#' @return
#' @export
#'
#' @examples
nunknown<-function(x,mu)
{
  2*pt(unknownpopvarsmallsamplesizetvalue(x,mu),degfree(x),lower.tail = F)
}

#' Calculates t test for paired vectors
#'
#' @param x1 First Vector
#' @param x2 Second Vector
#'
#' @return
#' @export
#'
#' @examples
paireddatasmallsample<-function(x1,x2)
{
  t.test(x1,x2,paired = T)
}






