#NewHypothesisTesting
#' Calculates basic test statistic Z or T
#'
#' @param xbar Sample mean
#' @param mu Hypothesized sample mean
#' @param sd Standard deviation
#' @param n Sample size
#' @param large Treat the sample as large?
#' @param comparison notequal, greaterthan or lessthan
#'
#' @return The test statistic
#' @export
#'
#' @examples
test.stat <- function(xbar, mu, sd, n, large = TRUE, comparison="notequal")
{
  RFunc = ""
  if (comparison == "notequal" |
      comparison == "ne" |
      comparison == "<>" | comparison == "!=") {
    RFunc = "2x + low F"
  }
  else if (comparison == "lessthan" |
           comparison == "lt" | comparison == "<") {
    RFunc = NA
  }
  else if (comparison == "greaterthan" |
           comparison == "gt" | comparison == ">") {
    RFunc = "low F"
  }
  else {
    stop("Unknown comparison")
  }

  if (large)
  {
    list(method = "pnorm",
         R.Func = RFunc,
         test.stat.Z = (xbar - mu) / (sd / sqrt(n)))
  }
  else
  {
    list(
      method = "pt",
      R.Func = RFunc,
      test.stat.T = (xbar - mu) / (sd / sqrt(n)),
      degrees.free = (n - 1))
  }
}

another_function<-function()
{

}

#' Hypothesis test for a VECTOR that is population mean difference. Outputs values needed to put into pt
#'
#' @param x Vector
#' @param mu Hypothesized sample mean
#'
#' @return
#' @export
#'
#' @examples
H.diff.small.nopopvariance<-function(x,mu)
{
  list((mean(x)-mu)/(sd(x)/sqrt(length(x))),length(x)-1)
}

#' Calculates Z test stat for a difference between two (paired) population proportions, population variances unknown and not assumed equal. Gives inputs for pnorm
#'
#' @param xbar1 Average for bigger data
#' @param xbar2 Average for smaller data
#' @param ssd1 Sample Standard Deviation of the bigger data
#' @param ssd2 Sample Standard Deviation of the smaller data
#' @param n1 Bigger sample size
#' @param n2 Smaller sample size
#'
#' @return
#' @export
#'
#' @examples
H.diff.prop.nopopvariance<-function(xbar1,xbar2,ssd1,ssd2,n1,n2,islarge=TRUE)
{
  if (islarge)
  {
    list(method="pnorm",
         Z=(xbar1-xbar2)/sqrt((ssd1^2/n1)+(ssd2^2/n2)))
  }
  else
  {
    list(method="pt",
         T=(xbar1-xbar2)/sqrt(ssd1^2/n1+ssd2^2/n2),
         Student.T.Dist=(ssd1^2/n1+ssd2^2/n2)^2/((ssd1^2/n1)^2/(n1-1)+(ssd2^2/n2)^2/(n2-1)))
  }
}

#' Calculates Z test statistic for a population proportion, variance is unknown, N is large
#'
#' @param p1 The numerator of the actual proportion
#' @param n The sample size (denominator of the actual proportion)
#' @param p0 The wanted proportion as a decimal
#'
#' @return
#' @export
#'
#' @examples
H.prop.l.nv<-function(p1,n,p0)
{
  ((p1/n)-p0)/sqrt(p0*(1-p0)/n)
}

#' Calculates Z test statistic for a difference between two (paired) popultation proportions
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
H.paired.prop.l.nv<-function(p1,s1,p2,s2)
{
  ((p2/s2)-(p1/s1))/(sqrt(((p2+p1)/(s2+s1))*(1-((p2+p1)/(s2+s1)))*(1/s1+1/s2)))
}

#' Calculates t test for paired small vectors
#'
#' @param x1 First Vector
#' @param x2 Second Vector
#'
#' @return
#' @export
#'
#' @examples
H.Vectorsgiven.paired.small.nopopvariance<-function(x1,x2)
{
  t.test(x1,x2,paired = T)
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
z.ignore.calculates.p.hat.sd.for.detsiglevel<-function(p,n)
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
z.ignore.calculates.p.hat.for.detsiglevel<-function(H0,n)
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
  pnorm(z.ignore.calculates.p.hat.for.detsiglevel(H0,n),p,z.ignore.calculates.p.hat.sd.for.detsiglevel(p,n))
}

