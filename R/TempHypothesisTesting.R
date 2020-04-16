#Hypothesis testing functions!

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
Z <- function(xbar, mu, sd, n)
{
  (xbar-mu)/(sd/sqrt(n))
}

#' Calculates test statistic for a difference between paired data population proportions
#'
#' @param p1 How many out of the overall sample size for first sample
#' @param s1 Sample size for first sample
#' @param p2 How many out of the overall sample size for second sample
#' @param s2 Sample size for second sample
#'
#' @return
#' @export
#'
#' @examples
HTest.diff.pop.proportions<-function(p1,s1,p2,s2)
{
  ((p2/s2)-(p1/s1))/(sqrt(((p2+p1)/(s2+s1))*(1-((p2+p1)/(s2+s1)))*(1/s1+1/s2)))
}
#' Hypothesis test z test stat inputs for pnorm for paired data, n is large, population variances unknown and not assumed equal
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
HTest.diff.prop.nopopvariance<-function(xbar1,xbar2,ssd1,ssd2,n1,n2,islarge=TRUE)
{
  if (islarge)
  {
    list(method="pnorm",
         Z=(xbar1-xbar2)/sqrt((ssd1^2/n1)+(ssd2^2/n2)))
  }
  else
  {
    list(method="pt",
         T=(xbar1-xbar2)/sqrt(ssd1/4+ssd2/6),
         Student.T.Dist=(ssd1/4+ssd2/6)^2/((ssd1/4)^2/(4-1)+(ssd2/6)^2/(6-1)))
  }
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
H.prop.l.nv<-function(p1,n,p0)
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
  2*pnorm(abs(Z(xbar, mu, sd, n)),lower.tail = F)
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
  pnorm(Z(xbar, mu, sd, n))
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
  pnorm(Z(xbar, mu, sd, n),lower.tail = F)
}

#' Hypothesis test for a population mean difference VECTOR. Outputs values needed to put into pt
#'
#' @param x Vector
#' @param mu Hypothesizedsample mean
#'
#' @return
#' @export
#'
#' @examples
H.d.s.nv<-function(x,mu)
{
  list((mean(x)-mu)/(sd(x)/sqrt(length(x))),length(x)-1)
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

#' Computes P value when alt =/ null, n is small,variance is known
#' ITS PT SO YOU KNOW ITS A T VALUE NOT Z
#' @param xbar population mean
#' @param mu hypothesized mean
#' @param sd standard dev
#' @param n sample size
#'
#' @return
#' @export
#'
#' @examples
n.small.p<-function(xbar, mu, sd, n)
{
  2*pt(abs(Z(xbar, mu, sd, n)),(n-1),lower.tail = F)
}

#' Computes P value when alt is greater than null, n is small, variance is known
#' ITS PT SO YOU KNOW ITS A T VALUE NOT Z
#' @param xbar pop mean
#' @param mu hyp mean
#' @param sd standard dev
#' @param n sample size
#'
#' @return
#' @export
#'
#' @examples
g.small.pvalue<-function(xbar, mu, sd, n)
{
  pt(Z(xbar, mu, sd, n),(n-1),lower.tail = F)
}

#' Computes P value when alt is less than null, n is small,variance is known
#' ITS PT SO YOU KNOW ITS A T VALUE NOT Z
#' @param xbar pop mean
#' @param mu hyp mean
#' @param sd standard dev
#' @param n sample size
#'
#' @return
#' @export
#'
#' @examples
l.small.pvalue<-function(xbar, mu, sd, n)
{
  pt(Z(xbar, mu, sd, n),(n-1))
}

#' Computes P value for PAIRED DATA when alt greater than null, n is small,variance is unknown
#' ITS PT SO YOU KNOW ITS A T VALUE NOT Z
#' @param x1 Vector with larger values
#' @param x2 Vector with smaller values
#' @param mu hypothesized mean
#'
#' @return
#' @export
#'
#' @examples
g.p.small.p<-function(x1,x2,mu)
{
  pt((mean(x1-x2)-mu)/(sd(x1-x2)/sqrt(length(x1))),length(x1)-1,lower.tail = F)
}
#' Computes P value for PAIRED DATA when alt less than null, n is small,variance is unknown
#' ITS PT SO YOU KNOW ITS A T VALUE NOT Z
#' @param x1 Vector with larger values
#' @param x2 Vector with smaller values
#' @param mu hypothesized mean
#'
#' @return
#' @export
#'
#' @examples
l.p.small.p<-function(x1,x2,mu)
{
  pt((mean(x1-x2)-mu)/(sd(x1-x2)/sqrt(length(x1))),length(x1)-1)
}
#' Computes P value for PAIRED DATA when alt =/ null, n is small, variance is unknown
#' ITS PT SO YOU KNOW ITS A T VALUE NOT Z
#' @param x1 Vector with larger values
#' @param x2 Vector with smaller values
#' @param mu hypothesized mean
#'
#' @return
#' @export
#'
#' @examples
n.p.small.p<-function(x1,x2,mu)
{
  2*pt((mean(x1-x2)-mu)/(sd(x1-x2)/sqrt(length(x1))),length(x1)-1,lower.tail = F)
}
#' Computes effect sice plus p
#'
#' @param es effect size
#' @param p p
#'
#' @return
#' @export
#'
#' @examples
throwawayforcomputingpower<-function(es,p)
{
  es+p
}
#' Threshold for rejecting the null hypothesis
#'
#' @param p p number given
#' @param n sample size
#' @param a significance level
#'
#' @return
#' @export
#'
#' @examples
rejectnullif<-function(p,n,a)
{
  sqrt(p*(1-p)/n)*qnorm(a,lower.tail = F)+p
}
#' Standard deviation of p hat for calculating power level
#'
#' @param es effect size
#' @param p p number given
#' @param n significance level
#'
#' @return
#' @export
#'
#' @examples
sdphatpowerlvl<-function(es,p,n)
{
  sqrt((es+p)*(1-(es+p))/n)
}
#' Determines Power of a hypothesis test given p, effect size, sample size, and significance level
#'
#' @param p p number given
#' @param es effect size
#' @param n sample mean
#' @param a significance level
#'
#' @return
#' @export
#'
#' @examples
detpowerlevel<-function(p,es,n,a)
{
  1-pnorm(rejectnullif(p,n,a),throwawayforcomputingpower(es,p),sdphatpowerlvl(es,p,n))
}
#' Calculates 1. test statistic t and 2. t distribution needed for a p value when calculating for a difference in two population means where the sample size is small and a two vectors are given
#'
#' @param xbar1 Sample mean 1 (smaller)
#' @param xbar2 Sample mean 2 (bigger)
#' @param ssd1 Sample variance 1
#' @param ssd2 Sample variance 2
#'
#' @return
#' @export
#'
#' @examples Then, to obtain p value, plug the two numbers into pt() based on if they want the alt hypothesis to be greater,less than, or unequal to
H.d.p.s.nv<-function(xbar1,xbar2,ssd1,ssd2)
{

}

