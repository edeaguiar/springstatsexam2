#' Component of another function for non normal distribution determining power
#'
#' @param nnd ADD THE ^2 not normal distribution given ADD THE ^2
#' @param n sample size
#'
#' @return
#' @export
#'
#' @examples
zignore<-function(nnd,n)
{
  (sqrt(nnd)/sqrt(n))
}
#' NOT NORMAL DISTRIBUTION Threshold for rejecting the null hypothesis for determining power level
#'
#' @param nnd not normal distribution
#' @param n sample size
#' @param a significance level
#' @param p p value given
#'
#' @return
#' @export
#'
#' @examples
NOTREGDISTrejectnullif<-function(nnd,n,a,p)
{
  (sqrt(nnd)/sqrt(n))*qnorm(a,lower.tail = F)+p
}
#' Determines power for hypothesis test with a NOT NORMAL DISTRIBUTION o
#' given the p value, effect size, not normal distribution SQUARE IT, sample size, effect size, and significance level
#' @param p
#' @param nnd
#' @param n
#' @param es
#' @param a
#'
#' @return
#' @export
#'
#' @examples
detpowerlevelNOTNORMALDIST<-function(p,nnd,n,es,a)
{
  1-pnorm(NOTREGDISTrejectnullif(nnd,n,a,p),throwawayforcomputingpower(es,p),zignore(nnd,n))
}
#' Calculates error margin with a sample standard dev and a small sample size
#'
#' @param CI Confidence Interval, as a decimal
#' @param n Sample size
#' @param samplesd Sample standard deviation
#'
#' @return
#' @export
#'
#' @examples
s.nv.errormargin<-function(CI,n,samplesd)
{
  (qt((1-CI)/2,n-1,lower.tail = F))*samplesd/sqrt(n)
}
#' Calculates upper bound for a confidence interval with no variance and a small sample
#'
#' @param xbar population mean aka middle of normal distribution
#' @param CI Confidence interval, as a decimal
#' @param n Sample size
#' @param samplesd Sample standard deviation
#'
#' @return
#' @export
#'
#' @examples
s.nv.upperboundCI<-function(xbar,CI,n,samplesd)
{
  xbar+s.nv.errormargin(CI,n,samplesd)
}
#' Calculates lower bound for a confidence interval with no variance and a small sample
#'
#' @param xbar Population mean aka middle of normal distribution
#' @param CI Confidence interval, as a decimal
#' @param n Sample size
#' @param samplesd Sample standard deviation
#'
#' @return
#' @export
#'
#' @examples
s.nv.lowerboundCI<-function(xbar,CI,n,samplesd)
{
  xbar-s.nv.errormargin(CI,n,samplesd)
}
#' Calculates upper bound for CI with a non normal distribution
#' sample is large
#' @param xbar
#' @param CI
#' @param ssd
#' @param n
#'
#' @return
#' @export
#'
#' @examples
nnd.l.v.upperbound<-function(xbar,CI,ssd,n)
{
  xbar+nnd.l.v.errormargin(CI,ssd,n)
}
#' Calculated error margin of a confidence interval with a non normal distribution
#'
#' @param CI Confidence interval, as a decimal
#' @param ssd Sample standard deviation
#' @param n Sample size
#'
#' @return
#' @export
#'
#' @examples
nnd.l.v.errormargin<-function(CI,ssd,n)
{
  abs(qnorm((CI+(1-CI)/2),lower.tail = F)*ssd/sqrt(n))
}
#' Calculates lower bound for CI with a non normal distribution
#'
#' @param xbar Population mean
#' @param CI Confidence interval, expressed as a decimal
#' @param ssd Sample standard deviation
#' @param n Sample size
#'
#' @return
#' @export
#'
#' @examples
nnd.l.v.lowerbound<-function(xbar,CI,ssd,n)
{
  xbar-nnd.l.v.errormargin(CI,ssd,n)
}
#' Normal distribution of a confidence interval with a population standard deviation
#'
#' @param CI Confidence interval, expressed as a decimal
#' @param n Sample size
#' @param population.sd Population standard deviation
#'
#' @return
#' @export
#'
#' @examples
ndist.v.errormargin<-function(CI,n,population.sd)
{
  abs(qnorm((CI+(1-CI)/2),lower.tail = F)*population.sd/sqrt(n))
}
#' Calculates how big a sample size needs to be, given the # a parameter is not to exceed, a confidence interval, and a sample standard dev.
#'
#' @param CI Confidence interval, expressed as a decimal
#' @param ssd Sample standard deviation
#' @param no.larger.than The other (and probably new) number given in the word problem
#'
#' @return
#' @export
#'
#' @examples
CI.sample.size<-function(CI,ssd,no.larger.than)
{
  abs(qnorm((CI+(1-CI)/2),lower.tail=F)*ssd/no.larger.than)^2
}
pooledsd<-function(X,Y)
{
  sqrt((length(X)-1)*(sd(X)^2)+(length(Y)-1)*sd(Y)^2/(length(X)+length(Y)-2))
}
t.p.l.nv<-function(X,Y)
{
  (mean(X-Y))/(pooledsd(X,Y)*sqrt((1/length(X))+(1/length(Y))))
}
df.p.s.nv<-function(ssd1,ssd2,n1,n2)
  + {
  (ssd1^2/n1+ssd2^2/n2)^2/((ssd1^2/n1)^2/(n1-1)+(ssd2^2/n2)^2/(n2-1))
}
p.s.nv.CI.t.dist<-function(CI,ssd1,ssd2,n1,n2)
{
  abs(qt(abs((CI+(1-CI)/2)),df.p.s.nv(ssd1,ssd2,n1,n2),low=F))
}
#' Margin of error for the difference between two population means when sample size is small
#'
#' @param CI Confidence interval, expressed as a decimal
#' @param ssd1 Sample Standard Deviation 1
#' @param ssd2 Sample Standard Deviation 2
#' @param n1 Sample size 1
#' @param n2 Sample size 2
#'
#' @return
#' @export
#'
#' @examples
p.s.nv.T.errormargin<-function(CI,ssd1,ssd2,n1,n2)
{
  p.s.nv.CI.t.dist(CI,ssd1,ssd2,n1,n2)*(sqrt((ssd1^2/n1)+(ssd2^2/n2)))
}
#' Margin of Error for the difference between two population means
#'
#' @param CI Confidence interval, expressed as a decimal
#' @param ssd1 Sample Standard Deviation 1
#' @param ssd2 Sample Standard Deviation 2
#' @param n1 Sample size 1
#' @param n2 Sample size 2
#'
#' @return
#' @export
#'
#' @examples
CI.p.l.nv.errormargin<-function(CI,ssd1,ssd2,n1,n2)
{
  qnorm(CI+(1-CI)/2)*sqrt(ssd1^2/n1+ssd2^2/n2)
}
#' Margin of error for a single population proportion
#'
#' @param CI
#' @param p1
#' @param n
#'
#' @return
#' @export
#'
#' @examples To calculate bounds, subtract/add this functions output to p=(p1+2)/(n+4)
CI.l.nv.prop.errormargin<-function(CI,p1,n)
{
  qnorm(CI+(1-CI)/2)*sqrt((p1+2)/(n+4)*(1-(p1+2)/(n+4))/(n+4))
}
#' Margin of error for the difference between two population proportions
#'
#' @param CI
#' @param S1
#' @param S2
#' @param n1
#' @param n2
#'
#' @return
#' @export
#'
#' @examples To compute p1 and p2 from the S and n for calculating the lower and upper bounds: p1=(s1+1)/(n1+2) and p2=(s2+1)/(n2+2). Subtract p1 (the bigger proportion) from the smaller p2 and then subtract/add the margin of error
CI.p.l.nv.prop.errormargin<-function(CI,S1,S2,n1,n2)
{
  qnorm(CI+(1-CI)/2)*sqrt((S1+1)/(n1+2)*(1-(S1+1)/(n1+2))/(n1+2)+(S2+1)/(n2+2)*(1-(S2+1)/(n2+2))/(n2+2))
}

