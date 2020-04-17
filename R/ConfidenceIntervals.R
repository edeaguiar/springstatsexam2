# Confidence Intervals functions
#' Error margin for normal distribution CI with a population standard deviation where n is large
#'
#' @param CI Confidence interval expressed as a decimal
#' @param n Sample size
#' @param sd Population standard deviation
#'
#' @return
#' @export
#'
#' @examples
error.margin <- function(CI, n, sd, large = TRUE)
{
  if (large)
  {
    list(method = "proceed with upper and lower bound as normal",
         errormargin = abs(qnorm((CI + (
           1 - CI
         ) / 2), lower.tail = F) * sd / sqrt(n)))
  }
  else
  {
    list(
      method = "small.sample.size",
      errormargin.studentTdistribution = (qt((1 - CI) / 2, n - 1, lower.tail = F)) *
        sd / sqrt(n)
    )
  }
}
#' Margin of error for difference between population means
#'
#' @param CI Confidence interval expressed as a decimal
#' @param x1 Vector (bigger values)
#' @param x2 Vector (smaller values)
#'
#' @return
#' @export
#'
#' @examples
CI.diff <- function(CI, x1, x2)
{
  (qt((1 - CI) / 2, length(x1) - 1, lower.tail = F)) * sd(x1 - x2) / sqrt(length(x1))
}
#' Margin of Error for difference between paired population means
#'
#' @param CI Confidence interval expressed as a decimal
#' @param ssd1 Sample Standard Deviation 1
#' @param ssd2 Sample Standard Deviation 2
#' @param n1 Sample size 1
#' @param n2 Sample size 2
#'
#' @return
#' @export
#'
#' @examples
CI.paired.diff<-function(CI,ssd1,ssd2,n1,n2,large=TRUE)
{
  if (large)
  {
    large = qnorm(CI + (1 - CI) / 2) * sqrt(ssd1 ^ 2 / n1 + ssd2 ^ 2 / n2)
  }
  else
  {
    small = abs(qt(abs((CI + (
      1 - CI
    ) / 2)), (ssd1 ^ 2 / n1 + ssd2 ^ 2 / n2) ^ 2 / ((ssd1 ^ 2 / n1) ^ 2 / (n1 -
                                                                             1) + (ssd2 ^ 2 / n2) ^ 2 / (n2 - 1)
    ), low = F))
  }
}
#' Margin of error for a single large population proportion
#'
#' @param CI Confidence interval expressed as a decimal
#' @param p1 Number of subjects 1 (numerator of proportion)
#' @param n Sample size
#'
#' @return
#' @export
#'
#' @examples To calculate bounds, subtract/add this functions output to p=(p1+2)/(n+4)
CI.prop<-function(CI,p1,n)
{
  qnorm(CI+(1-CI)/2)*sqrt((p1+2)/(n+4)*(1-(p1+2)/(n+4))/(n+4))
}
#' Margin of error for the difference between two population (paired) proportions, large
#'
#' @param CI Confidence interval, expressed as a decimal
#' @param S1 Number of subjects 1 (numerator of proportion)
#' @param S2 Number of subjects 2 (numerator of proportion)
#' @param n1 Sample size 1
#' @param n2 Sample size 2
#'
#' @return
#' @export
#'
#' @examples To compute p1 and p2 from the S and n for calculating the lower and upper bounds: p1=(s1+1)/(n1+2) and p2=(s2+1)/(n2+2). Subtract p1 (the bigger proportion) from the smaller p2 and then subtract/add the margin of error
CI.paired.prop.diff<-function(CI,S1,S2,n1,n2)
{
  qnorm(CI+(1-CI)/2)*sqrt((S1+1)/(n1+2)*(1-(S1+1)/(n1+2))/(n1+2)+(S2+1)/(n2+2)*(1-(S2+1)/(n2+2))/(n2+2))
}

