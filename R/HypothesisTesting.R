#NewHypothesisTesting
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
