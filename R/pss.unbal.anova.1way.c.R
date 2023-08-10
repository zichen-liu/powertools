#' Power calculations for one-way unbalanced analysis of variance omnibus F test
#'
#' @param nvec A vector of group sample sizes c(n1, n2, ...).
#' @param mvec A vector of group mvec c(mu1, mu2, ...).
#' @param cvec A vector of contrast cvecicients c(c1, c2, ...).
#' @param sd The estimated standard deviation within each group; defaults to 1.
#' @param alpha The significance level or type 1 error rate; defaults to 0.05.
#'
#' @return A list of the arguments (including the computed power).
#' @export
#'
#' @examples
#' # Example 5.7
#' pss.anova.1way.c(n = 20, mvec = c(5, 10, 12), cvec = c(1, -1, 0), sd = 10, alpha = 0.025)
#' pss.anova.1way.c(n = 20, mvec = c(5, 10, 12), cvec = c(1, 0, -1), sd = 10, alpha = 0.025)

pss.unbal.anova.1way.c <- function (nvec = NULL, mvec = NULL, sd = NULL, alpha = 0.05) {

  # Check if the arguments are specified correctly
  a <- length(mvec)
  if (a < 2)
    stop("number of a must be at least 2")
  if (any(nvec < 2))
    stop("number of observations in each group must be at least 2")
  if(is.null(nvec) | is.null(mvec) | is.null(cvec))
    stop("sample size vector, mvec vector, and contrast coefficients vector must all be specified")
  if(a != length(nvec))
    stop("number of sample sizes must equal to the number of groups")
  if(a != length(cvec))
    stop("number of contrast coefficients must be equal to the number of groups")

  # Get weighted sum
  N <- sum(nvec)
  mu <- mean(mvec)
  mmA <- mvec - mu

  # Get f effect size
  sdA <- sqrt(sum(mmA^2) / a)
  f <- sdA / sd

  # Get ncp
  props <- nvec / N
  ws <- props %*% mmA
  l <- sapply(X = 1:a, FUN = function(i) nvec[i] * ((mmA[i] - ws) / sd)^2)
  Lambda <- sum(l)

  # Calculate power
  power <- stats::pf(stats::qf(alpha, a - 1, N - a, lower.tail = FALSE),
                     a - 1, N - a, Lambda, lower.tail = FALSE)

  # Generate output text
  METHOD <- "Unbalanced one-way analysis of variance\n     contrast test power calculation"

  # Print output as a power.htest object
  structure(list(a = a, nvec = nvec, mvec = mvec, cvec = cvec,
                 sd = sd, alpha = alpha, power = power,
                 method = METHOD), class = "power.htest")
}
