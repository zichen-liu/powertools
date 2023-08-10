#' Power calculations for one-way unbalanced analysis of variance omnibus F test
#'
#' @param nvec A vector of group sample sizes c(n1, n2, ...).
#' @param mvec A vector of group mvec c(mu1, mu2, ...).
#' @param sd The estimated standard deviation within each group.
#' @param alpha The significance level or type 1 error rate; defaults to 0.05.
#'
#' @return A list of the arguments (including the computed power).
#' @export
#'
#' @examples
#' # Example 5.2
#' pss.unbal.anova.1way(nvec = c(20, 20, 20), mvec = c(5, 10, 12), sd = 10)

pss.unbal.anova.1way <- function (nvec = NULL, mvec = NULL, sd = NULL, alpha = 0.05) {

  # Check if the arguments are specified correctly
  a <- length(mvec)
  if (a < 2)
    stop("number of a must be at least 2")
  if (any(nvec < 2))
    stop("number of observations in each group must be at least 2")
  if(is.null(nvec) | is.null(mvec))
    stop("sample size vector and mvec vector must both be specified")
  if(is.null(sd))
    stop("sd must be specified")
  if(a != length(nvec))
    stop("number of sample sizes must equal to the number of groups")

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
  temp <- sapply(X = 1:a, FUN = function(i) nvec[i] * ((mmA[i] - ws) / sd)^2)
  Lambda <- sum(temp)

  # Calculate power
  power <- stats::pf(stats::qf(alpha, a - 1, N - a, lower.tail = FALSE),
                     a - 1, N - a, Lambda, lower.tail = FALSE)

  # Generate output text
  METHOD <- "Unbalanced one-way analysis of variance\n     omnibus F test power calculation"

  # Print output as a power.htest object
  structure(list(a = a, nvec = nvec, mvec = mvec,
                 sd = sd, f = f, alpha = alpha, power = power,
                 method = METHOD), class = "power.htest")
}
