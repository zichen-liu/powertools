#' Power calculations for one-way unbalanced analysis of variance omnibus F test
#'
#' @param nvec A vector of group sample sizes.
#' @param means A vector of group means c(mu1, mu2, ...).
#' @param sd The estimated standard deviation within each group; defaults to 1.
#' @param alpha The significance level or type 1 error rate; defaults to 0.05.
#'
#' @return A list of the arguments (including the computed power).
#' @export
#'
#' @examples
#' # Example 5.2
#' pss.unbal.anova.1way(nvec = c(20, 20, 20), means = c(5, 10, 12), sd = 10)

pss.unbal.anova.1way <- function (nvec = NULL, means = NULL, sd = NULL, alpha = 0.05) {

  # Check if the arguments are specified correctly
  groups <- length(means)
  if (groups < 2)
    stop("number of groups must be at least 2")
  if (any(nvec < 2))
    stop("number of observations in each group must be at least 2")
  if(is.null(nvec) | is.null(means))
    stop("sample size vector and means vector must both be specified")
  if(groups != length(nvec))
    stop("number of sample sizes must equal to the number of groups")

  # Get weighted sum
  N <- sum(nvec)
  mu <- mean(means)
  alphas <- means - mu

  # Get f effect size
  sdA <- sqrt(sum(alphas^2) / groups)
  f <- sdA / sd

  # Get ncp
  props <- nvec / N
  ws <- props %*% alphas
  l <- sapply(X = 1:groups, FUN = function(i) nvec[i] * ((alphas[i] - ws) / sd)^2)
  Lambda <- sum(l)

  # Calculate power
  power <- stats::pf(stats::qf(alpha, groups - 1, N - groups, lower.tail = FALSE),
                     groups - 1, N - groups, Lambda, lower.tail = FALSE)

  # Generate output text
  METHOD <- "Unbalanced one-way analysis of variance\n     omnibus F test power calculation"

  # Print output as a power.htest object
  structure(list(groups = groups, nvec = nvec, means = means,
                 sd = sd, f = f, alpha = alpha, power = power,
                 method = METHOD), class = "power.htest")
}
