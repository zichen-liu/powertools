#' Power calculations for one-way balanced analysis of variance omnibus F test
#'
#' @param n The sample size per group.
#' @param means A vector of group means c(mu1, mu2, ...).
#' @param sd The estimated standard deviation within each group.
#' @param alpha The significance level or type 1 error rate; defaults to 0.05.
#' @param power The specified level of power.
#'
#' @return A list of the arguments (including the computed one).
#' @export
#'
#' @examples
#' # Example 5.2
#' pss.anova.1way(n = 20, means = c(5, 10, 12), sd = 10)

pss.anova.1way <- function (n = NULL, means = NULL, sd = NULL,
                             alpha = 0.05, power = NULL) {

  # Check if the arguments are specified correctly
  groups <- length(means)
  if (sum(vapply(list(n, sd, alpha, power), is.null, NA)) != 1)
    stop("exactly one of 'n', 'sd', 'alpha', and 'power' must be NULL")
  if (groups < 2)
    stop("number of groups must be at least 2")
  if (!is.null(n) && n < 2)
    stop("number of observations in each group must be at least 2")

  # Get between group variance; sd is within group standard deviation
  within <- sd^2
  between <- var(means)

  # Calculate df and ncp
  p.body <- quote({
    lambda <- (groups - 1) * n * (between / within)
    stats::pf(stats::qf(alpha, groups - 1, (n - 1) * groups, lower.tail = FALSE),
              groups - 1, (n - 1) * groups, lambda, lower.tail = FALSE)
  })

  # Use uniroot function to calculate missing argument
  if (is.null(power))
    power <- eval(p.body)
  else if (is.null(n))
    n <- uniroot(function(n) eval(p.body) - power, c(2, 1e+05))$root
  else if (is.null(sd)) {
    within <- uniroot(function(within) eval(p.body) - power, between * c(1e-07, 1e+07))$root
    sd <- sqrt(within)
  }
  else if (is.null(alpha))
    alpha <- uniroot(function(alpha) eval(p.body) - power, c(1e-10, 1 - 1e-10))$root
  else stop("internal error", domain = NA)

  # Generate output text
  NOTE <- "n is the number in each group"
  METHOD <- "Balanced one-way analysis of variance\n     omnibus F test power calculation"

  # Print output as a power.htest object
  structure(list(groups = groups, n = n, means = means,
                 sd = sd, alpha = alpha, power = power,
                 note = NOTE, method = METHOD), class = "power.htest")
}
