#' Power calculations for one-way balanced analysis of variance omnibus F test
#'
#' @param n The sample size per group.
#' @param mvec A vector of group means c(mu1, mu2, ...).
#' @param sd The estimated standard deviation within each group; defaults to 1.
#' @param alpha The significance level or type 1 error rate; defaults to 0.05.
#' @param power The specified level of power.
#'
#' @return A list of the arguments (including the computed one).
#' @export
#'
#' @examples
#' # Example 5.2
#' pss.anova.1way(n = 20, mvec = c(5, 10, 12), sd = 10)

pss.anova.1way <- function (n = NULL, mvec = NULL, sd = 1,
                             alpha = 0.05, power = NULL) {

  # Check if the arguments are specified correctly
  a <- length(mvec)
  if (sum(vapply(list(n, alpha, power), is.null, NA)) != 1)
    stop("exactly one of 'n', 'alpha', and 'power' must be NULL")
  if (a < 2)
    stop("number of groups must be at least 2")
  if (!is.null(n) && n < 2)
    stop("number of observations in each group must be at least 2")
  if(is.null(sd))
    stop("sd must be specified")

  # Get between group variance; sd is within group standard deviation
  within <- sd^2
  between <- var(mvec)
  f <- sqrt(between / within)

  # Calculate df and ncp
  p.body <- quote({
    lambda <- (a - 1) * n * f^2
    stats::pf(stats::qf(alpha, a - 1, (n - 1) * a, lower.tail = FALSE),
              a - 1, (n - 1) * a, lambda, lower.tail = FALSE)
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
  structure(list(a = a, n = n, mvec = mvec,
                 sd = sd, f = f, alpha = alpha, power = power,
                 note = NOTE, method = METHOD), class = "power.htest")
}
