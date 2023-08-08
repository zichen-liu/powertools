#' Power calculations for one sample t tests
#'
#' @param n The sample size.
#' @param delta muA (the true mean) - mu0 (the mean under the null).
#' @param sd The estimated standard deviation; defaults to 1.
#' @param alpha The significance level or type 1 error rate; defaults to 0.05.
#' @param power The specified level of power.
#' @param sided Either "one" or "two" (default) to specify a one- or two- sided hypothesis test.
#' @param strict Use strict interpretation in two-sided case; defaults to TRUE.
#'
#' @return A list of the arguments (including the computed one).
#' @export
#'
#' @examples
#' # Example 3.2
#' pss.t.test.1samp(n = 36, delta = 4.9 - 5.7, sd = 2, sided = "one")
#' # Example 3.3
#' pss.t.test.1samp(n = 36, delta = 6.3 - 5.7, sd = 2, sided = "one")
#' # Example 3.5
#' pss.t.test.1samp(n = 36, delta = 4.9 - 5.7, sd = 2, sided = "two")
#' # Example 3.6
#' pss.t.test.1samp(delta = 0.6, sd = 1, power = 0.8, sided = "one")

pss.t.test.1samp <- function (n = NULL, delta = NULL, sd = 1,
                              alpha = 0.05, power = NULL,
                              sided = c("two", "one"), strict = TRUE) {

  # Check if the arguments are specified correctly
  if (sum(sapply(list(n, delta, sd, power, alpha), is.null)) != 1)
    stop("exactly one of n, delta, sd, alpha, and power must be NULL")

  # Assign number of sides
  sided <- match.arg(sided)
  side <- switch(sided, one = 1, two = 2)

  # Use absolute value of the effect size
  if (!is.null(delta))
    delta <- abs(delta)

  # Calculate df and ncp
  p.body <- quote({
    stats::pt(stats::qt(alpha / side, n - 1, lower.tail = FALSE), n - 1,
              sqrt(n) * delta / sd, lower.tail = FALSE)
  })

  # The strict two-sided case uses the F distribution
  if (strict && side == 2)
    p.body <- quote({
      ncp <- sqrt(n) * delta / sd
      stats::pf(stats::qf(alpha, 1, n - 1, lower.tail = FALSE),
                1, n - 1, ncp^2, lower.tail = FALSE)
    })

  # Use uniroot function to calculate missing argument
  if (is.null(power))
    power <- eval(p.body)
  else if (is.null(n))
    n <- uniroot(function(n) eval(p.body) - power, c(2, 1e+07))$root
  else if (is.null(sd))
    sd <- uniroot(function(sd) eval(p.body) - power, delta * c(1e-07, 1e+07))$root
  else if (is.null(delta))
    delta <- uniroot(function(delta) eval(p.body) - power,  c(1e-07, 1e+07))$root
  else if (is.null(alpha))
    alpha <- uniroot(function(alpha) eval(p.body) - power, c(1e-10, 1 - 1e-10))$root
  else stop("internal error")

  # Generate output text
  METHOD <- "One-sample t test power calculation"

  # Print output as a power.htest object
  structure(list(n = n, delta = delta, sd = sd, alpha = alpha,
                 power = power, sided = sided,
                 method = METHOD), class = "power.htest")
}
