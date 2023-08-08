#' Power calculations for paired z tests
#'
#' @param n The sample size; the number of pairs.
#' @param delta DeltaA (the true mean difference) - Delta0 (the mean difference under the null).
#' @param sd.pre The estimated pre standard deviation; defaults to 1.
#' @param sd.post The estimated post standard deviation; defaults to 1.
#' @param rho The estimated correlation between pre and post measurements on the same individual.
#' @param alpha The significance level or type 1 error rate; defaults to 0.05.
#' @param power The specified level of power.
#' @param sided Either "one" or "two" (default) to specify a one- or two- sided hypothesis test.
#' @param strict Use strict interpretation in two-sided case; defaults to TRUE.
#'
#' @return A list of the arguments (including the computed one).
#' @export
#'
#' @examples
#' # Example 3.15 approximation
#' pss.z.test.paired(n = NULL, delta = 4, sd.pre = 10, sd.post = 10, rho = 0.4, power = 0.8, sided = "two")

pss.z.test.paired <- function (n = NULL, delta = NULL,
                               sd.pre = 1, sd.post = 1, rho = NULL,
                               alpha = 0.05, power = NULL,
                               sided = c("two", "one"), strict = TRUE) {

  # Check if the arguments are specified correctly
  if (sum(sapply(list(n, delta, power, alpha), is.null)) != 1)
    stop("exactly one of n, delta, alpha, and power must be NULL")
  if (is.null(sd.pre) | is.null(sd.post) | is.null(rho))
    stop("please specify sd.pre, sd.post, and rho")

  # Assign number of sides
  sided <- match.arg(sided)
  side <- switch(sided, one = 1, two = 2)

  # Use absolute value of the effect size
  if (!is.null(delta))
    delta <- abs(delta)

  # Calculate the standard deviation of differences within pairs
  sd.d <- sqrt(sd.pre^2 + sd.post^2 - 2 * rho * sd.pre * sd.post)

  # Calculate test statistic
  p.body <- quote({
    stats::pnorm(stats::qnorm(alpha / side) + sqrt(n) * delta / sd.d)
  })

  if (strict && side == 2)
    p.body <- quote({
      stats::pnorm(stats::qnorm(alpha / side) + sqrt(n) * delta / sd.d) +
        stats::pnorm(stats::qnorm(alpha / side) - sqrt(n) * delta / sd.d)
    })

  # Use uniroot function to calculate missing argument
  if (is.null(power))
    power <- eval(p.body)
  else if (is.null(n))
    n <- uniroot(function(n) eval(p.body) - power, c(2, 1e+07))$root
  else if (is.null(delta))
    delta <- uniroot(function(delta) eval(p.body) - power,  c(1e-07, 1e+07))$root
  else if (is.null(alpha))
    alpha <- uniroot(function(alpha) eval(p.body) - power, c(1e-10, 1 - 1e-10))$root
  else stop("internal error")

  # Generate output text
  METHOD <- "Paired z test power calculation"
  NOTE <- "n is the number of pairs"
  sd <- c(sd.pre, sd.post)

  # Print output as a power.htest object
  structure(list(n = n, delta = delta, sd = sd, rho = rho,
                 alpha = alpha, power = power, sided = sided,
                 method = METHOD, note = NOTE), class = "power.htest")
}
