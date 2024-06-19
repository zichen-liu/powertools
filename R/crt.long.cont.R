#' Power for test of treatment effect in longitudinal cluster randomized trial with baseline measurement
#'
#' @param m The number of subjects measured during each cluster-period.
#' @param J1 The number of clusters in arm 1.
#' @param J.ratio The ratio J2/J1 between the number of clusters in the two arms; defaults to 1 (equal clusters per arm).
#' @param delta The difference between the intervention and control means under the alternative minus the difference under the null hypothesis.
#' @param sd The total standard deviation of the outcome variable; defaults to 1.
#' @param icc The within-cluster, within-period intraclass correlation coefficient; defaults to 0.
#' @param cac The cluster autocorrelation; defaults to 0.
#' @param sac The subject autocorrelation; defaults to 0.
#' @param alpha The significance level or type 1 error rate; defaults to 0.05.
#' @param power The specified level of power.
#' @param sides Either 1 or 2 (default) to specify a one- or two- sided hypothesis test.
#' @param v Either TRUE for verbose output or FALSE to output computed argument only.
#' @return A list of the arguments (including the computed one).
#' @export
#'
#' @examples
#' crt.long.cont(m = 30, J1 = 8, delta = 0.3, icc = 0.05, cac = 0.4, sac = 0.5)


crt.long.cont <- function (m = NULL, J1 = NULL, J.ratio = 1, delta = NULL, sd = 1,
                           icc = 0, cac = 0, sac = 0,
                           alpha = 0.05, power = NULL, sides = 2,
                           v = FALSE) {

  # Check if the arguments are specified correctly
  check.many(list(m, J1, delta, alpha, power), "oneof")
  check(m, "pos")
  check(J.ratio, "req"); check(J.ratio, "pos")
  check(J1, "min", min = 2)
  if (!is.null(J1))
    check(J1 * J.ratio, "min", min = 2)
  check(delta, "num")
  check(sd, "req"); check(sd, "pos")
  check(icc, "req"); check(icc, "uniti")
  check(cac, "req"); check(cac, "uniti")
  check(sac, "req"); check(sac, "uniti")
  check(alpha, "unit")
  check(power, "unit")
  check(sides, "req"); check(sides, "vals", valslist = c(1, 2))
  check(v, "req"); check(v, "bool")

  # Calculate power
  p.body <- quote({
    J <- J1 + J1 * J.ratio
    N <- m * J
    df <- J - 2
    d <- delta / sd
    de.pa <- 1 + (m - 1) * icc
    r <- m * icc * cac / de.pa + (1 - icc) * sac / de.pa
    de <- 4 * (1 - r^2) * de.pa
    ncp <- d / sqrt(de / N)
    crit <- stats::qt(1 - alpha / sides, df)
    1 - stats::pt(crit, df, ncp)
  })

  # Use uniroot to calculate missing argument
  if (is.null(alpha)) {
    alpha <- stats::uniroot(function(alpha) eval(p.body) - power, c(1e-10, 1 - 1e-10))$root
    if (!v) return(alpha)
  }
  else if (is.null(power)) {
    power <- eval(p.body)
    if (!v) return(power)
  }
  else if (is.null(J1)) {
    J1 <- stats::uniroot(function(J1) eval(p.body) - power, c(2 + 1e-10, 1e+07))$root
    if (!v) return(J1)
  }
  else if (is.null(m)) {
    m <- stats::uniroot(function(m) eval(p.body) - power, c(2 + 1e-10, 1e+07))$root
    if (!v) return(m)
  }
  else if (is.null(delta)) {
    delta <- stats::uniroot(function(delta) eval(p.body) - power, c(1e-07, 1e+07))$root
    if (!v) return(delta)
  }

  # Generate output text
  METHOD <- "Power for test of treatment effect in a longitudinal cluster randomized trial with baseline measurement"
  J <- c(J1, J1 * J.ratio)
  de.pa <- 1 + (m - 1) * icc
  r <- m * icc * cac / de.pa + (1 - icc) * sac / de.pa
  ccs <- c(icc, cac, sac)
  out <- list(m = m, `J1, J2` = J, delta = delta, sd = sd, `icc, cac, sac` = ccs, r = r,
              alpha = alpha, power = power, sides = sides, method = METHOD)

  # Print output as a power.htest object
  structure(out, class = "power.htest")

}
