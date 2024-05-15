#' Power for test of treatment effect in 2x2 crossover cluster randomized trial
#'
#' @param m The number of subjects measured during each cluster-period.
#' @param J.arm The number of clusters in each arm.
#' @param delta The difference between the intervention and control means under the alternative minus the difference under the null hypothesis.
#' @param sd The total standard deviation of the outcome variable; defaults to 1.
#' @param icc The within-cluster, within-period intraclass correlation coefficient; defaults to 0.
#' @param iccb The within-cluster, between-period intraclass correlation coefficient. Either iccb and xi, OR cac and sac must be specified.
#' @param xi The  within-cluster, within-subject correlation (correlation between two measurements within the same subject). Either iccb and xi, OR cac and sac must be specified.
#' @param cac The cluster autocorrelation. Either iccb and xi, OR cac and sac must be specified.
#' @param sac The subject autocorrelation. Either iccb and xi, OR cac and sac must be specified.
#' @param alpha The significance level or type 1 error rate; defaults to 0.05.
#' @param power The specified level of power.
#' @param sides Either 1 or 2 (default) to specify a one- or two- sided hypothesis test.
#' @return A list of the arguments (including the computed one).
#' @export
#'
#' @examples
#' pss.crt.xo.cont(m = 30, J.arm = 4, delta = 0.3, icc = 0.05, cac = 0.8, sac = 0.4)
#' pss.crt.xo.cont(m = 30, J.arm = 4, delta = 0.3, icc = 0.05, iccb = 0.04, xi = 0.42)


pss.crt.xo.cont <- function (m = NULL, J.arm = NULL, delta = NULL, sd = 1,
                             icc = 0, iccb = NULL, xi = NULL, cac = NULL, sac = NULL,
                             alpha = 0.05, power = NULL, sides = 2) {

  # Check if the arguments are specified correctly
  if ((is.null(iccb) | is.null(xi)) & (is.null(cac) | is.null(sac)))
    stop("iccb and xi, OR cac and sac must be specified")
  pss.check.many(list(m, J.arm, delta, alpha, power), "oneof")
  pss.check(m, "int")
  pss.check(J.arm, "min", min = 2)
  pss.check(delta, "num")
  pss.check(sd, "req"); pss.check(sd, "pos")
  pss.check(icc, "req"); pss.check(icc, "uniti")
  pss.check(iccb, "uniti")
  pss.check(xi, "uniti")
  pss.check(cac, "uniti")
  pss.check(sac, "uniti")
  pss.check(alpha, "unit")
  pss.check(power, "unit")
  pss.check(sides, "req"); pss.check(sides, "vals", valslist = c(1, 2))

  if (is.null(iccb) & is.null(xi)) {
    iccb <- icc * cac
    xi <- ifelse(sac == 0, 0, (1 - icc) * sac + icc * cac)
  }

  # Calculate power
  p.body <- quote({
    J <- 2 * J.arm
    N <- m * J
    df <- 2 * J - 3
    d <- delta / sd
    de <- ifelse(sac == 0, 2 * (1 + (m - 1) * icc - m * iccb) , 2 * (1 - xi + (m - 1) * (icc - iccb)))
    ncp <- d / sqrt(de / N)
    crit <- stats::qt(1 - alpha / sides, df)
    1 - stats::pt(crit, df, ncp)
  })

  # Use uniroot to calculate missing argument
  if (is.null(alpha))
    alpha <- stats::uniroot(function(alpha) eval(p.body) - power, c(1e-10, 1 - 1e-10))$root
  else if (is.null(power))
    power <- eval(p.body)
  else if (is.null(J.arm))
    J.arm <- stats::uniroot(function(J.arm) eval(p.body) - power, c(2 + 1e-10, 1e+07))$root
  else if (is.null(m))
    m <- stats::uniroot(function(m) eval(p.body) - power, c(2 + 1e-10, 1e+07))$root
  else if (is.null(delta))
    delta <- stats::uniroot(function(delta) eval(p.body) - power, c(1e-07, 1e+07))$root

  # Generate output text
  METHOD <- "Power for test of treatment effect in a 2x2 cluster randomized crossover trial"
  m <- m
  J <- c(J.arm, J.arm)

  # Print output as a power.htest object depending on which inputs were given
  if (!is.null(cac) & !is.null(sac)) {
    acs <- c(cac, sac)
    out <- list(m = m, `J.arm1, J.arm2` = J, delta = delta, sd = sd, icc = icc, `CAC, SAC` = acs,
                alpha = alpha, power = power, sides = sides, method = METHOD)
    structure(out, class = "power.htest")
  } else {
    iccs <- c(icc, iccb, xi)
    out <- list(m = m, `J.arm1, J.arm2` = J, delta = delta, sd = sd, `icc, iccb, xi` = iccs,
                alpha = alpha, power = power, sides = sides, method = METHOD)
    structure(out, class = "power.htest")
  }
}



