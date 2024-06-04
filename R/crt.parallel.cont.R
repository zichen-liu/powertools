#' Power for a cluster randomized trial with a continuous outcome
#'
#' @param m The number of subjects per cluster or the mean cluster size (if unequal number of participants per cluster).
#' @param m.sd The standard deviation of cluster sizes (provide if unequal number of participants per cluster); defaults to 0.
#' @param J1 The number of clusters in arm 1.
#' @param J.ratio The ratio J2/J1 between the number of clusters in the two arms; defaults to 1 (equal clusters per arm).
#' @param delta The difference between the intervention and control means under the alternative minus the difference under the null hypothesis.
#' @param sd The total standard deviation of the outcome variable; defaults to 1.
#' @param icc1 The intraclass correlation coefficient in arm 1; defaults to 0.
#' @param icc2 The intraclass correlation coefficient in arm 2; defaults to 0.
#' @param RsqB The estimated proportion of total variance explained by cluster-level covariates; defaults to 0.
#' @param RsqW The estimated proportion of total variance explained by individual-level covariates; defaults to 0.
#' @param ncov The number of cluster-level and individual-level covariates; defaults to 0.
#' @param alpha The significance level or type 1 error rate; defaults to 0.05.
#' @param power The specified level of power.
#' @param sides Either 1 or 2 (default) to specify a one- or two- sided hypothesis test.
#' @param v Either TRUE for verbose output or FALSE to output computed argument only.
#'
#' @return A list of the arguments (including the computed one).
#' @export
#'
#' @examples
#' crt.parallel.cont(m = 30, J1 = 8, delta = 0.4, icc1 = 0.05, icc2 = 0.05)
#' crt.parallel.cont(m = NULL, J1 = 6, delta = 0.5, icc1 = 0.05, icc2 = 0.05, power = 0.8)
#' crt.parallel.cont(m = 25, m.sd = 15, J1 = NULL, delta = 0.3, icc1 = 0.05,
#' icc2 = 0.05, power = 0.8)
#' crt.parallel.cont(m = 20, J1 = 15, delta = 0.3, icc1 = 0.05, icc2 = 0.05,
#' RsqB = 0.1, ncov = 1, sides = 1)
#' crt.parallel.cont(m = 20, J1 = 15, delta = 0.3, icc1 = 0.05, icc2 = 0.05,
#' RsqW = 0.5^2, ncov = 1, sides = 1)
#' crt.parallel.cont(m = NULL, J1 = 6, delta = 0.5, icc1 = 0.05, icc2 = 0.05, power = 0.8)

crt.parallel.cont <- function (m = NULL, m.sd = 0, J1 = NULL, J.ratio = 1, delta = NULL, sd = 1,
                                   icc1 = 0, icc2 = 0, ncov = 0, RsqB = 0, RsqW = 0,
                                   alpha = 0.05, power = NULL, sides = 2,
                                   v = TRUE) {

  # Check if the arguments are specified correctly
  check.many(list(m, J1, delta, alpha, power), "oneof")
  check(m, "int")
  check(m.sd, "req"); check(m.sd, "min", min = 0)
  check(J.ratio, "req"); check(J.ratio, "pos")
  check(J1, "min", min = 2)
  if (!is.null(J1))
    check(J1 * J.ratio, "min", min = 2)
  check(delta, "num")
  check(sd, "req"); check(sd, "pos")
  check(icc1, "req"); check(icc1, "uniti")
  check(icc2, "req"); check(icc2, "uniti")
  check(ncov, "req"); check(ncov, "int")
  check(RsqB, "req"); check(RsqB, "uniti")
  check(RsqW, "req"); check(RsqW, "uniti")
  check(alpha, "unit")
  check(power, "unit")
  check(sides, "req"); check(sides, "vals", valslist = c(1, 2))
  check(v, "req"); check(v, "bool")

  if ((RsqB > 0 | RsqW > 0) & ncov == 0)
    stop("please specify ncov or set RsqB & RsqW to 0")

  # feasability check J > rho * Nindep
  if (is.null(m)) {
    out <- t.test.2samp(n1 = NULL, delta = delta, sd1 = sd, power = power, sides = sides)
    Nindep <- out[[1]][1] + out[[1]][2]
    iccavg <- mean(icc1, icc2)
    if (J1 + J1 * J.ratio <= iccavg * Nindep)
      stop("desired power unable to be achieved with given conditions")
  }

  # Calculate power
  p.body <- quote({
    J <- J1 + J1 * J.ratio
    N <- m * J
    df <- J - 2 - ncov
    d <- delta / sd

    RE <- varym.re(m = m, m.sd = m.sd, icc = (icc1 + icc2)/2)

    w <- 1 / (1 + J.ratio)
    de1 <- 1 + (m - 1) * (1 - RsqB) * icc1 -
           m * RsqW * (1 - 2 * icc1) / (m - 1)
    de2 <- 1 + (m - 1) * (1 - RsqB) * icc2 -
           m * RsqW * (1 - 2 * icc2) / (m - 1)
    detot <- de1 / w + de2 / (1 - w)
    ncp <- d / sqrt(detot / N / RE)
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
  METHOD <- "Power for a cluster randomized trial with a continuous outcome"
  m <- ifelse(m.sd == 0, m, paste0(m, " (sd = ", m.sd, ")"))
  J <- c(J1, J1 * J.ratio)
  icc <- c(icc1, icc2)
  Rsq <- c(RsqB, RsqW)
  out <- list(m = m, `J1, J2` = J, delta = delta, sd = sd, `icc1, icc2` = icc,
              ncov = ncov, `RsqB, RsqW` = Rsq,
              alpha = alpha, power = power, sides = sides, method = METHOD)

  # Print output as a power.htest object
  if (RsqB < 0.0000000001 & RsqW < 0.0000000001)
    out <- out[!names(out) %in% c("ncov", "RsqB, RsqW")]
  structure(out, class = "power.htest")

}

