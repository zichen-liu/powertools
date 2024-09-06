#' Power for detecting treatment effect heterogeneity in a cluster randomized trial with a continuous outcome
#'
#' @description
#' This function performs power and sample size calculations for detecting a treatment-by-covariate
#' interaction effect in a two-arm cluster randomized trial
#' with a continuous, normal outcome when the data will be analyzed using a linear mixed effects model
#' (random intercept for cluster and fixed effect for the treatment-by-covariate interaction).
#'
#'
#' @param m The number of subjects per cluster or the mean cluster size (if unequal number of participants per cluster).
#' @param m.sd The standard deviation of cluster sizes (provide if unequal number of participants per cluster); defaults to 0.
#' @param J1 The number of clusters in arm 1.
#' @param J.ratio The ratio J2/J1 between the number of clusters in the two arms; defaults to 1 (equal clusters per arm).
#' @param delta The difference between the intervention and control means under the alternative minus the difference under the null hypothesis.
#' @param sd The total standard deviation of the outcome variable; defaults to 1.
#' @param icc1 The intraclass correlation coefficient in arm 1; defaults to 0.
#' @param alpha The significance level (type 1 error rate); defaults to 0.05.
#' @param power The specified level of power.
#' @param sides Either 1 or 2 (default) to specify a one- or two- sided hypothesis test.
#' @param v Either TRUE for verbose output or FALSE (default) to output computed argument only.
#'
#' @return A list of the arguments (including the computed one).
#' @export
#'
#' @examples
#' crt.parallel.hte(m = 30, J1 = 8, delta = 0.4, icc1 = 0.05, icc2 = 0.05)


crt.parallel.hte <- function (m = NULL, m.sd = 0, J1 = NULL, J.ratio = 1, delta = NULL, sd = 1,
                               icc1 = 0, icc2 = 0,
                               alpha = 0.05, power = NULL, sides = 2,
                               v = FALSE) {



  # Calculate power
  if (sides == 1)
    p.body <- quote({
      J <- J1 + J1 * J.ratio
      N <- m * J
      df <- J - 2 - ncov
      d <- delta / sd

      RE <- re.clustsize.cont(m = m, m.sd = m.sd, icc = (icc1 + icc2)/2)

      w <- 1 / (1 + J.ratio)
      de1 <- 1 + (m * (1 - RsqB) - 1) * icc1 -
             m * RsqW * (1 - 2 * icc1) / (m - 1)
      de2 <- 1 + (m * (1 - RsqB) - 1) * icc2 -
             m * RsqW * (1 - 2 * icc2) / (m - 1)
      detot <- de1 / w + de2 / (1 - w)
      ncp <- d / sqrt(detot / N / RE)
      crit <- stats::qt(1 - alpha, df)
      1 - stats::pt(crit, df, ncp)
    })
  else if (sides == 2)
    p.body <- quote({
      J <- J1 + J1 * J.ratio
      N <- m * J
      df2 <- J - 2 - ncov
      d <- delta / sd

      RE <- re.clustsize.cont(m = m, m.sd = m.sd, icc = (icc1 + icc2)/2)

      w <- 1 / (1 + J.ratio)
      de1 <- 1 + (m * (1 - RsqB) - 1) * icc1 -
        m * RsqW * (1 - 2 * icc1) / (m - 1)
      de2 <- 1 + (m * (1 - RsqB) - 1) * icc2 -
        m * RsqW * (1 - 2 * icc2) / (m - 1)
      detot <- de1 / w + de2 / (1 - w)
      ncp <- d / sqrt(detot / N / RE)
      crit <- stats::qf(1 - alpha, 1, df2)
      1 - stats::pf(crit, 1, df2, ncp^2)
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
  else if (is.null(J.ratio)) {
    J.ratio <- stats::uniroot(function(J.ratio) eval(p.body) - power, c(1e-10, 1e+07))$root
    if (!v) return(J.ratio)
  }
  else if (is.null(m)) {
    m <- stats::uniroot(function(m) eval(p.body) - power, c(m.sd/2 + 2, 1e+07))$root
    if (!v) return(m)
  }
  else if (is.null(delta)) {
    delta <- stats::uniroot(function(delta) eval(p.body) - power, c(1e-07, 1e+07))$root
    if (!v) return(delta)
  }

  # Generate output text
  METHOD <- "Power for treatment-by-covariate interaction in a cluster randomized trial with a continuous outcome"
  m <- ifelse(m.sd == 0, m, paste0(m, " (sd = ", m.sd, ")"))
  J <- c(J1, J1 * J.ratio)
  icc <- c(icc1, icc2)
  out <- list(m = m, `J1, J2` = J, delta = delta, sd = sd, `icc1, icc2` = icc,
              alpha = alpha, power = power, sides = sides, method = METHOD)

  # Print output as a power.htest object

  structure(out, class = "power.htest")

}

