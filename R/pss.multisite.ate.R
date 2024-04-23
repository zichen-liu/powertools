#' Power for test of average treatment effect in a multisite trial
#'
#' @param m The number of subjects per site or the mean cluster size (if unequal number of participants per site).
#' @param m.sd The standard deviation of cluster sizes (provide if unequal number of participants per site); defaults to 0.
#' @param alloc.ratio The allocation ratio of intervention/control per site; defaults to 1.
#' @param J The number of sites.
#' @param delta The difference between the intervention and control means in the outcome variable.
#' @param sd The total standard deviation of the outcome variable; defaults to 1.
#' @param icc0 The proportion of total variance of the outcome attributable to variation in site-level means.
#' @param icc1 The proportion of total variance of the outcome attributable to variation in the treatment effect across sites.
#' @param Rsq The estimated R^2 for regressing the outcome on the covariates; defaults to 0.
#' @param alpha The significance level or type 1 error rate; defaults to 0.05.
#' @param power The specified level of power.
#' @param sides Either 1 or 2 (default) to specify a one- or two- sided hypothesis test.
#' @return A list of the arguments (including the computed one).
#' @export
#'
#' @examples
#' pss.multisite.ate(m = 20, J = 10, delta = 3, sd = sqrt(40), icc0 = 0.1, icc1 = 0)
#' pss.multisite.ate(m = 20, J = 10, delta = 3, sd = sqrt(48), icc0 = 0.095, icc1 = 0.048)
#' pss.multisite.ate(m = 20, alloc.ratio = 1.5, J = 10, delta = 0.43, icc0 = 0.095, icc1 = 0.048)
#' pss.multisite.ate(m = 10, J = NULL, delta = 0.5, sd = 1, icc0 = 0, icc1 = 0.05, power = 0.8)
#' pss.multisite.ate(m = 20, m.sd = 5, J = 10, delta = 3, sd = sqrt(48), icc0 = 0.095, icc1 = 0.048)
#' pss.multisite.ate(m = 20, J = 10, delta = 3, sd = sqrt(48), icc0 = 0.095, icc1 = 0.048, Rsq = 0.5^2)

pss.multisite.ate <- function (m = NULL, m.sd = 0, alloc.ratio = 1, J = NULL,
                               delta = NULL, sd = 1,
                               icc0 = NULL, icc1 = NULL, Rsq = 0,
                               alpha = 0.05, power = NULL, sides = 2) {

  # Check if the arguments are specified correctly
  pss.check.many(list(m, J, delta, alpha, power), "oneof")
  pss.check(m, "int")
  pss.check(m.sd, "req"); pss.check(m.sd, "min", min = 0)
  pss.check(alloc.ratio, "req"); pss.check(alloc.ratio, "pos")
  pss.check(J, "min", min = 2)
  pss.check(delta, "num")
  pss.check(sd, "req"); pss.check(sd, "pos")
  pss.check(icc0, "req"); pss.check(icc0, "uniti")
  pss.check(icc1, "req"); pss.check(icc1, "uniti")
  pss.check(Rsq, "req"); pss.check(Rsq, "uniti")
  pss.check(alpha, "unit")
  pss.check(power, "unit")
  pss.check(sides, "req"); pss.check(sides, "vals", valslist = c(1, 2))

  # Calculate power
  p.body <- quote({
    N <- m * J
    df <- J - 1
    d <- delta / (sd * sqrt((1 - Rsq)))

    RE <- pss.re(m = m, m.sd = m.sd, icc = icc1)

    c <- (1 + alloc.ratio)^2 / alloc.ratio
    ncp <- d / sqrt(c * (1 - icc0 + (4 * m / c - 1) * icc1) / N / RE)
    crit <- stats::qt(1 - alpha / sides, df)
    1 - stats::pt(crit, df, ncp)
  })

  # Use uniroot to calculate missing argument
  if (is.null(alpha))
    alpha <- stats::uniroot(function(alpha) eval(p.body) - power, c(1e-10, 1 - 1e-10))$root
  else if (is.null(power))
    power <- eval(p.body)
  else if (is.null(J))
    J <- stats::uniroot(function(J) eval(p.body) - power, c(2 + 1e-10, 1e+07))$root
  else if (is.null(m))
    m <- stats::uniroot(function(m) eval(p.body) - power, c(2 + 1e-10, 1e+07))$root
  # else if (is.null(alloc.ratio)) {
  #  alloc.ratio <- stats::uniroot(function(alloc.ratio) eval(p.body) - power, c(1 + 1e-10, 1e+07))$root
  # }
  else if (is.null(delta))
    delta <- stats::uniroot(function(delta) eval(p.body) - power, c(1e-07, 1e+07))$root

  # Generate output text
  METHOD <- "Power for test of average treatment effect in multisite trials"
  NOTE <- "NOTE: m1, m2 are the number of subjects within site in\ncondition 1, condition 2 (total of m1 + m2 per site)"
  icc <- c(icc0, icc1)
  c <- m / (alloc.ratio + 1)
  t <- alloc.ratio * c
  m <- ifelse(m.sd == 0, paste0(t, ", ", c), paste0(t, ", ",  c, " (sd = ", m.sd, ")"))

  out <- list(`m1, m2` = m, J = J, delta = delta, sd = sd,
              `icc0, icc1` = icc, Rsq = Rsq,
              alpha = alpha, power = power, sides = sides,
              method = METHOD, note = NOTE)

  # Print output as a power.htest object
  if (Rsq < 0.00000000001) out <- out[!names(out) %in% c("Rsq")]
  structure(out, class = "power.htest")

}

