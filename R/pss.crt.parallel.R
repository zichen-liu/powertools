#' Power for test of treatment effect in cluster randomized trials
#'
#' @param m The number of subjects per cluster or the mean cluster size (if unequal number of participants per cluster).
#' @param m.sd The standard deviation of cluster sizes (provide if unequal number of participants per cluster); defaults to 0.
#' @param w The proportion of clusters allocated to arm 1; defaults to 0.5.
#' @param J The number of clusters.
#' @param delta The difference between the intervention and control means in the outcome variable.
#' @param sd The total standard deviation of the outcome variable; defaults to 1.
#' @param rho1 The intraclass correlation coefficient in arm 1; defaults to 0.
#' @param rho2 The intraclass correlation coefficient in arm 2; defaults to 0.
#' @param rhoBsq The estimated proportion of total variance explained by cluster-level covariates; defaults to 0.
#' @param rhoWsq The estimated proportion of total variance explained by variation of the covariate within the cluster; defaults to 0.
#' @param q The number of covariates; defaults to 0.
#' @param alpha The significance level or type 1 error rate; defaults to 0.05.
#' @param power The specified level of power.
#' @param sides Either 1 or 2 (default) to specify a one- or two- sided hypothesis test.
#' @return A list of the arguments (including the computed one).
#' @export
#'
#' @examples
#' pss.crt.parallel(m = 30, J = 16, delta = 0.4, rho1 = 0.05, rho2 = 0.05)
#' pss.crt.parallel(m = NULL, J = 12, delta = 0.5, rho1 = 0.05, rho2 = 0.05, power = 0.8)
#' pss.crt.parallel(m = 25, m.sd = 15, J = NULL, delta = 0.3, rho1 = 0.05, rho2 = 0.05, power = 0.8)

pss.crt.parallel <- function (m = NULL, m.sd = 0, w = 0.5, J = NULL, delta = NULL, sd = 1,
                              rho1 = 0, rho2 = 0, Rsq = 0, q = 0, alpha = 0.05, power = NULL, sides = 2) {

  # Calculate power
  p.body <- quote({
    N <- m * J
    df <- J - 2 - q
    d <- delta / sd

    cv <- m.sd / m
    K <- (m * rho1) / (1 + (m - 1) * rho1)
    RE <- 1 - cv^2 * K * (1 - K)

    # where do the rhoBsq & rhoWsq go? do the q's add up?
    de <- (((1 + (m - 1) * rho1) / w) + ((1 + (m - 1) * rho2) / (1 - w)))
    ncp <- d / sqrt(de / N / RE) # does RE go in the sqrt for multisite too?
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
  else if (is.null(delta))
    delta <- stats::uniroot(function(delta) eval(p.body) - power, c(1e-07, 1e+07))$root

  # Generate output text
  METHOD <- "Power for test of treatment effect in cluster randomized trials"
  NOTE <- "m is the subjects per site"

  # Print output as a power.htest object
  structure(list(m = m, m.sd = m.sd, J = J, delta = delta, sd = sd,
                 rho1 = rho1, Rsq = Rsq, q = q, alpha = alpha, power = power, sides = sides,
                 method = METHOD, note = NOTE), class = "power.htest")

}

