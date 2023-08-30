#' Power calculations for one-way balanced analysis of variance omnibus F test
#'
#' @param n The sample size per group.
#' @param mvec A vector of group means c(mu1, mu2, ...).
#' @param sd The estimated standard deviation within each group; defaults to 1.
#' @param rho The estimated correlation between covariates and the outcome; defaults to 0.
#' @param ncov The number of covariates adjusted for in the model; defaults to 0.
#' @param alpha The significance level or type 1 error rate; defaults to 0.05.
#' @param power The specified level of power.
#'
#' @return A list of the arguments (including the computed one).
#' @export
#'
#' @examples
#' # Example 5.2
#' pss.anova.1w(n = 20, mvec = c(5, 10, 12), sd = 10)

pss.anova.1w <- function (n = NULL, mvec = NULL, sd = 1, rho = 0, ncov = 0,
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

  # Get f effect size
  f <- pss.effect.size(means = mvec, sd = sd)$fA

  # Calculate df's and ncp
  p.body <- quote({
    Lambda <- a * n * f^2 / (1 - rho^2)
    df1 <- a - 1
    df2 <- (n - 1) * a - ncov
    stats::pf(stats::qf(alpha, df1, df2, lower.tail = FALSE),
              df1, df2, Lambda, lower.tail = FALSE)
  })

  # Use uniroot function to calculate missing argument
  if (is.null(power))
    power <- eval(p.body)
  else if (is.null(n))
    n <- uniroot(function(n) eval(p.body) - power, c(2, 1e+05))$root
  else if (is.null(alpha))
    alpha <- uniroot(function(alpha) eval(p.body) - power, c(1e-10, 1 - 1e-10))$root
  else stop("internal error", domain = NA)

  # Generate output text
  METHOD <- paste0("Balanced one-way analysis of ", ifelse(ncov < 1, "", "co"),
                   "variance\n     omnibus F test power calculation")
  out <- list(a = a, mvec = mvec, n = n, sd = sd,
              f = f, ncov = ncov, rho = rho,
              alpha = alpha, power = power,
              method = METHOD)

  # Print output as a power.htest object
  if (ncov < 1) out <- out[!names(out) %in% c("ncov", "rho")]
  structure(out, class = "power.htest")

}
