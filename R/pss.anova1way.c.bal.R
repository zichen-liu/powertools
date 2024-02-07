#' Power calculations for one-way balanced analysis of variance contrast test
#'
#' @param n The sample size per group.
#' @param mvec A vector of group means c(mu1, mu2, ...).
#' @param cvec A vector of contrast coefficients c(c1, c2, ...).
#' @param sd The estimated standard deviation within each group; defaults to 1.
#' @param Rsq The estimated R^2 for regressing the outcome on the covariates; defaults to 0.
#' @param ncov The number of covariates adjusted for in the model; defaults to 0.
#' @param alpha The significance level or type 1 error rate; defaults to 0.05.
#' @param power The specified level of power.
#'
#' @return A list of the arguments (including the computed one).
#' @export
#'
#' @examples
#' pss.anova1way.c.bal(n = 20, mvec = c(5, 10, 12), cvec = c(1, -1, 0), sd = 10, alpha = 0.025)
#' pss.anova1way.c.bal(n = 20, mvec = c(5, 10, 12), cvec = c(1, 0, -1), sd = 10, alpha = 0.025)

pss.anova1way.c.bal <- function (n = NULL, mvec = NULL, cvec = NULL, sd = 1,
                                Rsq = 0, ncov = 0, alpha = 0.05, power = NULL) {

  # Check if the arguments are specified correctly
  a <- length(mvec)
  if (sum(vapply(list(n, alpha, power), is.null, NA)) != 1)
    stop("exactly one of 'n', 'alpha', and 'power' must be NULL")
  if (a < 2)
    stop("number of groups must be at least 2")
  if (!is.null(n) && n < 2)
    stop("number of observations in each group must be at least 2")
  if (a != length(cvec))
    stop("number of contrast coefficients must be equal to the number of groups")
  if(is.null(sd))
    stop("sd must be specified")

  # Calculate df and ncp
  p.body <- quote({
    lambda <- cvec %*% mvec / sd / sqrt(sum(cvec^2 / n)) / sqrt(1 - Rsq)
    df2 <- n * a - a - ncov
    stats::pf(q = stats::qf(alpha, 1, df2, lower.tail = FALSE),
              1, df2, lambda^2, lower.tail = FALSE)
  })

  # Use stats::uniroot function to calculate missing argument
  if (is.null(power))
    power <- eval(p.body)
  else if (is.null(n))
    n <- stats::uniroot(function(n) eval(p.body) - power, c(2, 1e+05))$root
  else if (is.null(alpha))
    alpha <- stats::uniroot(function(alpha) eval(p.body) - power, c(1e-10, 1 - 1e-10))$root
  else stop("internal error", domain = NA)

  # Generate output text
  METHOD <- paste0("Balanced one-way analysis of ", ifelse(ncov < 1, "", "co"),
                   "variance\n     contrast test power calculation")
  out <- list(a = a, mvec = mvec, cvec = cvec, n = n, sd = sd,
              ncov = ncov, Rsq = Rsq, alpha = alpha, power = power,
              method = METHOD)

  # Print output as a power.htest object
  if (ncov < 1) out <- out[!names(out) %in% c("ncov", "Rsq")]
  structure(out, class = "power.htest")

}
