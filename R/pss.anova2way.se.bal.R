#' Power calculations for two-way balanced analysis of variance simple effects test
#'
#' @param n The sample size per group.
#' @param mmatrix A matrix of group means (see example).
#' @param cmatrix A matrix of contrast coefficients (see example).
#' @param sd The estimated standard deviation within each group; defaults to 1.
#' @param Rsq The estimated R^2 for regressing the outcome on the covariates; defaults to 0.
#' @param ncov The number of covariates adjusted for in the model; defaults to 0.
#' @param alpha The significance level or type 1 error rate; defaults to 0.05.
#' @param power The specified level of power.
#' @param v Either TRUE for verbose output or FALSE to output computed argument only.
#'
#' @return A list of the arguments (including the computed one).
#' @export
#'
#' @examples
#' mmatrix <- matrix(c(9.3, 8.9, 8.5, 8.7, 8.3, 7.3), nrow = 2, byrow = TRUE)
#' cmatrix <- matrix(c(-1, 0, 0, 1, 0, 0), nrow = 2, byrow = TRUE)
#' pss.anova2way.se.bal(n = 30, mmatrix = mmatrix, cmatrix = cmatrix, sd = 2, alpha = 0.025)

pss.anova2way.se.bal <- function (n = NULL, mmatrix = NULL, cmatrix = NULL,
                                 sd = 1, Rsq = 0, ncov = 0,
                                 alpha = 0.05, power = NULL, v = TRUE) {

  # Check if the arguments are specified correctly
  pss.check.many(list(n, alpha, power), "oneof")
  pss.check(n, "int"); pss.check(n, "min", min = 2)
  pss.check(mmatrix, "req"); pss.check(mmatrix, "mat")
  pss.check(cmatrix, "req"); pss.check(cmatrix, "mat")
  pss.check(sd, "req"); pss.check(sd, "pos")
  pss.check(Rsq, "req"); pss.check(Rsq, "uniti")
  pss.check(ncov, "req"); pss.check(ncov, "int")
  pss.check(alpha, "unit")
  pss.check(power, "unit")
  pss.check(v, "req"); pss.check(v, "bool")

  a <- nrow(mmatrix)
  b <- ncol(mmatrix)
  if (a != nrow(cmatrix) | b != ncol(cmatrix))
    stop("number of contrast coefficients must be equal to the number of groups")

   if (Rsq > 0 & ncov == 0)
    stop("please specify ncov or set Rsq to 0")

  # See if there is an interaction
  fAB <- pss.es.anova.f(means = mmatrix, sd = sd)$fAB
  intx <- ifelse(fAB == 0, FALSE, TRUE)

  # Get test statistic
  p.body <- quote({
    lambda <- sum(cmatrix * mmatrix) / sd / sqrt(sum(cmatrix^2 / n)) /
      sqrt(1 - Rsq)
    N <- a * b * n
    df <- ifelse(intx, N - a * b - ncov, N - a - b + 1 - ncov)
    stats::pt(q = stats::qt(alpha, df), df, lambda)
  })

  # Use stats::uniroot function to calculate missing argument
  if (is.null(power)) {
    power <- eval(p.body)
    if (!v) return(power)
  }
  else if (is.null(n)) {
    n <- stats::uniroot(function(n) eval(p.body) - power, c(2, 1e+05))$root
    if (!v) return(n)
  }
  else if (is.null(alpha)) {
    alpha <- stats::uniroot(function(alpha) eval(p.body) - power, c(1e-10, 1 - 1e-10))$root
    if (!v) return(alpha)
  }
  else stop("internal error", domain = NA)

  # Generate output text
  METHOD <- paste0("Balanced two-way analysis of ", ifelse(ncov < 1, "", "co"),
                   "variance\n     simple effects power calculation",
                   ifelse(intx, " with interaction", ""))
  out <- list(n = n, mmatrix = pss.matrix.format(mmatrix),
              cmatrix = pss.matrix.format(cmatrix),
              sd = sd, ncov = ncov, Rsq = Rsq, alpha = alpha, power = power,
              method = METHOD)

  # Print output as a power.htest object
  if (ncov < 1) out <- out[!names(out) %in% c("ncov", "Rsq")]
  structure(out, class = "power.htest")

}

