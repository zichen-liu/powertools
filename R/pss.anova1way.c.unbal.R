#' Power calculations for one-way unbalanced analysis of variance contrast test
#'
#' @param nvec A vector of group sample sizes c(n1, n2, ...).
#' @param mvec A vector of group mvec c(mu1, mu2, ...).
#' @param cvec A vector of contrast coefficients c(c1, c2, ...).
#' @param sd The estimated standard deviation within each group.
#' @param Rsq The estimated R^2 for regressing the outcome on the covariates; defaults to 0.
#' @param ncov The number of covariates adjusted for in the model; defaults to 0.
#' @param alpha The significance level or type 1 error rate; defaults to 0.05.
#' @param v Either TRUE for verbose output or FALSE to output computed argument only.
#'
#' @return A list of the arguments (including the computed power).
#' @export
#'
#' @examples
#' pss.anova1way.c.unbal(nvec = c(20, 20, 20), mvec = c(5, 10, 12), cvec = c(1, -1, 0),
#' sd = 10, alpha = 0.025)
#' pss.anova1way.c.unbal(nvec = c(20, 20, 20), mvec = c(5, 10, 12), cvec = c(1, 0, -1),
#' sd = 10, alpha = 0.025)

pss.anova1way.c.unbal <- function (nvec = NULL, mvec = NULL, cvec = NULL,
                                   sd = NULL, Rsq = 0, ncov = 0, alpha = 0.05,
                                   v = TRUE) {

  # Check if the arguments are specified correctly
  pss.check(nvec, "req"); pss.check(nvec, "vec")
  pss.check(mvec, "req"); pss.check(mvec, "vec")
  pss.check(cvec, "req"); pss.check(cvec, "vec")
  pss.check(sd, "req"); pss.check(sd, "pos")
  pss.check(Rsq, "req"); pss.check(Rsq, "uniti")
  pss.check(ncov, "req"); pss.check(ncov, "int")
  pss.check(alpha, "req"); pss.check(alpha, "unit")
  pss.check(v, "req"); pss.check(v, "bool")

  a <- length(mvec)
  if (a != length(nvec))
    stop("number of sample sizes must equal to the number of groups")
  if (a != length(cvec))
    stop("number of contrast coefficients must be equal to the number of groups")

  if (any(nvec < 2))
    stop("number of observations in each group must be at least 2")

  if (Rsq > 0 & ncov == 0)
    stop("please specify ncov or set Rsq to 0")

  # Get lambda (Lambda = lambda^2)
  num <- cvec %*% mvec
  temp <- sapply(X = 1:a, FUN = function(i) cvec[i]^2 / nvec[i])
  den <- sd * sqrt(sum(temp))
  lambda <- num / den / sqrt(1 - Rsq)

  # Calculate power
  df2 <- sum(nvec) - a - ncov
  power <- stats::pf(stats::qf(alpha, 1, df2, lower.tail = FALSE),
                     1, df2, lambda^2, lower.tail = FALSE)
  if (!v) return(power)

  # Generate output text
  METHOD <- paste0("Unbalanced one-way analysis of ", ifelse(ncov < 1, "", "co"),
                   "variance\n     contrast test power calculation")
  out <- list(nvec = nvec, mvec = mvec, cvec = cvec,
              sd = sd, ncov = ncov, Rsq = Rsq, alpha = alpha, power = power,
              method = METHOD)

  # Print output as a power.htest object
  if (ncov < 1) out <- out[!names(out) %in% c("ncov", "Rsq")]
  structure(out, class = "power.htest")

}
