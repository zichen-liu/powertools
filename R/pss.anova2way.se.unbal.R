#' Power calculations for two-way unbalanced analysis of variance simple effects test
#'
#' @param nmatrix A matrix of sample sizes (see example).
#' @param mmatrix A matrix of group means (see example).
#' @param cmatrix A matrix of contrast coefficients (see example).
#' @param sd The estimated standard deviation within each group.
#' @param Rsq The estimated R^2 for regressing the outcome on the covariates; defaults to 0.
#' @param ncov The number of covariates adjusted for in the model; defaults to 0.
#' @param alpha The significance level or type 1 error rate; defaults to 0.05.
#'
#' @return A list of the arguments (including the computed power).
#' @export
#'
#' @examples
#' nmatrix <- matrix(c(30, 30, 30, 30, 30, 30), nrow = 2, byrow = TRUE)
#' mmatrix <- matrix(c(9.3, 8.9, 8.5, 8.7, 8.3, 7.3), nrow = 2, byrow = TRUE)
#' cmatrix <- matrix(c(-1, 0, 0, 1, 0, 0), nrow = 2, byrow = TRUE)
#' pss.anova2way.se.unbal(nmatrix = nmatrix, mmatrix = mmatrix, cmatrix = cmatrix,
#' sd = 2, alpha = 0.025)

pss.anova2way.se.unbal <- function (nmatrix = NULL, mmatrix = NULL, cmatrix = NULL,
                                   sd = 0, Rsq = 0, ncov = 0, alpha = 0.05) {

  # Check if the arguments are specified correctly
  pss.check(nmatrix, "req"); pss.check(nmatrix, "mat")
  pss.check(mmatrix, "req"); pss.check(mmatrix, "mat")
  pss.check(cmatrix, "req"); pss.check(cmatrix, "mat")
  pss.check(sd, "req"); pss.check(sd, "pos")
  pss.check(Rsq, "req"); pss.check(Rsq, "uniti")
  pss.check(ncov, "req"); pss.check(ncov, "int")
  pss.check(alpha, "req"); pss.check(alpha, "unit")

  a <- nrow(mmatrix)
  b <- ncol(mmatrix)
  if (a != nrow(cmatrix) | b != ncol(cmatrix))
    stop("number of contrast coefficients must be equal to the number of groups")
  if (a != nrow(nmatrix) | b != ncol(nmatrix))
    stop("number of sample sizes must be equal to the number of groups")

  if (any(nmatrix < 2))
    stop("number of observations in each group must be at least 2")

  if (Rsq > 0 & ncov == 0)
    stop("please specify ncov or set Rsq to 0")

  # See if there is an interaction
  fAB <- pss.es.anova.f(means = mmatrix, sd = sd)$fAB
  intx <- ifelse(fAB == 0, FALSE, TRUE)

  # Get lambda
  lambda <- sum(cmatrix * mmatrix) / sd / sqrt(sum(cmatrix^2 / nmatrix)) /
    sqrt(1 - Rsq)

  # Calculate power
  N <- sum(nmatrix)
  df <- ifelse(intx, N - a * b - ncov, N - a - b + 1 - ncov)
  power <- stats::pt(q = stats::qt(alpha, df), df, lambda)

  # Generate output text
  METHOD <- paste0("Unbalanced two-way analysis of ", ifelse(ncov < 1, "", "co"),
                   "variance\n     simple effects power calculation",
                   ifelse(intx, " with interaction", ""))
  out <- list(nmatrix = pss.matrix.format(nmatrix),
              mmatrix = pss.matrix.format(mmatrix),
              cmatrix = pss.matrix.format(cmatrix),
              sd = sd, ncov = ncov, Rsq = Rsq, alpha = alpha, power = power,
              method = METHOD)

  # Print output as a power.htest object
  if (ncov < 1) out <- out[!names(out) %in% c("ncov", "Rsq")]
  structure(out, class = "power.htest")

}

