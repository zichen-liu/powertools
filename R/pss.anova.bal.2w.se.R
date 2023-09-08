#' Power calculations for two-way balanced analysis of variance simple effects test
#'
#' @param n The sample size per group.
#' @param mmatrix A matrix of group means (see example).
#' @param cmatrix A matrix of contrast coefficients (see example).
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
#' # Example 5.12
#' mmatrix <- matrix(c(9.3, 8.9, 8.5, 8.7, 8.3, 7.3), nrow = 2, byrow = TRUE)
#' cmatrix <- matrix(c(-1, 0, 0, 1, 0, 0), nrow = 2, byrow = TRUE)
#' pss.anova.bal.2w.se(n = 30, mmatrix = mmatrix, cmatrix = cmatrix, sd = 2, alpha = 0.025)

pss.anova.bal.2w.se <- function (n = NULL, mmatrix = NULL, cmatrix = NULL,
                                 sd = 1, rho = 0, ncov = 0,
                                 alpha = 0.05, power = NULL) {

  # Check if the arguments are specified correctly
  a <- nrow(mmatrix)
  b <- ncol(mmatrix)
  if (sum(vapply(list(n, alpha, power), is.null, NA)) != 1)
    stop("exactly one of 'n', 'alpha', and 'power' must be NULL")
  if (a < 2 | b < 2)
    stop("number of groups per intervention must be at least 2")
  if (!is.null(n) && n < 2)
    stop("number of observations in each group must be at least 2")
  if (a != nrow(cmatrix) | b != ncol(cmatrix))
    stop("number of contrast coefficients must be equal to the number of groups")
  if(is.null(sd))
    stop("sd must be specified")

  # See if there is an interaction
  fAB <- pss.effect.size(means = mmatrix, sd = sd)$fAB
  intx <- ifelse(fAB == 0, FALSE, TRUE)

  # Get test statistic
  p.body <- quote({
    lambda <- sum(cmatrix * mmatrix) / sd / sqrt(sum(cmatrix^2 / n)) /
      sqrt(1 - rho^2)
    N <- a * b * n
    df <- ifelse(intx, N - a * b - ncov, N - a - b + 1 - ncov)
    stats::pt(q = stats::qt(alpha, df), df, lambda)
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
  ab <- c(a, b)
  METHOD <- paste0("Balanced two-way analysis of ", ifelse(ncov < 1, "", "co"),
                   "variance\n     simple effects power calculation",
                   ifelse(intx, " with interaction", ""))
  out <- list(`a, b` = ab, mmatrix = pss.matrix.format(mmatrix),
              cmatrix = pss.matrix.format(cmatrix), n = n,
              sd = sd, ncov = ncov, rho = rho, alpha = alpha, power = power,
              method = METHOD)

  # Print output as a power.htest object
  if (ncov < 1) out <- out[!names(out) %in% c("ncov", "rho")]
  structure(out, class = "power.htest")

}

