#' Power calculations for two-way unbalanced analysis of variance omnibus F test
#'
#' @param nmatrix A matrix of group sample sizes (see example).
#' @param mmatrix A matrix of group means (see example).
#' @param sd The estimated standard deviation within each group.
#' @param rho The estimated correlation between covariates and the outcome; defaults to 0.
#' @param ncov The number of covariates adjusted for in the model; defaults to 0.
#' @param alpha The significance level or type 1 error rate; defaults to 0.05.
#'
#' @return A list of the arguments (including the computed power).
#' @export
#'
#' @examples
#' # Example 5.8
#' nmatrix <- matrix(c(30, 30, 30, 30, 30, 30), nrow = 2, byrow = TRUE)
#' mmatrix <- matrix(c(9.3, 8.9, 8.5, 8.7, 8.3, 7.9), nrow = 2, byrow = TRUE)
#' pss.anova.unbal.2w(nmatrix = nmatrix, mmatrix = mmatrix, sd = 2, alpha = 0.05)
#' # Example 5.10
#' nmatrix <- matrix(c(30, 30, 30, 30, 30, 30), nrow = 2, byrow = TRUE)
#' mmatrix <- matrix(c(9.3, 8.9, 8.5, 8.7, 8.3, 7.3), nrow = 2, byrow = TRUE)
#' pss.anova.unbal.2w(nmatrix = nmatrix, mmatrix = mmatrix, sd = 2, alpha = 0.05)
#' # Example 5.14
#' nmatrix <- matrix(c(30, 30, 30, 30, 30, 30), nrow = 2, byrow = TRUE)
#' mmatrix <- matrix(c(9.3, 8.9, 8.5, 8.7, 8.3, 7.9), nrow = 2, byrow = TRUE)
#' pss.anova.unbal.2w(nmatrix = nmatrix, mmatrix = mmatrix, sd = 2, rho = 0.4, ncov = 1, alpha = 0.05)

pss.anova.unbal.2w <- function (nmatrix = NULL, mmatrix = NULL, sd = NULL,
                                rho = 0, ncov = 0, alpha = 0.05) {

  # Check if the arguments are specified correctly
  a <- nrow(mmatrix)
  b <- ncol(mmatrix)
  if (a < 2 | b < 2)
    stop("number of groups per factor must be at least 2")
  if (any(nmatrix < 2))
    stop("number of observations in each group must be at least 2")
  if(is.null(nmatrix) | is.null(mmatrix))
    stop("sample size matrix and means matrix must both be specified")
  if(is.null(sd))
    stop("sd must be specified")
  if(nrow(nmatrix) != a | ncol(nmatrix) != b)
    stop("number of sample sizes must equal to the number of groups")

  # Get marginal means
  es <- pss.anova.f.es(means = mmatrix, sd = sd)
  mmA <- es$mmA
  mmB <- es$mmB
  ints <- es$ints

  # Get f effect sizes
  fA <- es$fA
  fB <- es$fB
  fAB <- es$fAB
  intx <- ifelse(fAB == 0, FALSE, TRUE)

  # Get Lambdas
  LambdaA <- 0
  for (i in 1:a) {
    for (j in 1:b) {
      temp <- mmA[i] / (sd / sqrt(nmatrix[i, j]))
      LambdaA <- LambdaA + temp^2
    }
  }
  LambdaA <- LambdaA / (1 - rho^2)
  LambdaB <- 0
  for (j in 1:b) {
    for (i in 1:a) {
      temp <- mmB[j] / (sd / sqrt(nmatrix[i, j]))
      LambdaB <- LambdaB + temp^2
    }
  }
  LambdaB <- LambdaB / (1 - rho^2)
  LambdaAB <- 0
  for (i in 1:a) {
    for (j in 1:b) {
      temp <- ints[i, j] / (sd / sqrt(nmatrix[i, j]))
      LambdaAB <- LambdaAB + temp^2
    }
  }
  LambdaAB <- LambdaAB / (1 - rho^2)

  # Calculate power
  N <- sum(nmatrix)
  df1A <- a - 1; df1B <- b - 1; df1AB <- (a - 1) * (b - 1)
  df2 <- ifelse(intx, N - a * b - ncov, N - a - b + 1 - ncov)
  powerA <- stats::pf(stats::qf(alpha, df1A, df2, lower.tail = FALSE),
                      df1A, df2, LambdaA, lower.tail = FALSE)
  powerB <- stats::pf(stats::qf(alpha, df1B, df2, lower.tail = FALSE),
                      df1B, df2, LambdaB, lower.tail = FALSE)
  if (intx) {
    powerAB <- stats::pf(stats::qf(alpha, df1AB, df2, lower.tail = FALSE),
                         df1AB, df2, LambdaAB, lower.tail = FALSE)
  }
  else { powerAB <- 0}

  # Generate output text
  ab <- c(a, b)
  power <- c(powerA, powerB)
  f <- c(fA, fB)
  METHOD <- paste0("Unalanced two-way analysis of ", ifelse(ncov < 1, "", "co"),
                   "variance\n     omnibus f test power calculation",
                   ifelse(intx, " with interaction", ""))
  out <- list(`a, b` = ab, mmatrix = pss.matrix.format(mmatrix),
              nmatrix = pss.matrix.format(nmatrix),
              sd = sd, ncov = ncov, rho = rho,
              alpha = alpha, f = f, power = power,
              f.int = fAB, power.int = powerAB,
              method = METHOD)

  # Print output as a power.htest object
  if (ncov < 1) out <- out[!names(out) %in% c("ncov", "rho")]
  if (!intx) out <- out[!names(out) %in% c("n.int", "f.int", "power.int")]
  structure(out, class = "power.htest")

}

