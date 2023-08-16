#' Power calculations for two-way unbalanced analysis of variance omnibus F test
#'
#' @param nmatrix A matrix of group sample sizes (see example).
#' @param mmatrix A matrix of group means (see example).
#' @param sd The estimated standard deviation within each group.
#' @param alpha The significance level or type 1 error rate; defaults to 0.05.
#'
#' @return A list of the arguments (including the computed power).
#' @export
#'
#' @examples
#' # Example 5.8
#' nmatrix <- matrix(c(30, 30, 30, 30, 30, 30), nrow = 2, byrow = TRUE)
#' mmatrix <- matrix(c(9.3, 8.9, 8.5, 8.7, 8.3, 7.9), nrow = 2, byrow = TRUE)
#' pss.unbal.anova.2way.main(nmatrix = nmatrix, mmatrix = mmatrix, sd = 2, alpha = 0.05)

pss.unbal.anova.2way.main <- function (nmatrix = NULL, mmatrix = NULL, sd = NULL, alpha = 0.05) {

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
  mu <- mean(mmatrix)
  mmA <- rowMeans(mmatrix - mu)
  mmB <- colMeans(mmatrix - mu)

  # Get LambdaA
  LambdaA <- 0
  for (i in 1:a) {
    for (j in 1:b) {
      temp <- mmA[i] / (sd / sqrt(nmatrix[i, j]))
      LambdaA <- LambdaA + temp^2
    }
  }

  # Get LambdaB
  LambdaB <- 0
  for (j in 1:b) {
    for (i in 1:a) {
      temp <- mmB[j] / (sd / sqrt(nmatrix[i, j]))
      LambdaB <- LambdaB + temp^2
    }
  }

  # Get f effect sizes
  sdA <- sqrt(sum(mmA^2) / a)
  sdB <- sqrt(sum(mmB^2) / b)
  fA <- sdA / sd
  fB <- sdB / sd

  # Calculate power
  N <- sum(nmatrix)
  df2 <- N - a - b + 1
  powerA <- stats::pf(stats::qf(alpha, a - 1, df2, lower.tail = FALSE),
                      a - 1, df2, LambdaA, lower.tail = FALSE)
  powerB <- stats::pf(stats::qf(alpha, b - 1, df2, lower.tail = FALSE),
                      b - 1, df2, LambdaB, lower.tail = FALSE)
  power <- min(powerA, powerB)

  # Generate output text
  NOTE <- "power is the minimum power among two factors"
  METHOD <- "Unbalanced two-way analysis of variance omnibus f test\n     power calculation for main effects only"
  nrows <- c(); mrows <- c()
  for (i in 1:a) nrows <- c(nrows, paste(nmatrix[i,], collapse = ', '))
  for (i in 1:a) mrows <- c(mrows, paste(mmatrix[i,], collapse = ', '))

  # Print output as a power.htest object
  structure(list(`a, b` = c(a, b),
                 nmatrix = paste(nrows, collapse = " | "),
                 mmatrix = paste(mrows, collapse = " | "),
                 sd = sd, `fA, fB` = c(fA, fB), alpha = alpha,
                 `powerA, powerB` = c(powerA, powerB), power = power,
                 note = NOTE, method = METHOD), class = "power.htest")
}

