#' Power calculations for two-way unbalanced analysis of variance test for main effects and interaction effect
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
#' # Example 5.10
#' nmatrix <- matrix(c(30, 30, 30, 30, 30, 30), nrow = 2, byrow = TRUE)
#' mmatrix <- matrix(c(9.3, 8.9, 8.5, 8.7, 8.3, 7.3), nrow = 2, byrow = TRUE)
#' pss.anova.2w.intx(nmatrix = nmatrix, mmatrix = mmatrix, sd = 2, alpha = 0.05)

pss.anova.2w.intx <- function (nmatrix = NULL, mmatrix = NULL, sd = 1,
                               alpha = 0.05, power = NULL) {

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
  means <- mmatrix - mu
  mmA <- rowMeans(means)
  mmB <- colMeans(means)
  temp <- sweep(x = means, MARGIN = 2, STATS = mmB, FUN = "-")
  ints <- sweep(x = temp, MARGIN = 1, STATS = mmA, FUN = "-")

  # Get Lambdas
  LambdaA <- 0
  for (i in 1:a) {
    for (j in 1:b) {
      temp <- mmA[i] / (sd / sqrt(nmatrix[i, j]))
      LambdaA <- LambdaA + temp^2
    }
  }
  LambdaB <- 0
  for (j in 1:b) {
    for (i in 1:a) {
      temp <- mmB[j] / (sd / sqrt(nmatrix[i, j]))
      LambdaB <- LambdaB + temp^2
    }
  }
  LambdaAB <- 0
  for (i in 1:a) {
    for (j in 1:b) {
      temp <- ints[i, j] / (sd / sqrt(nmatrix[i, j]))
      LambdaAB <- LambdaAB + temp^2
    }
  }

  # Get f effect sizes
  sdA <- sqrt(sum(mmA^2) / a)
  sdB <- sqrt(sum(mmB^2) / b)
  sdAB <- sqrt(sum(ints^2) / (a * b))
  fA <- sdA / sd
  fB <- sdB / sd
  fAB <- sdAB / sd

  # Calculate power
  N <- sum(nmatrix)
  df2 <- N - a * b
  powerA <- stats::pf(stats::qf(alpha, a - 1, df2, lower.tail = FALSE),
                      a - 1, df2, LambdaA, lower.tail = FALSE)
  powerB <- stats::pf(stats::qf(alpha, b - 1, df2, lower.tail = FALSE),
                      b - 1, df2, LambdaB, lower.tail = FALSE)
  powerAB <- stats::pf(stats::qf(alpha, (a - 1) * (b - 1), df2, lower.tail = FALSE),
                       (a - 1) * (b - 1), df2, LambdaAB, lower.tail = FALSE)

  # Generate output text
  power <- c(powerA, powerB)
  ab <- c(a, b)
  f <- c(fA, fB)
  METHOD <- "Unalanced two-way analysis of variance power calculation\n     for main effects and interaction effect"
  nrows <- c(); mrows <- c()
  for (i in 1:a) nrows <- c(nrows, paste(nmatrix[i,], collapse = ', '))
  nmatrix <- paste(nrows, collapse = " | ")
  for (i in 1:a) mrows <- c(mrows, paste(mmatrix[i,], collapse = ', '))
  mmatrix <- paste(mrows, collapse = " | ")

  # Print output as a power.htest object
  structure(list(`a, b` = ab, nmatrix = nmatrix, mmatrix = mmatrix,
                 sd = sd, alpha = alpha, f = f, power = power,
                 f.int = fAB, power.int = powerAB,
                 method = METHOD), class = "power.htest")
}

