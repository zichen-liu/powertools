#' Power calculations for two-way unbalanced analysis of variance contrast test
#'
#' @param nmatrix A matrix of group sample sizes (see example).
#' @param mmatrix A matrix of group means (see example).
#' @param cvec A vector of contrast coefficients c(c1, c2, ...).
#' @param factor Either "a" or "b" depending on which factor the contrast test is being made on.
#' @param sd The estimated standard deviation within each group.
#' @param indx Whether there is an interaction between the two factors.
#' @param alpha The significance level or type 1 error rate; defaults to 0.05.
#' @param power The specified level of power.
#'
#' @return A list of the arguments (including the computed one).
#' @export
#'
#' @examples
#' # Example 5.11
#' nmatrix <- matrix(c(30, 30, 30, 30, 30, 30), nrow = 2, byrow = TRUE)
#' mmatrix <- matrix(c(9.3, 8.9, 8.5, 8.7, 8.3, 7.3), nrow = 2, byrow = TRUE)
#' pss.anova.unbal.2w.c(nmatrix = nmatrix, mmatrix = mmatrix, cvec = c(1, 0, -1), factor = "b", sd = 2, intx = TRUE, alpha = 0.05)

pss.anova.unbal.2w.c <- function (nmatrix = nmatrix, mmatrix = NULL, cvec = NULL,
                            factor = c("a", "b"), sd = NULL, intx = FALSE,
                            alpha = 0.05) {

  # Check if the arguments are specified correctly
  a <- nrow(mmatrix)
  b <- ncol(mmatrix)
  factor <- match.arg(factor)
  if (a < 2 | b < 2)
    stop("number of groups per intervention must be at least 2")
  if (!is.null(n) && n < 2)
    stop("number of observations in each group must be at least 2")
  if (switch(factor, "a" = a, "b" = b) != length(cvec))
    stop("number of contrast coefficients must be equal to the number of groups")
  if(is.null(sd))
    stop("sd must be specified")

  # Get grand mean and marginal means
  mu <- mean(mmatrix)
  mmA <- rowMeans(mmatrix - mu)
  mmB <- colMeans(mmatrix - mu)

  # Get lambda (Lambda = lambda^2)
  num <- switch(factor, "a" = mmA, "b" = mmB) %*% cvec
  temp <- switch(factor,
                 "a" = sapply(X = 1:a, FUN = function(i)
                       cvec[i]^2 / sum(nmatrix[i, ])),
                 "b" = sapply(X = 1:b, FUN = function(j)
                       cvec[j]^2 / sum(nmatrix[, j])))
  den <- sd * sqrt(sum(temp))
  lambda <- num / den
  print(lambda^2)

  # Calculate power
  N <- sum(nmatrix)
  df2 <- ifelse(intx, N - a * b, N - a - b + 1)
  power <- stats::pf(q = stats::qf(alpha, 1, df2, lower.tail = FALSE),
                     1, df2, lambda^2, lower.tail = FALSE)

  # Generate output text
  ab <- c(a, b)
  METHOD <- "Unbalanced two-way analysis of variance\n     contrast test power calculation"
  mrows <- c()
  for (i in 1:a) mrows <- c(mrows, paste(mmatrix[i,], collapse = ', '))
  mmatrix <- paste(mrows, collapse = " | ")

  # Print output as a power.htest object
  structure(list(`a, b` = ab, mmatrix = mmatrix,
                 factor = factor, cvec = cvec, n = n,
                 sd = sd, alpha = alpha, power = power,
                 method = METHOD), class = "power.htest")
}

