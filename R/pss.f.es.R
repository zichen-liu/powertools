#' Cohen's f effect size calculation for one- or two- way analysis of variance
#'
#' @param means A vector or matrix of group means.
#' @param sd The estimated standard deviation within each group.
#'
#' @return A list of the arguments and various f effect sizes.
#' @export
#'
#' @examples
#' pss.f.es(means = c(5, 10, 12), sd = 10)
#' mmatrix <- matrix(c(9.3, 8.9, 8.5, 8.7, 8.3, 7.9), nrow = 2, byrow = TRUE)
#' pss.f.es(means = mmatrix, sd = 2)

pss.f.es <- function (means = NULL, sd = NULL) {

  # Check if the arguments are specified correctly
  if(is.null(sd))
    stop("sd must be specified")
  if (!is.vector(means) & !is.matrix(means))
    stop("a vector or matrix of means must be specified (see examples)")
  if (is.vector(means)) {
    a <- length(means)
    b <- 1
    if (a < 2)
      stop("number of groups must be at least 2")
    mmatrix <- matrix(means)
  }
  else if (is.matrix(means)) {
    a <- nrow(means)
    b <- ncol(means)
    if (a < 2 | b < 2)
      stop("number of groups per intervention must be at least 2")
    mmatrix <- means
  }

  # Get grand mean, marginal means, and interaction effects
  mu <- mean(mmatrix)
  temp1 <- mmatrix - mu
  mmA <- rowMeans(temp1)
  mmB <- colMeans(temp1)
  temp2 <- sweep(x = temp1, MARGIN = 2, STATS = mmB, FUN = "-")
  ints <- sweep(x = temp2, MARGIN = 1, STATS = mmA, FUN = "-")

  # Get sds and effect sizes
  sdA <- sqrt(sum(mmA^2) / a)
  sdB <- sqrt(sum(mmB^2) / b)
  sdAB <- sqrt(sum(ints^2) / (a * b))
  fA <- sdA / sd
  fB <- sdB / sd
  fAB <- sdAB / sd
  if (fAB < 0.000001) fAB <- 0

  # Print output as a power.htest object
  if (is.vector(means)) {
    METHOD <- "Cohen's f effect size calculation for\n     one-way analysis of variance"
    structure(list(a = a, means = means, sd = sd,
                   mmA = mmA, fA = fA,
                   method = METHOD), class = "power.htest")

  } else if (is.matrix(means)) {
    METHOD <- "Cohen's f effect size calculation for\n     two-way analysis of variance"
    mrows <- c()
    for (i in 1:a) mrows <- c(mrows, paste(means[i,], collapse = ', '))
    means <- paste(mrows, collapse = "\n                  ")
    structure(list(a = a, b = b, means = means, sd = sd,
                   mmA = mmA, mmB = mmB, ints = ints,
                   fA = fA, fB = fB, fAB = fAB,
                   method = METHOD), class = "power.htest")
  }

}
