#' Power calculations for two-way balanced analysis of variance contrast test
#'
#' @param n The sample size per group.
#' @param mmatrix A matrix of group means (see example).
#' @param cvec A vector of contrast coefficients c(c1, c2, ...).
#' @param factor Either "a" for rows or "b" for columns depending on which factor the contrast test is being made on.
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
#' anova2way.c.bal(n = 30, mmatrix = mmatrix, cvec = c(1, 0, -1), factor = "b",
#' sd = 2, alpha = 0.05)

anova2way.c.bal <- function (n = NULL, mmatrix = NULL, cvec = NULL,
                                factor = c("a", "b"), sd = 1, Rsq = 0, ncov = 0,
                                alpha = 0.05, power = NULL, v = FALSE) {

  # Check if the arguments are specified correctly
  check.many(list(n, alpha, power), "oneof")
  check(n, "pos"); check(n, "min", min = 2)
  check(mmatrix, "req"); check(mmatrix, "mat")
  check(cvec, "req"); check(cvec, "vec")
  check(factor, "req"); check(factor, "vals", valslist = c("a", "b"))
  check(sd, "req"); check(sd, "pos")
  check(Rsq, "req"); check(Rsq, "uniti")
  check(ncov, "req"); check(ncov, "int")
  check(alpha, "unit")
  check(power, "unit")
  check(v, "req"); check(v, "bool")

  a <- nrow(mmatrix)
  b <- ncol(mmatrix)
  factor <- match.arg(factor)
  if (switch(factor, "a" = a, "b" = b) != length(cvec))
    stop("number of contrast coefficients must be equal to the number of groups")

  if (Rsq > 0 & ncov == 0)
    stop("please specify ncov or set Rsq to 0")

  # See if there is an interaction
  fAB <- es.anova.f(means = mmatrix, sd = sd)$fAB
  intx <- ifelse(fAB == 0, FALSE, TRUE)

  # Get grand mean and marginal means
  mu <- mean(mmatrix)
  temp1 <- mmatrix - mu
  mmA <- rowMeans(temp1)
  mmB <- colMeans(temp1)

  # Calculate df's and ncp's
  p.body <- quote({
    temp <- switch(factor, "a" = mmA, "b" = mmB) %*% cvec
    nj <- n * switch(factor, "a" = b, "b" = a)
    Lambda <- temp^2 / (sd^2 * (1 / nj + 1 / nj)) / (1 - Rsq)
    N <- a * b * n
    df2 <- ifelse(intx, N - a * b - ncov, N - a - b + 1 - ncov)
    stats::pf(q = stats::qf(alpha, 1, df2, lower.tail = FALSE),
              1, df2, Lambda, lower.tail = FALSE)
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
                   "variance\n     contrast test power calculation")
  out <- list(n = n, mmatrix = matrix.format(mmatrix),
              factor = factor, cvec = cvec,
              sd = sd, ncov = ncov, Rsq = Rsq,
              alpha = alpha, power = power,
              method = METHOD)

  # Print output as a power.htest object
  if (ncov < 1) out <- out[!names(out) %in% c("ncov", "Rsq")]
  structure(out, class = "power.htest")

}

