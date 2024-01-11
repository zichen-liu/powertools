#' Power calculations for two-way balanced analysis of variance contrast test
#'
#' @param n The sample size per group.
#' @param mmatrix A matrix of group means (see example).
#' @param cvec A vector of contrast coefficients c(c1, c2, ...).
#' @param factor Either "a" or "b" depending on which factor the contrast test is being made on.
#' @param sd The estimated standard deviation within each group; defaults to 1.
#' @param Rsq The estimated R^2 for regressing the outcome on the covariates; defaults to 0.
#' @param ncov The number of covariates adjusted for in the model; defaults to 0.
#' @param alpha The significance level or type 1 error rate; defaults to 0.05.
#' @param power The specified level of power.
#'
#' @return A list of the arguments (including the computed one).
#' @export
#'
#' @examples
#' mmatrix <- matrix(c(9.3, 8.9, 8.5, 8.7, 8.3, 7.3), nrow = 2, byrow = TRUE)
#' pss.anova2way.c.bal(n = 30, mmatrix = mmatrix, cvec = c(1, 0, -1), factor = "b",
#' sd = 2, alpha = 0.05)

pss.anova2way.c.bal <- function (n = NULL, mmatrix = NULL, cvec = NULL,
                                factor = c("a", "b"), sd = 1, Rsq = 0, ncov = 0,
                                alpha = 0.05, power = NULL) {

  # Check if the arguments are specified correctly
  a <- nrow(mmatrix)
  b <- ncol(mmatrix)
  factor <- match.arg(factor)
  if (sum(vapply(list(n, alpha, power), is.null, NA)) != 1)
    stop("exactly one of 'n', 'alpha', and 'power' must be NULL")
  if (a < 2 | b < 2)
    stop("number of groups per intervention must be at least 2")
  if (!is.null(n) && n < 2)
    stop("number of observations in each group must be at least 2")
  if (switch(factor, "a" = a, "b" = b) != length(cvec))
    stop("number of contrast coefficients must be equal to the number of groups")
  if(is.null(sd))
    stop("sd must be specified")

  # Get grand mean and marginal means
  es <- pss.es.anova.f(means = mmatrix, sd = sd)
  mmA <- es$mmA
  mmB <- es$mmB

  # See if there is an interaction
  fAB <- es$fAB
  intx <- ifelse(fAB == 0, FALSE, TRUE)

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
                   "variance\n     contrast test power calculation")
  out <- list(`a, b` = ab, mmatrix = pss.matrix.format(mmatrix),
              factor = factor, cvec = cvec, n = n,
              sd = sd, ncov = ncov, Rsq = Rsq,
              alpha = alpha, power = power,
              method = METHOD)

  # Print output as a power.htest object
  if (ncov < 1) out <- out[!names(out) %in% c("ncov", "Rsq")]
  structure(out, class = "power.htest")

}

