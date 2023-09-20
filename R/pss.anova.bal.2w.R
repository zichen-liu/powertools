#' Power calculations for two-way balanced analysis of variance omnibus F test
#'
#' @param n The sample size per group.
#' @param mmatrix A matrix of group means (see example).
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
#' # Example 5.8
#' mmatrix <- matrix(c(9.3, 8.9, 8.5, 8.7, 8.3, 7.9), nrow = 2, byrow = TRUE)
#' pss.anova.bal.2w(n = 30, mmatrix = mmatrix, sd = 2, alpha = 0.05)
#' # Example 5.10
#' mmatrix <- matrix(c(9.3, 8.9, 8.5, 8.7, 8.3, 7.3), nrow = 2, byrow = TRUE)
#' pss.anova.bal.2w(n = 30, mmatrix = mmatrix, sd = 2, alpha = 0.05)
#' # Example 5.14
#' mmatrix <- matrix(c(9.3, 8.9, 8.5, 8.7, 8.3, 7.9), nrow = 2, byrow = TRUE)
#' pss.anova.bal.2w(n = 30, mmatrix = mmatrix, sd = 2, Rsq = 0.4, ncov = 1, alpha = 0.05)

pss.anova.bal.2w <- function (n = NULL, mmatrix = NULL, sd = 1,
                              Rsq = 0, ncov = 0, alpha = 0.05, power = NULL) {

  # Check if the arguments are specified correctly
  a <- nrow(mmatrix)
  b <- ncol(mmatrix)
  if (sum(vapply(list(n, alpha, power), is.null, NA)) != 1)
    stop("exactly one of 'n', 'alpha', and 'power' must be NULL")
  if (a < 2 | b < 2)
    stop("number of groups per factor must be at least 2")
  if (!is.null(n) && n < 2)
    stop("number of observations in each group must be at least 2")
  if(is.null(sd))
    stop("sd must be specified")

  # Set default values if given
  nA <- n; nB <- n; nAB <- n
  powerA <- power; powerB <- power; powerAB <- power

  # Get f effect sizes
  es <- pss.anova.f.es(means = mmatrix, sd = sd)
  fA <- es$fA
  fB <- es$fB
  fAB <- es$fAB
  intx <- ifelse(fAB == 0, FALSE, TRUE)

  # Calculate df's and ncp's
  p.body.A <- quote({
    N <- n * a * b
    LambdaA <- N * fA^2 / (1 - Rsq)
    df1 <- a - 1
    df2 <- ifelse(intx, N - a * b - ncov, N - a - b + 1 - ncov)
    stats::pf(q = stats::qf(alpha, df1, df2, lower.tail = FALSE),
              df1, df2, LambdaA, lower.tail = FALSE)
  })
  p.body.B <- quote({
    N <- n * a * b
    LambdaB <- N * fB^2 / (1 - Rsq)
    df1 <- b - 1
    df2 <- ifelse(intx, N - a * b - ncov, N - a - b + 1 - ncov)
    stats::pf(q = stats::qf(alpha, df1, df2, lower.tail = FALSE),
              df1, df2, LambdaB, lower.tail = FALSE)
  })
  p.body.AB <- quote({
    N <- n * a * b
    LambdaAB <- N * fAB^2 / (1 - Rsq)
    df1 <- (a - 1) * (b - 1)
    df2 <- N - a * b - ncov
    stats::pf(q = stats::qf(alpha, df1, df2, lower.tail = FALSE),
              df1, df2, LambdaAB, lower.tail = FALSE)
  })

  # Use uniroot function to calculate missing argument
  if (is.null(power)) {
    powerA <- round(eval(p.body.A), 4)
    powerB <- round(eval(p.body.B), 4)
    if (intx) powerAB <- round(eval(p.body.AB), 4)
  }
  else if (is.null(n)){
    nA <- round(uniroot(function(n) eval(p.body.A) - power, c(2, 1e+05))$root, 4)
    nB <- round(uniroot(function(n) eval(p.body.B) - power, c(2, 1e+05))$root, 4)
    if (intx)
      nAB <- round(uniroot(function(n) eval(p.body.AB) - power, c(2, 1e+05))$root, 4)
  }
  else if (is.null(alpha))
    alpha <- uniroot(function(alpha) eval(p.body.A) - power, c(1e-10, 1 - 1e-10))$root
  else stop("internal error", domain = NA)

  # Generate output text
  ab <- c(a, b)
  if (is.null(power))
    if (intx) power <- c(powerA, powerB, powerAB) else power <- c(powerA, powerB)
  else if (is.null(n))
    if (intx) n <- c(nA, nB, nAB) else n <- c(nA, nB)
  if (intx) f <- c(round(fA, 4), round(fB, 4), round(fAB, 4))
  else f <- c(round(fA, 4), round(fB, 4))
  METHOD <- paste0("Balanced two-way analysis of ", ifelse(ncov < 1, "", "co"),
                   "variance\n     omnibus F test power calculation",
                   ifelse(intx, " with interaction", ""))
  NOTE <- "The 3rd value for f and power or n is for the interaction"
  out <- list(`a, b` = ab, mmatrix = pss.matrix.format(mmatrix),
              n = n, sd = sd, ncov = ncov, Rsq = Rsq,
              alpha = alpha, f = f, power = power,
              method = METHOD, note = NOTE)

  # Print output as a power.htest object
  if (ncov < 1) out <- out[!names(out) %in% c("ncov", "Rsq")]
  if (!intx) out <- out[!names(out) == "note"]
  structure(out, class = "power.htest")

}

