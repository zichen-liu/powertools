#' Power calculations for two-way balanced analysis of variance test for main effects
#'
#' @param n The sample size per group.
#' @param mmatrix A matrix of group means (see example).
#' @param sd The estimated standard deviation within each group; defaults to 1.
#' @param alpha The significance level or type 1 error rate; defaults to 0.05.
#' @param power The specified level of power.
#'
#' @return A list of the arguments (including the computed one).
#' @export
#'
#' @examples
#' # Example 5.8
#' mmatrix <- matrix(c(9.3, 8.9, 8.5, 8.7, 8.3, 7.9), nrow = 2, byrow = TRUE)
#' pss.anova.2way.main(n = 30, mmatrix = mmatrix, sd = 2, alpha = 0.05)

pss.anova.2way.main <- function (n = NULL, mmatrix = NULL, sd = 1,
                            alpha = 0.05, power = NULL) {

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
  nA <- n; nB <- n
  powerA <- power; powerB <- power

  # Get grand mean and marginal mmatrix
  mu <- mean(mmatrix)
  mmA <- rowMeans(mmatrix - mu)
  mmB <- colMeans(mmatrix - mu)

  # Get sds and f's
  sdA <- sqrt(sum(mmA^2) / a)
  sdB <- sqrt(sum(mmB^2) / b)
  fA <- sdA / sd
  fB <- sdB / sd

  # Calculate df and ncp for both factors
  p.body.A <- quote({
    N <- n * a * b
    df2 <- N - a - b + 1
    LambdaA <- N * fA^2
    stats::pf(q = stats::qf(alpha, a - 1, df2, lower.tail = FALSE),
              a - 1, df2, LambdaA, lower.tail = FALSE)
  })
  p.body.B <- quote({
    N <- n * a * b
    df2 <- N - a - b + 1
    LambdaB <- N * fB^2
    stats::pf(q = stats::qf(alpha, b - 1, df2, lower.tail = FALSE),
              b - 1, df2, LambdaB, lower.tail = FALSE)
  })

  # Use uniroot function to calculate missing argument
  if (is.null(power)) {
    powerA <- eval(p.body.A)
    powerB <- eval(p.body.B)
  }
  else if (is.null(n)){
    nA <- uniroot(function(n) eval(p.body.A) - power, c(2, 1e+05))$root
    nB <- uniroot(function(n) eval(p.body.B) - power, c(2, 1e+05))$root
  }
  else if (is.null(alpha))
    alpha <- uniroot(function(alpha) eval(p.body.A) - power, c(1e-10, 1 - 1e-10))$root
  else stop("internal error", domain = NA)

  # Generate output text
  NOTE <- "power is the minimum power among two factors;\n      n is the maximum required n among two factors"
  METHOD <- "Balanced two-way analysis of variance omnibus f test\n     power calculation for main effects only"
  mrows <- c()
  for (i in 1:a) mrows <- c(mrows, paste(mmatrix[i,], collapse = ', '))

  # Print output as a power.htest object
  structure(list(`a, b` = c(a, b), `nA, nB` = c(nA, nB),
                 mmatrix = paste(mrows, collapse = " | "),
                 sd = sd, `fA, fB` = c(fA, fB), alpha = alpha,
                 `powerA, powerB` = c(powerA, powerB),
                 note = NOTE, method = METHOD), class = "power.htest")
}

