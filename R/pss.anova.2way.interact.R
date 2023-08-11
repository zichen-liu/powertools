#' Power calculations for two-way balanced analysis of variance test for main effects and interaction effect
#'
#' @param n The sample size per group.
#' @param means A matrix of group means (see example).
#' @param sd The estimated standard deviation within each group; defaults to 1.
#' @param alpha The significance level or type 1 error rate; defaults to 0.05.
#' @param power The specified level of power.
#'
#' @return A list of the arguments (including the computed one).
#' @export
#'
#' @examples
#' # Example 5.10
#' means <- matrix(c(9.3, 8.9, 8.5, 8.7, 8.3, 7.3), nrow = 2, byrow = TRUE)
#' pss.anova.2way.interact(n = 30, means = means, sd = 2, alpha = 0.05)

pss.anova.2way.interact <- function (n = NULL, means = NULL, sd = 1,
                            alpha = 0.05, power = NULL) {

  # Check if the arguments are specified correctly
  a <- nrow(means)
  b <- ncol(means)
  if (sum(vapply(list(n, alpha, power), is.null, NA)) != 1)
    stop("exactly one of 'n', 'alpha', and 'power' must be NULL")
  if (a < 2 | b < 2)
    stop("number of groups per intervention must be at least 2")
  if (!is.null(n) && n < 2)
    stop("number of observations in each group must be at least 2")
  if(is.null(sd))
    stop("sd must be specified")

  # Set default values if given
  nA <- n; nB <- n; nAB <- n
  powerA <- power; powerB <- power; powerAB <- power

  # Get grand mean, marginal means, and interaction effects
  mu <- mean(means)
  means2 <- means - mu
  mmA <- rowMeans(means2)
  mmB <- colMeans(means2)
  ints <- sweep(x = means2, MARGIN = 2, STATS = mmB, FUN = "-")
  ints <- sweep(x = ints, MARGIN = 1, STATS = mmA, FUN = "-")

  # Get sds and f's
  sdA <- sqrt(sum(mmA^2) / a)
  sdB <- sqrt(sum(mmB^2) / b)
  sdAB <- sqrt(sum(ints^2) / (a * b))
  fA <- sdA / sd
  fB <- sdB / sd
  fAB <- sdAB / sd

  # Calculate df and ncp for both factors and the interaction
  p.body.A <- quote({
    N <- n * a * b
    df2 <- N - a  * b
    LambdaA <- N * fA^2
    stats::pf(q = stats::qf(alpha, a - 1, df2, lower.tail = FALSE),
              a - 1, df2, LambdaA, lower.tail = FALSE)
  })
  p.body.B <- quote({
    N <- n * a * b
    df2 <- N - a * b
    LambdaB <- N * fB^2
    stats::pf(q = stats::qf(alpha, b - 1, df2, lower.tail = FALSE),
              b - 1, df2, LambdaB, lower.tail = FALSE)
  })
  p.body.AB <- quote({
    N <- n * a * b
    df2 <- N - a * b
    LambdaAB <- N * fAB^2
    stats::pf(q = stats::qf(alpha, (a - 1) * (b - 1), df2, lower.tail = FALSE),
              (a - 1) * (b - 1), df2, LambdaAB, lower.tail = FALSE)
  })

  # Use uniroot function to calculate missing argument
  if (is.null(power)) {
    powerA <- eval(p.body.A)
    powerB <- eval(p.body.B)
    powerAB <- eval(p.body.AB)
  }
  else if (is.null(n)){
    nA <- uniroot(function(n) eval(p.body.A) - power, c(2, 1e+05))$root
    nB <- uniroot(function(n) eval(p.body.B) - power, c(2, 1e+05))$root
    nAB <- uniroot(function(n) eval(p.body.AB) - power, c(2, 1e+05))$root
  }
  else if (is.null(alpha))
    alpha <- uniroot(function(alpha) eval(p.body.A) - power, c(1e-10, 1 - 1e-10))$root
  else stop("internal error", domain = NA)

  # Generate output text
  METHOD <- "Balanced two-way analysis of variance power calculation\n     for main effects and interaction effect"
  mrows <- c()
  for (i in 1:a) mrows <- c(mrows, paste(means[i,], collapse = ', '))

  # Print output as a power.htest object
  structure(list(`a, b` = c(a, b), `nA, nB, nAB` = c(nA, nB, nAB),
                 means = paste(mrows, collapse = " | "),
                 sd = sd, `fA, fB, fAB` = c(fA, fB, fAB), alpha = alpha,
                 `powerA, powerB, powerAB` = c(powerA, powerB, powerAB),
                 method = METHOD), class = "power.htest")
}

