pss.anova.2way <- function (n = NULL, means = NULL, sigma = NULL,
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
  if(is.null(sigma))
    stop("sigma must be specified")

  powerA <- power
  powerB <- power

  # Get degrees of freedom
  N <- n * a * b
  df <- N - a - b + 1

  # Get grand mean and marginal means
  mu <- mean(means)
  mmA <- rowMeans(means - mu)
  mmB <- colMeans(means - mu)

  # Get sigmas and f's
  sigmaA <- sqrt(sum(mmA^2) / a)
  sigmaB <- sqrt(sum(mmB^2) / b)
  fA <- sigmaA / sigma
  fB <- sigmaB / sigma

  # Copied from Example 5.8
  p.body.A <- quote({
    LambdaA <- N * fA^2
    stats::pf(q = stats::qf(alpha, a - 1, df, lower.tail = FALSE),
              a - 1, df, LambdaA, lower.tail = FALSE)
  })
  p.body.B <- quote({
    LambdaB <- N * fB^2
    stats::pf(q = stats::qf(alpha, b - 1, df, lower.tail = FALSE),
              b - 1, df, LambdaB, lower.tail = FALSE)
  })

  # Use uniroot function to calculate missing argument
  if (is.null(power)) {
    powerA <- eval(p.body.A)
    powerB <- eval(p.body.B)
    power <- min(powerA, powerB)
  }
  else if (is.null(n)){
    nA <- uniroot(function(n) eval(p.body.A) - power, c(2, 1e+05))$root
    nB <- uniroot(function(n) eval(p.body.B) - power, c(2, 1e+05))$root
    n <- max(nA, nB)
  }
  else if (is.null(alpha))
    alpha <- uniroot(function(alpha) eval(p.body.A) - power, c(1e-10, 1 - 1e-10))$root
  else stop("internal error", domain = NA)
  NOTE <- "power is the minimum power among two factors"
  METHOD <- "Balanced two-way analysis of variance power calculation"
  structure(list(a = a, b = b, n = n, means = means,
                 sigma = sigma, alpha = alpha,
                 powerA = powerA, powerB = powerB, power = power,
                 note = NOTE, method = METHOD), class = "power.htest")
}
