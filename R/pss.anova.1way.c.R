pss.anova.1way.c <- function (n = NULL, means = NULL, coeff = NULL, sigma = NULL,
                            alpha = 0.05, power = NULL) {

  # Check if the arguments are specified correctly
  groups <- length(means)
  if (sum(vapply(list(n, coeff, alpha, power), is.null, NA)) != 1)
    stop("exactly one of 'n', 'coeff', 'alpha', and 'power' must be NULL")
  if (groups < 2)
    stop("number of groups must be at least 2")
  if (!is.null(n) && n < 2)
    stop("number of observations in each group must be at least 2")
  if (groups != length(coeff))
    stop("number of contrast coefficients must be equal to the number of groups")
  if(is.null(sigma))
    stop("sigma must be specified")

  # Create nvec
  nvec <- rep(n, groups)

  # Copied from Example 5.7
  p.body <- quote({
    lambda <- coeff %*% means / (sigma * sqrt(sum(coeff^2 / nvec)))
    df <- sum(nvec) - groups
    stats::pf(q = stats::qf(alpha, 1, df, lower.tail = FALSE),
              1, df, lambda^2, lower.tail = FALSE)
  })

  # Use uniroot function to calculate missing argument
  if (is.null(power))
    power <- eval(p.body)
  else if (is.null(n))
    n <- uniroot(function(n) eval(p.body) - power, c(2, 1e+05))$root
  else if (is.null(alpha))
    alpha <- uniroot(function(alpha) eval(p.body) - power, c(1e-10, 1 - 1e-10))$root
  else stop("internal error", domain = NA)
  NOTE <- "n is the number in each group"
  METHOD <- "Balanced one-way analysis of variance with contrasts power calculation"
  structure(list(groups = groups, n = n, means = means, coeff = coeff,
                 sigma = sigma, alpha = alpha, power = power,
                 note = NOTE, method = METHOD), class = "power.htest")
}

# pss.anova.1way.c(n = 20, means = c(5, 10, 20), c = c(1, -1, 0), sigma = 10, alpha = 0.025)
