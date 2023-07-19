pss.anova.1way <- function (n = NULL, means = NULL, sigma = 1,
                             alpha = 0.05, power = NULL) {

  # Check if the arguments are specified correctly
  groups <- length(means)
  if (sum(vapply(list(n, sigma, alpha, power), is.null, NA)) != 1)
    stop("exactly one of 'n', 'sigma', 'alpha', and 'power' must be NULL")
  if (!is.null(groups) && groups < 2)
    stop("number of groups must be at least 2")
  if (!is.null(n) && n < 2)
    stop("number of observations in each group must be at least 2")

  # Get between group variance; sigma is within group variance
  between <- var(means)

  # Copied from power.anova.test
  p.body <- quote({
    lambda <- (groups - 1) * n * (between / sigma)
    stats::pf(stats::qf(alpha, groups - 1, (n - 1) * groups, lower.tail = FALSE),
              groups - 1, (n - 1) * groups, lambda, lower.tail = FALSE)
  })

  # Use uniroot function to calculate missing argument
  if (is.null(power))
    power <- eval(p.body)
  else if (is.null(n))
    n <- uniroot(function(n) eval(p.body) - power, c(2, 1e+05))$root
  else if (is.null(sigma))
    sigma <- uniroot(function(sigma) eval(p.body) - power, between * c(1e-07, 1e+07))$root
  else if (is.null(between))
    between <- uniroot(function(between) eval(p.body) - power, sigma * c(1e-07, 1e+07))$root
  else if (is.null(alpha))
    alpha <- uniroot(function(alpha) eval(p.body) - power, c(1e-10, 1 - 1e-10))$root
  else stop("internal error", domain = NA)
  NOTE <- "n is the number in each group"
  METHOD <- "Balanced one-way analysis of variance power calculation"
  structure(list(groups = groups, n = n, between = between,
                 sigma = sigma, alpha = alpha, power = power,
                 note = NOTE, method = METHOD), class = "power.htest")
}
