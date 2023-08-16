#' # Example 6.1
#' pss.prop.1samp(n = NULL, p0 = 0.2, pA = 0.3, power = 0.8, sided = "one")
#' # Example 6.2
#' pss.prop.1samp(n = NULL, p0 = 0.4, pA = 0.5, power = 0.8, sided = "one")
#'

pss.prop.1samp <- function (n = NULL, p0 = NULL, pA = NULL, alpha = 0.05,
                            power = NULL, sided = c("two", "one")) {
  if (sum(sapply(list(n, power, alpha), is.null)) != 1)
    stop("exactly one of 'n', 'alpha', and 'power' must be NULL")

  sided <- match.arg(sided)
  side <- switch(sided, one = 1, two = 2)

  p.body <- quote({
    d <- abs(pA - p0)
    (stats::pnorm(sqrt(n) * d / sqrt(pA * (1 - pA)) -
                  stats::qnorm(alpha / side, lower = FALSE)))
  })

  if (is.null(power))
    power <- eval(p.body)
  else if (is.null(n))
    n <- uniroot(function(n) eval(p.body) - power, c(2 + 1e-10, 1e+09))$root
  else if (is.null(alpha))
    alpha <- uniroot(function(alpha) eval(p.body) - power, c(1e-10, 1 - 1e-10))$root
  else stop("internal error", domain = NA)

  METHOD <- "Two-sample comparison of proportions power calculation"
  structure(list(n = n, p0 = p0, pA = pA, alpha = alpha,
                 power = power, sided = sided,
                 method = METHOD), class = "power.htest")
}

