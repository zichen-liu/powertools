#' # Example 6.6
#' pss.prop.2samp(n = NULL, p1 = 0.6, p2 = 0.8, alpha = 0.025, power = 0.9, sided = "one")

pss.prop.2samp <- function (n = NULL, p1 = NULL, p2 = NULL, alpha = 0.05,
          power = NULL, n.ratio = 1, sided = c("two", "one")) {
  if (sum(sapply(list(n, power, alpha, n.ratio), is.null)) != 1)
    stop("exactly one of 'n', 'n.ratio', 'alpha', and 'power' must be NULL")
  if (!is.null(n.ratio) && n.ratio <= 0)
    stop("n.ratio between group sizes must be positive")
  sided <- match.arg(sided)
  side <- switch(sided, one = 1, two = 2)
  p.body <- quote({
    d <- abs(p1 - p2)
    q1 <- 1 - p1
    q2 <- 1 - p2
    ((stats::qnorm(alpha / side) + stats::qnorm(1 - power))^2 *
    (n.ratio * p1 * q1 + p2 * q2) / (n.ratio * d^2))
  })

  tol = .Machine$double.eps^0.25
  if (is.null(n))
    n <- eval(p.body)
  else if (is.null(power))
    power <- uniroot(function(power) eval(p.body) - n, c(1e-05, 0.99999), tol = tol, extendInt = "upX")$root
  else if (is.null(n.ratio))
    n.ratio <- uniroot(function(n.ratio) eval(p.body) - n, c(2/n, 1e+07))$root
  else if (is.null(alpha))
    alpha <- uniroot(function(alpha) eval(p.body) - n, c(1e-10, 1 - 1e-10), tol = tol, extendInt = "upX")$root
  else stop("internal error", domain = NA)

  NOTE <- "n is the number in each group"
  METHOD <- "Two-sample comparison of proportions power calculation"
  structure(list(n = c(n, n * n.ratio), p1 = p1, p2 = p2, alpha = alpha,
                 power = power, sided = sided, note = NOTE,
                 method = METHOD), class = "power.htest")
}

