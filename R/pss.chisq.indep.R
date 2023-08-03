#' # Example 8.2
#' pss.chisq.indep(p = matrix(c(0.050, 0.350, 0.100, 0.075, 0.250, 0.175), nrow = 2, byrow = TRUE), n = 230)
#' # Example 8.3
#' pss.chisq.indep(p = matrix(c(0.3, 0.2, 0.4, 0.1), nrow = 2, byrow = TRUE), n = 200)

pss.chisq.indep <- function (p = NULL, n = NULL, alpha = 0.05, power = NULL) {

  # Check if the arguments are specified correctly
  if (sum(sapply(list(n, alpha, power), is.null)) != 1)
    stop("exactly one of n, alpha, and power must be NULL")
  if (is.null(p))
    stop("please specify the proportion matrix p")
  if (!is.null(n) && any(n < 1))
    stop("number of observations must be at least 1")

  # Calculate effect size and df
  pi <- apply(p, 1, sum)
  pj <- apply(p, 2, sum)
  p0 <- pi %*% t(pj)
  es <- sqrt(sum((p - p0)^2/p0))
  df <- (nrow(p) - 1) * (ncol(p) - 1)

  # Calculate test statistic
  p.body <- quote({
    stats::pchisq(stats::qchisq(alpha, df, lower = FALSE),
                  df, n * es^2, lower = FALSE)
  })

  # Use uniroot function to calculate missing argument
  if (is.null(power))
    power <- eval(p.body)
  else if (is.null(n))
    n <- uniroot(function(n) eval(p.body) - power, c(1 + 1e-10, 1e+09))$root
  else if (is.null(alpha))
    alpha <- uniroot(function(sig.level) eval(p.body) - power, c(1e-10, 1 - 1e-10))$root
  else stop("internal error")

  # Generate output text
  METHOD <- "Chi squared power calculation"
  NOTE <- "n is the number of observations"
  structure(list(p = p, effect.size = es, n = n, alpha = alpha,
                 power = power, method = METHOD, note = NOTE), class = "power.htest")
}

