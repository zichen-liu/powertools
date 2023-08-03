#' # Example 8.1
#' pss.chisq.gof(p0 = c(0.5, 0.3, 0.2), p1 = c(0.7, 0.2, 0.1), n = 50)

pss.chisq.gof <- function(p0 = NULL, p1 = NULL,
                          n = NULL, alpha = 0.05, power = NULL) {

  # Check if the arguments are specified correctly
  if (sum(sapply(list(n, alpha, power), is.null)) != 1)
    stop("exactly one of n, alpha, and power must be NULL")
  if (is.null(p0) | is.null(p1))
    stop("please specify the two proportion vectors p0 and p1")
  if (!is.null(n) && any(n < 1))
    stop("number of observations must be at least 1")
  if (length(p0) != length(p1))
    stop("the two proportion vectors must have the same lengths")

  # Calculate effect size and df
  es <- sqrt(sum((p1 - p0)^2 / p0))
  df <- length(p0) - 1

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
  structure(list(p1 = p1, p0 = p0, effect.size = es, n = n, alpha = alpha,
                 power = power, method = METHOD, note = NOTE), class = "power.htest")
}
