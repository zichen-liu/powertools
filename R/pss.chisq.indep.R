#' Power calculations for chi-square test of independence
#'
#' @param pmatrix The two-way probability table under the alternative hypothesis.
#' @param N The number of total observations.
#' @param alpha The significance level or type 1 error rate; defaults to 0.05.
#' @param power The specified level of power.
#'
#' @return A list of the arguments (including the computed one).
#' @export
#'
#' @examples
#' pss.chisq.indep(pmatrix = matrix(c(0.050, 0.350, 0.100, 0.075, 0.250, 0.175),
#' nrow = 2, byrow = TRUE), N = 230)
#' pss.chisq.indep(pmatrix = matrix(c(0.3, 0.2, 0.4, 0.1), nrow = 2, byrow = TRUE), N = 200)

pss.chisq.indep <- function (pmatrix = NULL, N = NULL, alpha = 0.05, power = NULL) {

  # Check if the arguments are specified correctly
  if (sum(sapply(list(N, alpha, power), is.null)) != 1)
    stop("exactly one of N, alpha, and power must be NULL")
  if (is.null(pmatrix))
    stop("please specify the proportion matrix p")
  if (!is.null(N) && any(N < 1))
    stop("number of observations must be at least 1")

  # Calculate effect size and df
  pi <- apply(pmatrix, 1, sum)
  pj <- apply(pmatrix, 2, sum)
  p0 <- pi %*% t(pj)
  es <- sqrt(sum((pmatrix - p0)^2/p0))
  df <- (nrow(pmatrix) - 1) * (ncol(pmatrix) - 1)

  # Calculate test statistic
  p.body <- quote({
    stats::pchisq(stats::qchisq(alpha, df, lower = FALSE),
                  df, N * es^2, lower = FALSE)
  })

  # Use stats::uniroot function to calculate missing argument
  if (is.null(power))
    power <- eval(p.body)
  else if (is.null(N))
    N <- stats::uniroot(function(N) eval(p.body) - power, c(1 + 1e-10, 1e+09))$root
  else if (is.null(alpha))
    alpha <- stats::uniroot(function(sig.level) eval(p.body) - power, c(1e-10, 1 - 1e-10))$root
  else stop("internal error")

  # Generate output text
  METHOD <- "Chi-square test of independence power calculation"
  prows <- c()
  for (i in 1:nrow(pmatrix)) prows <- c(prows, paste(pmatrix[i,], collapse = ', '))

  # Print output as a power.htest object
  structure(list(pmatrix = paste(prows, collapse = " | "),
                 effect.size = es, N = N, alpha = alpha,
                 power = power, method = METHOD), class = "power.htest")
}
