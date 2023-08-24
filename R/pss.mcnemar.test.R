#' Power approximation for McNemar's test for two correlated proportions
#'
#' @param n The sample size; the number of pairs.
#' @param p1 The proportion in condition 1.
#' @param p2 The proportion in condition 2.
#' @param rho The estimated correlation between the two conditions.
#' @param paid The smaller of the two discordant probabilities. Either p1, p2, and rho, OR paid and psi must be specified.
#' @param psi The discordant proportion ratio. Either p1, p2, and rho, OR paid and psi must be specified.
#' @param alpha The significance level or type 1 error rate; defaults to 0.05.
#' @param power The specified level of power.
#' @param sides Either 1 or 2 (default) to specify a one- or two- sided hypothesis test.
#'
#' @return A list of the arguments (including the computed one).
#' @export
#'
#' @examples
#' # Example 6.10
#' pss.mcnemar.test(n = NULL, p1 = 0.8, p2 = 0.9, rho = 0, power = 0.9, sides = 2)
#' # Example 6.10
#' pss.mcnemar.test(n = NULL, paid = 0.08, psi = 0.18 / 0.08, power = 0.9, sides = 2)

pss.mcnemar.test <- function (n = NULL, p1 = NULL, p2 = NULL, rho = NULL,
                              paid = NULL, psi = NULL, alpha = 0.05,
                              power = NULL, sides = c(2, 1)) {

  # Check if the arguments are specified correctly
  if (sum(sapply(list(n, alpha, power), is.null)) !=  1) {
    stop("exactly one of 'n', 'alpha', and 'power' must be NULL")
  }
  if ((is.null(p1) | is.null(p2) | is.null(rho)) & (is.null(psi) | is.null(paid)))
    stop("p1, p2, and rho OR psi and paid must be specified")

  # Calculate paid and psi if not given
  if (is.null(paid) & is.null(psi)) {
    p01 <- p1 * (1 - p2) - rho * sqrt((1 - p2) * p1 * (1 - p1) * p2)
    p10 <- p01 + p2 - p1
    paid <- p01
    psi <- p10 / p01
  }

  # Calculate test statistic
  p.body <- quote(stats::pnorm((sqrt(n * paid) * (psi - 1) -
                  stats::qnorm(alpha / sides, lower.tail = FALSE) *
                  sqrt(psi + 1)) / sqrt((psi + 1) - paid * (psi - 1)^2)))

  # Use uniroot function to calculate missing argument
  if (is.null(power))
    power <- eval(p.body)
  else if (is.null(n))
    n <- uniroot(function(n) eval(p.body) - power, c(ceiling(log(alpha) / log(0.5)), 1e+07))$root
  else if (is.null(alpha))
    alpha <- uniroot(function(alpha) eval(p.body) - power, c(1e-10, 1 - 1e-10))$root
  else stop("internal error", domain = NA)

  # Generate output text
  NOTE <- "n is the number of pairs"
  METHOD <-"McNemar paired comparison of proportions\n     approximate power calculation"

  # Print output as a power.htest object depending on which inputs were given
  if (!is.null(p1) & !is.null(p2) & !is.null(rho)) {
    structure(list(n = n, p1 = p1, p2 = p2, rho = rho, alpha = alpha,
                   power = power, sides = sides, note = NOTE,
                   method = METHOD), class = "power.htest")
  }
  else if (!is.null(paid) & !is.null(psi))
    structure(list(n = n, paid = paid, psi = psi, alpha = alpha,
                   power = power, sides = sides, note = NOTE,
                   method = METHOD), class = "power.htest")
}

