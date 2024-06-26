#' Power approximation for McNemar's test for two correlated proportions
#'
#' @param N The sample size; the number of pairs.
#' @param p1 The proportion in condition 1.
#' @param p2 The proportion in condition 2.
#' @param phi The estimated correlation between the two conditions.
#' @param paid The smaller of the two discordant probabilities. Either p1, p2, and phi, OR paid and dpr must be specified.
#' @param dpr The discordant proportion ratio. Either p1, p2, and phi, OR paid and dpr must be specified.
#' @param alpha The significance level or type 1 error rate; defaults to 0.05.
#' @param power The specified level of power.
#' @param sides Either 1 or 2 (default) to specify a one- or two- sided hypothesis test.
#' @param v Either TRUE for verbose output or FALSE to output computed argument only.
#'
#' @return A list of the arguments (including the computed one).
#' @export
#'
#' @examples
#' prop.paired(N = NULL, p1 = 0.8, p2 = 0.9, phi = 0, power = 0.9, sides = 2)
#' prop.paired(N = NULL, paid = 0.08, dpr = 0.18 / 0.08, power = 0.9, sides = 2)

prop.paired <- function (N = NULL, p1 = NULL, p2 = NULL, phi = NULL,
                         paid = NULL, dpr = NULL, alpha = 0.05,
                         power = NULL, sides = 2, v = FALSE) {

  # Check if the arguments are specified correctly
  if ((is.null(p1) | is.null(p2) | is.null(phi)) & (is.null(dpr) | is.null(paid)))
    stop("p1, p2, and phi OR dpr and paid must be specified")
  check.many(list(N, alpha, power), "oneof")
  check(N, "pos")
  check(p1, "unit")
  check(p2, "unit")
  check(paid, "unit")
  check(dpr, "pos")
  check(alpha, "unit")
  check(power, "unit")
  check(sides, "req"); check(sides, "vals", valslist = c(1, 2))
  check(v, "req"); check(v, "bool")

  # Calculate paid and dpr if not given
  if (is.null(paid) & is.null(dpr)) {
    p01 <- p1 * (1 - p2) - phi * sqrt((1 - p2) * p1 * (1 - p1) * p2)
    p10 <- p01 + p2 - p1
    paid <- p01
    dpr <- p10 / p01
  }

  # Calculate test statistic
  p.body <- quote(stats::pnorm((sqrt(N * paid) * (dpr - 1) -
                  stats::qnorm(alpha / sides, lower.tail = FALSE) *
                  sqrt(dpr + 1)) / sqrt((dpr + 1) - paid * (dpr - 1)^2)))

  NOTE <- "N is the number of pairs"
  if (!v) print(paste("NOTE:", NOTE))

  # Use stats::uniroot function to calculate missing argument
  if (is.null(power)) {
    power <- eval(p.body)
    if (!v) return(power)
  }
  else if (is.null(N)) {
    N <- stats::uniroot(function(N) eval(p.body) - power, c(ceiling(log(alpha) / log(0.5)), 1e+07))$root
    if (!v) return(N)
  }
  else if (is.null(alpha)) {
    alpha <- stats::uniroot(function(alpha) eval(p.body) - power, c(1e-10, 1 - 1e-10))$root
    if (!v) return(alpha)
  }
  else stop("internal error", domain = NA)

  # Generate output text
  METHOD <-"McNemar paired comparison of proportions\n     approximate power calculation"

  # Print output as a power.htest object depending on which inputs were given
  if (!is.null(p1) & !is.null(p2) & !is.null(phi)) {
    p <- c(p1, p2)
    structure(list(N = N, `p1, p2` = p, phi = phi, alpha = alpha,
                   power = power, sides = sides, note = NOTE,
                   method = METHOD), class = "power.htest")
  }
  else if (!is.null(paid) & !is.null(dpr))
    structure(list(N = N, paid = paid, dpr = dpr, alpha = alpha,
                   power = power, sides = sides, note = NOTE,
                   method = METHOD), class = "power.htest")
}

