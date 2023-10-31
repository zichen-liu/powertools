#' Power calculations for relative risk
#'
#' @param n The sample size for group 1.
#' @param n.ratio The ratio n2/n1 between the sample sizes of two groups; defaults to 1 (equal group sizes).
#' @param p1 The proportion in group 1.
#' @param p2 The proportion in group 2.
#' @param RR0 The relative risk under the null; defaults to 1.
#' @param alpha The significance level or type 1 error rate; defaults to 0.05.
#' @param power The specified level of power.
#' @param sides Either 1 or 2 (default) to specify a one- or two- sided hypothesis test.
#'
#' @return A list of the arguments (including the computed one).
#' @export
#'
#' @examples
#' pss.rr(n = NULL, n.ratio = 6, p1 = 0.1, p2 = 0.1 * 2, power = 0.8)

pss.rr <- function (n = NULL, n.ratio = 1, p1 = NULL, p2 = NULL, RR0 = 1,
                    alpha = 0.05, power = NULL, sides = 2) {

  # Check if the arguments are specified correctly
  if (sides != 1 & sides != 2)
    stop("please specify 1 or 2 sides")
  if (sum(sapply(list(n, n.ratio, power, alpha), is.null)) != 1)
    stop("exactly one of 'n', 'n.ratio', 'alpha', and 'power' must be NULL")
  if (!is.null(n.ratio) && n.ratio <= 0)
    stop("n.ratio between group sizes must be positive")

  # Calculate test statistic
  p.body <- quote({
    RR <- p2 / p1
    d <- abs(log(RR) - log(RR0))
    q1 <- 1 - p1
    q2 <- 1 - p2
    denom <- (1 / n) * (q1 / p1 + q2 / (n.ratio * p2))
    (stats::qnorm(alpha / sides) + d / sqrt(denom))
  })

  # Use uniroot function to calculate missing argument
  if (is.null(power))
    power <- eval(p.body)
  else if (is.null(n))
    n <- uniroot(function(n) eval(p.body) - power, c(2, 1e+07))$root
  else if (is.null(n.ratio))
    n.ratio <- uniroot(function(n.ratio) eval(p.body) - power, c(2/n, 1e+07))$root
  else if (is.null(alpha))
    alpha <- uniroot(function(alpha) eval(p.body) - power, c(1e-10, 1 - 1e-10))$root
  else stop("internal error")

  # Generate output text
  NOTE <- "n is the number in each group"
  METHOD <- "Relative risk power calculation"
  n <- c(n, n * n.ratio)
  p <- c(p1, p2)
  RR <- c(p2 / p1, RR0)

  print(n.ratio)

  # Print output as a power.htest object
  structure(list(n = n, `p1, p2` = p, `RR, RR0` = RR,
                 alpha = alpha, power = power, sides = sides, note = NOTE,
                 method = METHOD), class = "power.htest")
}

