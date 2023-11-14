#' Power calculations for test of equivalence of two proportions
#'
#' @param n1 The sample size for group 1.
#' @param n.ratio The ratio n2/n1 between the sample sizes of two groups; defaults to 1 (equal group sizes).
#' @param p1 The outcome proportion in group 1.
#' @param p2 The outcome proportion in group 2.
#' @param delta The equivalence margin.
#' @param alpha The significance level or type 1 error rate; defaults to 0.05.
#' @param power The specified level of power.
#' @param sides Either 1 or 2 (default) to specify a one- or two- sided hypothesis test.
#'
#' @return A list of the arguments (including the computed one).
#' @export
#'
#' @examples
#' pss.prop.test.equiv(n1 = NULL, p1 = 0.5, p2 = 0.5, delta = 0.1, alpha = 0.05, power = 0.8,
#' sides = 1)

pss.prop.test.equiv <- function (n1 = NULL, n.ratio = 1, p1 = NULL, p2 = NULL, delta = NULL,
                                 alpha = 0.05, power = NULL, sides = 2) {

  # Check if the arguments are specified correctly
  if (sides != 1 & sides != 2)
    stop("please specify 1 or 2 sides")
  if (sum(sapply(list(n1, n.ratio, power, alpha), is.null)) != 1)
    stop("exactly one of 'n1', 'n.ratio', 'alpha', and 'power' must be NULL")

  # Calculate test statistic
  p.body <- quote({
    d <- abs(p1 - p2)
    var <- p1 * (1 - p1) + p2 * (1 - p2)
    beta <- 1 - power
    ((stats::qnorm(1 - alpha / sides) + stats::qnorm(1 - beta / 2))^2 *
    var / n.ratio / (delta - d)^2)
  })

  # Use uniroot function to calculate missing argument
  if (is.null(n1))
    n1 <- eval(p.body)
  else if (is.null(power))
    power <- uniroot(function(power) eval(p.body) - n1, c(1e-05, 0.99999))$root
  else if (is.null(n.ratio))
    n.ratio <- uniroot(function(n.ratio) eval(p.body) - n1, c(2/n1, 1e+07))$root
  else if (is.null(alpha))
    alpha <- uniroot(function(alpha) eval(p.body) - n1, c(1e-10, 1 - 1e-10))$root
  else stop("internal error", domain = NA)

  # Generate output text
  NOTE <- "n is the number in each group"
  METHOD <- "Test for equivalence of two proportions power calculation"
  n <- c(n1, n1 * n.ratio)

  # Print output as a power.htest object
  structure(list(n = n, p1 = p1, p2 = p2,
                 delta = delta, alpha = alpha,
                 power = power, sides = sides, note = NOTE,
                 method = METHOD), class = "power.htest")
}

