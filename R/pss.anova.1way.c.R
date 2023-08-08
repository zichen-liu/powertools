#' Power calculations for one-way balanced analysis of variance contrast test
#'
#' @param n The sample size per group.
#' @param means A vector of group means c(mu1, mu2, ...).
#' @param coeff A vector of contrast coefficients c(c1, c2, ...).
#' @param sd The estimated standard deviation within each group; defaults to 1.
#' @param alpha The significance level or type 1 error rate; defaults to 0.05.
#' @param power The specified level of power.
#'
#' @return A list of the arguments (including the computed one).
#' @export
#'
#' @examples
#' # Example 5.7
#' pss.anova.1way.c(n = 20, means = c(5, 10, 12), coeff = c(1, -1, 0), sd = 10, alpha = 0.025)
#' pss.anova.1way.c(n = 20, means = c(5, 10, 12), coeff = c(1, 0, -1), sd = 10, alpha = 0.025)

pss.anova.1way.c <- function (n = NULL, means = NULL, coeff = NULL, sd = 1,
                              alpha = 0.05, power = NULL) {

  # Check if the arguments are specified correctly
  groups <- length(means)
  if (sum(vapply(list(n, coeff, alpha, power), is.null, NA)) != 1)
    stop("exactly one of 'n', 'coeff', 'alpha', and 'power' must be NULL")
  if (groups < 2)
    stop("number of groups must be at least 2")
  if (!is.null(n) && n < 2)
    stop("number of observations in each group must be at least 2")
  if (groups != length(coeff))
    stop("number of contrast coefficients must be equal to the number of groups")
  if(is.null(sd))
    stop("sd must be specified")

  # Calculate df and ncp
  p.body <- quote({
    lambda <- coeff %*% means / (sd * sqrt(sum(coeff^2 / n)))
    df <- n * groups - groups
    stats::pf(q = stats::qf(alpha, 1, df, lower.tail = FALSE),
              1, df, lambda^2, lower.tail = FALSE)
  })

  # Use uniroot function to calculate missing argument
  if (is.null(power))
    power <- eval(p.body)
  else if (is.null(n))
    n <- uniroot(function(n) eval(p.body) - power, c(2, 1e+05))$root
  else if (is.null(alpha))
    alpha <- uniroot(function(alpha) eval(p.body) - power, c(1e-10, 1 - 1e-10))$root
  else stop("internal error", domain = NA)

  # Generate output text
  NOTE <- "n is the number in each group"
  METHOD <- "Balanced one-way analysis of variance\n     contrast test power calculation"

  # Print output as a power.htest object
  structure(list(groups = groups, n = n, means = means, coeff = coeff,
                 sd = sd, alpha = alpha, power = power,
                 note = NOTE, method = METHOD), class = "power.htest")
}
