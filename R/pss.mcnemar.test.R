#' # Example 6.10
#' pss.mcnemar.test(n = NULL, p1 = 0.8, p2 = 0.9, rho = 0.1, power = 0.9)
#' # Example 6.10
#' pss.mcnemar.test(n = NULL, paid = 0.08, psi = 0.18 / 0.08, power = 0.9)

pss.mcnemar.test <- function (n = NULL, p1 = NULL, p2 = NULL, rho = NULL,
                              paid = NULL, psi = NULL, alpha = 0.05,
          power = NULL, sided = c("two", "one")) {

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

  # Assign number of sides
  sided <- match.arg(sided)
  side <- switch(sided, one = 1, two = 2)

  # Calculate power
  f <- function(n, paid, psi, alpha, power) {
    bc <- ceiling(paid * n * (1 + psi))
    stats::pbinom(stats::qbinom(alpha / side, size = bc, prob = 0.5) - 1,
                  size = bc, prob = 1/(1 + psi)) + 1 -
                  pbinom(qbinom(1 - alpha/side, size = bc, prob = 0.5),
                  size = bc, prob = 1/(1 + psi))}
  p.body <- quote(stats::pnorm((sqrt(n * paid) * (psi - 1) -
                  stats::qnorm(alpha / side, lower.tail = FALSE) *
                  sqrt(psi + 1)) / sqrt((psi + 1) - paid * (psi - 1)^2)))

  # Use uniroot function to calculate missing argument
  if (is.null(power))
    power <- eval(p.body)
  else if (is.null(n))
    n <- uniroot(function(n) eval(p.body) - power, c(ceiling(log(alpha)/log(0.5)), 1e+07))$root
  else if (is.null(alpha))
    alpha <- uniroot(function(alpha) eval(p.body) - power, c(1e-10, 1 - 1e-10))$root
  else stop("internal error", domain = NA)
  NOTE <- "n is the number of pairs"
  METHOD <-"McNemar paired comparison of proportions\n     approximate power calculation"

  if (!is.null(p1) & !is.null(p2) & !is.null(rho)) {
    structure(list(n = n, p1 = p1, p2 = p2, rho = rho, alpha = alpha,
                   power = power, sided = sided, note = NOTE,
                   method = METHOD), class = "power.htest")
  }
  else if (!is.null(paid) & !is.null(psi))
    structure(list(n = n, paid = paid, psi = psi, alpha = alpha,
                   power = power, sided = sided, note = NOTE,
                   method = METHOD), class = "power.htest")
}
