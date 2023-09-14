#' Power calculations for a multiple linear regression partial F test
#'
#' @param n The sample size.
#' @param p The number of control predictors.
#' @param q The number of test predictors.
#' @param Rsq.red The squared sample multiple correlation coefficient in the reduced model.
#' @param Rsq.full The squared sample multiple correlation coefficient in the full model.
#' @param pc The partial correlation coefficient. Either both Rsq terms OR pc must be specified.
#' @param alpha The significance level or type 1 error rate; defaults to 0.05.
#' @param power The specified level of power.
#'
#' @return A list of the arguments (including the computed one).
#' @export
#'
#' @examples
#' # Example 10.10
#' pss.mlrF.partial(n = 80, p = 3, q = 2, Rsq.red = 0.25, Rsq.full = 0.35)
#' # Example 10.11
#' pss.mlrF.partial(n = 150, p = 4, pc = 0.2)

pss.mlrF.partial <- function (n = NULL, p = NULL, q = NULL, pc = NULL,
                              Rsq.red = NULL, Rsq.full = NULL,
                              alpha = 0.05, power = NULL) {

  # Check if the arguments are specified correctly
  if (sum(sapply(list(n, power, alpha), is.null)) != 1)
    stop("exactly one of n, alpha, and power must be NULL")
  if ((is.null(pc) & (is.null(p) | is.null(q))) | (!is.null(pc) & is.null(p)))
    stop("please specify the number of predictors")
  if ((is.null(Rsq.red) | is.null(Rsq.full)) & is.null(pc))
    stop("please specify Rsq.red and Rsq.full OR pc")
  if (!is.null(n) && any(n < 7))
    stop("number of observations must be at least 7")

  # Calculate power
  p.body <- quote({
    if (is.null(pc)) {
      ncp <- n * (Rsq.full - Rsq.red) / (1 - Rsq.full)
      df2 <- n - p - q - 1
      crit <- stats::qf(1 - alpha, q, df2)
      1 - stats::pf(crit, q, df2, ncp)
    } else {
      ncp <- n * pc^2 / (1 - pc^2)
      df2 <- n - p - 2
      crit <- stats::qf(1 - alpha, 1, df2)
      1 - stats::pf(crit, 1, df2, ncp)
    }
  })

  # Use uniroot function to calculate missing argument
  if (is.null(power))
    power <- eval(p.body)
  else if (is.null(n))
    n <- uniroot(function(n) eval(p.body) - power, c(7, 1e+09))$root
  else if (is.null(alpha))
    alpha <- uniroot(function(alpha) eval(p.body) - power, c(1e-10, 1 - 1e-10))$root
  else stop("internal error")

  # Generate output text
  METHOD <- "Power calculation for a multiple linear regression\n     partial F test"

  # Print output as a power.htest object
  if (is.null(pc)) {
    structure(list(n = n, p = p, q = q,
                   Rsq.red = Rsq.red, Rsq.full = Rsq.full,
                   alpha = alpha, power = power,
                   method = METHOD), class = "power.htest")
  } else {
    structure(list(n = n, p = p, pc = pc,
                   alpha = alpha, power = power,
                   method = METHOD), class = "power.htest")
  }

}
