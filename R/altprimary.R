#' Power calculations for multiple alternative (at least one) primary continuous endpoints assuming known covariance matrix
#'
#' @description THIS IS A TEST DESCRIPTION.
#'
#' @param K The number of endpoints.
#' @param n1 The sample size for group 1.
#' @param n.ratio The ratio n2/n1 between the sample sizes of two groups; defaults to 1 (equal group sizes).
#' @param delta A vector of length K of the true mean differences mu1k - mu2k; must all be positive.
#' @param Sigma The covariance matrix of the K outcomes, of dimension K x K.
#' @param sd A vector of length K of the standard deviations of the K outcomes.
#' @param rho A vector of length 0.5*K*(K-1) of the correlations among the K outcomes.
#' @param alpha The significance level or type 1 error rate; defaults to 0.025. A one-sided test is assumed.
#' @param power The specified level of power.
#' @param tol The desired accuracy (convergence tolerance) for uniroot.
#' @param v Either TRUE for verbose output or FALSE to output computed argument only.
#'
#' @return A list of the arguments (including the computed one).
#' @import mvtnorm
#' @export
#'
#' @examples
#' altprimary(K = 2, n1 = 100, delta = c(0.4, 0.5), sd = c(1, 1), rho = 0.3,
#' alpha = 0.025, power = NULL)

altprimary <- function(K, n1 = NULL, n.ratio = 1, delta = NULL, Sigma, sd, rho,
                       alpha = 0.025, power = NULL,
                       tol = .Machine$double.eps^0.25, v = FALSE){

  # Check if the arguments are specified correctly
  check.many(list(n1, power), "oneof")
  check(n1, "pos")
  check(power, "unit")
  check(alpha, "req"); check(alpha, "unit")
  check(K, "req"); check(K, "min", min = 1); check(K, "int")
  check(n.ratio, "req"); check(n.ratio, "pos")
  check(v, "req"); check(v, "bool")

  check(delta, "req"); check(delta, "vec")
  if(length(delta) != K)
    stop("length of 'delta' must be equal to 'K'")
  if(!all(delta > 0))
    stop("all effect sizes need to be positive")

  if(!missing(Sigma)){ # Sigma is given
    check(Sigma, "mat")
    if(nrow(Sigma) != ncol(Sigma))
      stop("covariance matrix 'Sigma' must be quadratic")
    if(nrow(Sigma) != K)
      stop("covariance matrix must have dimension 'K' x 'K'")
    if(max(abs(Sigma - t(Sigma))) > 1e-10)
      stop("matrix 'Sigma' must be symmetric")
  }
  if(missing(Sigma)){ # make Sigma with sd & rho
    check(sd, "vec")
    if(missing(sd) || missing(rho))
      stop("if 'Sigma' is missing 'sd' and 'rho' must be given.")
    if(length(sd) != K)
      stop("length of 'sd' must be equal to 'K'")
    if(!all(sd > 0))
      stop("all standard deviations need to be positive")
    if(length(rho) != 0.5*K*(K-1))
      stop("length of 'rho' must be equal to '0.5*K*(K-1)'")
    if(!all(rho >= 0 & rho < 1))
      stop("all correlations need to be between 0 and 1")
    Sigma <- matrix(0, nrow = K, ncol = K)
    iter <- 0
    for(i in 1:(K-1)){
      for(j in (i+1):K){
        iter <- iter + 1
        Sigma[i,j] <- rho[iter]*sd[i]*sd[j]
      }
    }
    Sigma <- Sigma + t(Sigma)
    diag(Sigma) <- sd^2
  }
  if(!all(eigen(Sigma)$values > 0))
    stop("matrix 'Sigma' must be positive definite")
  Sigma.cor <- stats::cov2cor(Sigma)

  ## calculations

  if(is.null(power)){
    std.effect <- delta/sqrt(diag(Sigma))
    z.alpha <- stats::qnorm(1-alpha)
    crit.vals <- z.alpha - sqrt(n1*(n.ratio/(1+n.ratio)))*std.effect
    power <- 1-mvtnorm::pmvnorm(lower = -crit.vals, sigma = Sigma.cor)
    if (!v) return(power[1])
  }
  if(is.null(n1)){
    std.effect <- delta/sqrt(diag(Sigma))
    z.alpha <- stats::qnorm(1-alpha)
    ssize.fct <- function(n1, n.ratio, std.effect, z.alpha, Sigma.cor, power){
      crit.vals <- z.alpha - sqrt(n1*(n.ratio/(1+n.ratio)))*std.effect

      mvtnorm::pmvnorm(lower = -crit.vals, sigma = Sigma.cor) - (1- power)
    }
    n1 <- stats::uniroot(ssize.fct, c(2, 1e+05), tol = tol, extendInt = "yes",
                         n.ratio = n.ratio, std.effect = std.effect, z.alpha = z.alpha, Sigma.cor = Sigma.cor,
                         power = power)$root
    if (!v) return(n1)
  }
  n <- c(n1, n1*n.ratio)
  METHOD <- "Power calculation for alternative (at least one) primary endpoints"
  structure(list(n = n, n.ratio = n.ratio, delta = delta, sd = sqrt(diag(Sigma)),
                 rho = Sigma.cor[lower.tri(Sigma.cor)], Sigma = Sigma,
                 alpha = alpha, sides = 1,
                 power = power, method = METHOD),
            class = "power.htest")
}
