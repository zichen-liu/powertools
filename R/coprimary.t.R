#' Power calculations for multiple co-primary continuous endpoints assuming unknown covariance matrix
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
#' @param M The number of simulation.
#' @param min.n Minimum value of n1; used in search for n1 to achieve desired power.
#' @param max.n Maximum value of n1; used in search for n1 to achieve desired power.
#' @param tol The desired accuracy (convergence tolerance) for uniroot.
#' @param use.uniroot Whether to use the uniroot function to calculate n1; defaults to TRUE.
#' @param v Either TRUE for verbose output or FALSE to output computed argument only.
#'
#' @return A list of the arguments (including the computed one).
#' @import mvtnorm
#' @export
#'
#' @examples
#' coprimary.t(K = 2, n1 = 100, delta = c(0.4, 0.5), sd = c(1, 1), rho = 0.3,
#' alpha = 0.025, power = NULL)

coprimary.t <- function(K, n1 = NULL, n.ratio = 1, delta = NULL, Sigma, sd, rho,
                        alpha = 0.025, power = NULL, M = 10000, min.n = NULL, max.n = NULL,
                        tol = .Machine$double.eps^0.25, use.uniroot = TRUE, v = FALSE) {

  ## check of input

  if(missing(K))
    stop("specify the number of co-primary endpoints")
  if(!is.numeric(K))
    stop("'K' must be a natural number > 1")
  K <- as.integer(K)
  if(is.null(n1) & is.null(power))
    stop("either 'n1' or 'power' must be specified")
  if(!is.null(n1) & !is.null(power))
    stop("either 'n1' or 'power' must be NULL")
  if(!is.null(n1)){
    if(length(n1) > 1){
      warning("length of 'n1' is greater than 1, only the first entry is used")
      n1 <- n1[1]
    }
    n1 <- as.integer(n1)
  }
  if(is.null(n.ratio)){
    stop("n.ratio cannot be NULL")}
  if(!n.ratio > 0){
    stop("n.ratio must be positive")
  }
  if(!is.null(power)){
    if(power <= 0 | power >= 1)
      stop("power must be in (0, 1)")
  }
  if(is.null(delta))
    stop("expected effect size 'delta' is missing")
  if(length(delta) < 2)
    stop("length of 'delta' < 2: effect for at least two co-primary endpoints is required")
  if(!all(delta > 0))
    stop("all effect sizes need to be positive")
  if(!missing(Sigma)){
    if(nrow(Sigma) != ncol(Sigma))
      stop("covariance matrix 'Sigma' must be square")
    if(nrow(Sigma) != K)
      stop("covariance matrix must have dimension 'K' x 'K'")
    if(max(abs(Sigma - t(Sigma))) > 1e-10)
      stop("matrix 'Sigma' must be symmetric")
  }
  if(missing(Sigma)){
    if(missing(sd) || missing(rho))
      stop("if 'Sigma' is missing 'sd' and 'rho' must be given.")
    if(length(sd) != K)
      stop("length of 'sd' must be equal to 'K'")
    if(length(rho) != 0.5*K*(K-1))
      stop("length of 'rho' must be equal to '0.5*K*(K-1)'")
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

  if(alpha <= 0 | alpha >= 1)
    stop("significance level must be in (0, 1)")

  if(length(M) > 1){
    warning("length of 'M' is greater than 1, only the first entry is used")
    M <- M[1]
  }
  M <- as.integer(M)
  if(is.null(n1)){
    if(length(min.n) > 1){
      warning("length of 'min.n' is greater than 1, only the first entry is used")
      min.n <- min.n[1]
    }
    min.n <- as.integer(min.n)
    if(length(max.n) > 1){
      warning("length of 'max.n' is greater than 1, only the first entry is used")
      max.n <- max.n[1]
    }
    max.n <- as.integer(max.n)
    if(min.n < 4)
      stop("'min.n' must be >= 4")
    if(min.n >= max.n)
      stop("'min.n' must be < 'max.n'")
  }
  check(v, "req"); check(v, "bool")
  ## calculations

  if(is.null(power)){
    std.effect <- delta/sqrt(diag(Sigma))
    probs <- numeric(M)
    Ws <- stats::rWishart(M, df = n1*(n.ratio+1)-2, Sigma = Sigma.cor)
    for(i in 1:M){
      Wi <- diag(Ws[,,i])
      ci <- stats::qt(1-alpha, df = n1*(n.ratio+1)-2)*sqrt(Wi/(n1*(n.ratio+1)-2)) - sqrt(n1*(n.ratio/(1+n.ratio)))*std.effect
      probs[i] <- mvtnorm::pmvnorm(upper = -ci, sigma = Sigma.cor)
    }
    power <- mean(probs)
    if (!v) return(power)
  }

  if(is.null(n1)){
    std.effect <- delta/sqrt(diag(Sigma))
    ssize.fct <- function(n1, n.ratio, std.effect, Sigma.cor, power, M, verbose = FALSE){
      probs <- numeric(M)
      Ws <- stats::rWishart(M, df = n1*(n.ratio+1)-2, Sigma = Sigma.cor)
      for(i in 1:M){
        Wi <- diag(Ws[,,i])
        ci <- stats::qt(1-alpha, df = n1*(n.ratio+1)-2)*sqrt(Wi/(n1*(n.ratio+1)-2)) - sqrt(n1*(n.ratio/(1+n.ratio)))*std.effect
        probs[i] <- mvtnorm::pmvnorm(upper = -ci, sigma = Sigma.cor) - power
      }
      if(verbose) cat("Current precision:\t", mean(probs), "\n")
      mean(probs)
    }
    if(use.uniroot){
      n1 <- stats::uniroot(ssize.fct, c(min.n, max.n), tol = tol, extendInt = "yes", n.ratio = n.ratio,
                    std.effect = std.effect, Sigma.cor = Sigma.cor, power = power, M = M,
                    verbose = TRUE)$root
    }else{
      ns <- min.n:max.n
      res <- sapply(ns, ssize.fct, std.effect = std.effect, Sigma.cor = Sigma.cor,
                    power = power, M = M)
      ns.pos <- ns[res > 0]
      res.pos <- res[res > 0]
      cat("Precision:\t", min(res.pos), "\n")
      n1 <- ns.pos[which.min(res.pos)]
    }

    if (!v) return(n1)
  }
  n <- c(n1, n1/n.ratio)
  METHOD <- "Power calculation for multiple co-primary endpoints (covariance matrix unknown)"
  structure(list(n = n, n.ratio = n.ratio, delta = delta, sd = sqrt(diag(Sigma)),
                 rho = Sigma.cor[lower.tri(Sigma.cor)], Sigma = Sigma,
                 alpha = alpha, sides = 1,
                 power = power, method = METHOD), class = "power.htest")
}

