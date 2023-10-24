#' Power calculations for one-way unbalanced analysis of variance omnibus F test
#'
#' @param nvec A vector of group sample sizes c(n1, n2, ...).
#' @param mvec A vector of group mvec c(mu1, mu2, ...).
#' @param sd The estimated standard deviation within each group.
#' @param Rsq The estimated R^2 for regressing the outcome on the covariates; defaults to 0.
#' @param ncov The number of covariates adjusted for in the model; defaults to 0.
#' @param alpha The significance level or type 1 error rate; defaults to 0.05.
#'
#' @return A list of the arguments (including the computed power).
#' @export
#'
#' @examples
#' pss.anova.unbal.1w(nvec = c(10, 20, 30), mvec = c(5, 10, 12), sd = 10)

pss.anova.unbal.1w <- function (nvec = NULL, mvec = NULL, sd = NULL,
                                Rsq = 0, ncov = 0, alpha = 0.05) {

  # Check if the arguments are specified correctly
  a <- length(mvec)
  if (a < 2)
    stop("number of a must be at least 2")
  if (any(nvec < 2))
    stop("number of observations in each group must be at least 2")
  if(is.null(nvec) | is.null(mvec))
    stop("sample size vector and means vector must both be specified")
  if(is.null(sd))
    stop("sd must be specified")
  if(a != length(nvec))
    stop("number of sample sizes must equal to the number of groups")

  # Get marginal mean
  es <- pss.anova.f.es(means = mvec, sd = sd)
  mmA <- es$mmA

  # Get f effect size
  f <- es$fA

  # Get ncp
  N <- sum(nvec)
  props <- nvec / N
  ws <- props %*% mmA
  temp <- sapply(X = 1:a, FUN = function(i) nvec[i] * ((mmA[i] - ws) / sd)^2)
  Lambda <- sum(temp) / (1 - Rsq)

  # Calculate power
  df1 <- a - 1
  df2 <- N - a - ncov
  power <- stats::pf(stats::qf(alpha, df1, df2, lower.tail = FALSE),
                     df1, df2, Lambda, lower.tail = FALSE)

  # Generate output text
  METHOD <- paste0("Unbalanced one-way analysis of ", ifelse(ncov < 1, "", "co"),
                   "variance\n     omnibus F test power calculation")
  out <- list(a = a, mvec = mvec, nvec = nvec, sd = sd,
              f = f, ncov = ncov, Rsq = Rsq,
              alpha = alpha, power = power,
              method = METHOD)

  # Print output as a power.htest object
  if (ncov < 1) out <- out[!names(out) %in% c("ncov", "Rsq")]
  structure(out, class = "power.htest")

}

