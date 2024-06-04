#' One proportion sample size calculation
#'
#' Calculate sample size requirements for one proportion using the conditional
#' or unconditional method.
#'
#' @param p0 Null hypothesis proportion.
#' @param pA Alternative hypothesis proportion.
#' @param alpha Significance level.
#' @param power Power level as a decimal (0.8 for 80 percent power).
#' @param method Either "conditional" or "unconditional" method for calculation. Default is conditional.
#' @param one.or.two.sided Either "one" or "two" to specify a one or two sided hypothesis test. Default is two-sided.
#'
#' @return Returns n, the sample size needed for a one proportion test.
#' @export
#'
#' @keywords internal
#'
#' @examples
#' # Example 5.1: Calculate sample size needed for experimental therapy
#' # to be considered promising (at least 20%) if the true proportion of
#' # responders is 0.3. Using a one-sided test with significance at 0.05
#' # and power at 0.8.
#' oneprop_ss(p0 = 0.2, pA = 0.3, alpha = 0.05, power = 0.8,
#'            method = "conditional", one.or.two.sided = "one")
#' oneprop_ss(p0 = 0.2, pA = 0.3, alpha = 0.05, power = 0.8,
#'            method = "unconditional", one.or.two.sided = "one")
#' # Example 5.2: Same set up as above, but now the criteria are promising
#' # at 0.4 and true response is 0.5.
#' oneprop_ss(p0 = 0.4, pA = 0.5, alpha = 0.05, power = 0.8,
#'            method = "unconditional", one.or.two.sided = "one")
#' oneprop_ss(p0 = 0.4, pA = 0.5, alpha = 0.05, power = 0.8,
#'            method = "conditional", one.or.two.sided = "one")
oneprop_ss <- function(p0, pA, alpha, power, method = "conditional", one.or.two.sided = "two"){

  if(one.or.two.sided == "two"){
    sig_level <- alpha/2
  }
  if(one.or.two.sided == "one"){
    sig_level <- alpha
  }

  zb <- stats::qnorm(power)
  za <- stats::qnorm(1 - sig_level)

  if(method == "conditional"){
    n <- ceiling((zb*sqrt(pA*(1-pA)) + za*sqrt(p0*(1-p0)))^2 / (pA-p0)^2)
  }

  if(method == "unconditional"){
    n <- ceiling((zb+za)^2*pA*(1-pA)/(pA-p0)^2)
  }
  return(n)
}


