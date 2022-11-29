#' Search through multiple variables to achieve adequate power
#'
#' This function will search through any power function with up to 4 variables to search through. Because the search vectors are supplied without names, they must be provided in the default order of the power function.
#'
#' @param powerfunction The power function of interest. Any function can be specified.
#' @param searchvector1 One vector of variable values to search through. This is a required argument.
#' @param searchvector2 Optional second vector of variable values to search through. This must either be the same length as the searchvector1 or can be a multiple of searchvector1. If searchvector1 has 4 elements then searchvector2 needs to have either 2, 4, 8, etc elements.
#' @param searchvector3 Optional third vector of variable values to search through. Similarly to searchvector2, this must either be the same length as other searchvectors or a multiple.
#' @param searchvector4 Optional fourth vector of variable values to search through. This must either be the same length as other searchvectors or a multiple.
#' @param ... Additional arguments to pass to the specified powerfunction.
#'
#' @return A data frame with the resulting power at the specified search values.
#' @export
#'
#' @examples
#' # Example 1: Search through group sizes 3 and 4 (k) and group sample sizes 80
#' # and 100 (n) in one way ANOVA. alpha = 0.05, delta = 3, sigma = 8.5 are all
#' # arguments that are being passed to the pwr.1way function. These arguments are
#' # required for this function.
#' library(pwr2)
#' search_power(powerfunction = pwr.1way, searchvector1 = c(3, 3, 4, 4),
#'                      searchvector2 = c(80, 100, 80, 100),
#'                      searchvector3 = NULL, searchvector4= NULL,
#'                      alpha = 0.05, delta = 3, sigma = 8.5)
#'
#' # Example 2: Search through total number of observations 80, 100, or 120 (N) and
#' # allocation percents of 0.4 and 0.5 to group B (percent_B)
#' library(pwrAB)
#' search_power(powerfunction = AB_t2n, searchvector1 = c(80, 80, 100, 100, 120, 120),
#'                      searchvector2 = c(0.4, 0.5), mean_diff = -2,
#'                      sd_A = 4 , sd_B = 6, sig_level = 0.05,
#'                      power = NULL, alt = c("less"))
search_power <- function(powerfunction, searchvector1,
                               searchvector2 = NULL, searchvector3 = NULL,
                               searchvector4 = NULL, ...){

  if(!is.null(searchvector2) & !is.null(searchvector3) & !is.null(searchvector4)){
    out <- mapply(searchvector1, searchvector2, searchvector3, searchvector4,
           FUN = function(var1, var2, var3, var4)
             powerfunction(var1, var2, var3, var4, ...)) %>% t() %>% as.data.frame()
  }

  if(!is.null(searchvector2) & !is.null(searchvector3) & is.null(searchvector4)){
    out <- mapply(searchvector1, searchvector2, searchvector3,
           FUN = function(var1, var2, var3)
             powerfunction(var1, var2, var3, ...)) %>% t() %>% as.data.frame()
  }

  if(!is.null(searchvector2) & is.null(searchvector3) & is.null(searchvector4)){
    out <- mapply(searchvector1, searchvector2,
           FUN = function(var1, var2)
             powerfunction(var1, var2, ...)) %>% t() %>% as.data.frame()

  }

  if(is.null(searchvector2) & is.null(searchvector3) & is.null(searchvector4)){
   out <- mapply(searchvector1, FUN = function(var1)
             powerfunction(var1,  ...)) %>% t() %>% as.data.frame()
  }
  return(out)
}


