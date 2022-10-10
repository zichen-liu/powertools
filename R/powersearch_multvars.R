#' Search through different variables to achieve adequate power
#'
#' @param searchvector1 One vector of variable values to search through. This is a required argument.
#' @param powerfunction The power function of interest. Any function can be specified.
#' @param searchvector2 Optional second vector of variable values to search through.
#' @param searchvector3 Optional third vector of variable values to search through.
#' @param searchvector4 Optional fourth vector of variable values to search through.
#' @param ... Additional arguments to pass to the specified powerfunction.
#'
#' @return A data frame with the resulting power at the specified search values.
#' @export
#'
#' @examples
#' # Example 1: Search through group sizes 3 and 4 (k) and group sample sizes 80
#' # and 100 (n) in one way ANOVA.
#' library(pwr2)
#' powersearch_multvars(searchvector1 = c(3, 3, 4, 4), searchvector2 = c(80, 100, 80, 100),
#'                    powerfunction = pwr.1way, searchvector3 = NULL, searchvector4= NULL,
#'                    alpha = 0.05, delta = 0.3, sigma = 8.5)
#'
#' # Example 2: Search through total number of observations 80, 100, or 120 (N) and
#' # allocation percents of 0.4 and 0.5 to group B (percent_B)
#' library(pwrAB)
#' powersearch_multvars(searchvector1 = c(80, 80, 100, 100, 120, 120), searchvector2 = c(0.4, 0.5),
#'                 powerfunction = AB_t2n, mean_diff = -2,
#'                 sd_A = 4 , sd_B = 6, sig_level = 0.05,
#'                 power = NULL, alt = c("less"))
powersearch_multvars <- function(searchvector1, powerfunction,
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


