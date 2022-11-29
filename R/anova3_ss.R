#' One way ANOVA sample size calculation for omnibus F test with 3 group levels
#'
#' This function searches through different group sample sizes until the specified power level is achieved.
#' It allows for unequal standard deviations between groups and unequal group allocation.
#' This function is a wrapper for the 'anova1f_3' function in the 'pwr2ppl' package.
#'
#' @param means A vector of group means.
#' @param sds A vector of group standard deviations. This function allows unequal standard deviations between groups.
#' @param alphalevel The significance level.
#' @param powerlevel The level of power to achieve as a decimal.
#' @param groupratios A 3 element vector of allocation ratios to each group. Default is a 1:1:1 ratio, supplied as c(1,1,1).
#' @param Nmin The minimum total sample size N to search through. Default start is 30.
#' @param Nmax The maximum total sample size N to search through. Default ending is 300.
#' @param upperbound The largest integer to multiply the allocation ratios by. Default is 500 meaning if there is even allocation, this will allow a search of up to 500 people per group.
#'
#' @return Returns a data frame with sample sizes for each group indicated in n1, n2, and n3, actual power, and any notes from the computation.
#' @export
#'
#' @examples
#' # Sample Size for 3 groups with unequal allocation at 0.05 significance and 80% power
#' # explores range itself from N = 30 to 300
#' anova3_ss(means = c(5, 10, 12), sds = c(10, 8, 10), alphalevel = 0.05,
#' powerlevel = 0.8, groupratios = c(1, 1, 3),
#' Nmin = NULL, Nmax = NULL)
anova3_ss <- function(means, sds, alphalevel, powerlevel, groupratios = NULL,
                              Nmin = NULL, Nmax = NULL, upperbound = 500){
  # If no values supplied, use default values.
  # minimum of 30 ppl total, max of 300, with equal allocation.
  if(is.null(Nmin)){
    Nmin <- 30
  }
  if(is.null(Nmax)){
    Nmax <- 300
  }
  if(is.null(groupratios)){
    groupratios <- c(1,1,1)
  }

  # for exact ratio search
  i_vec <- seq(5, upperbound, 1)
  i_vec_mat <- matrix(i_vec, nrow=length(i_vec), ncol = 1)
  ratio_mat <- matrix(groupratios, nrow=1, ncol=length(groupratios))

  N_vecs <- i_vec_mat %*% ratio_mat
  Ntot <- rowSums(N_vecs)
  N_vec_4search <- N_vecs[Ntot >= Nmin & Ntot <= Nmax, ]
  # now this is the N vector to loop through

  # If this search vector is too long, then do 5 rough searches initially and then do more steps.
  N_rough_search <- NULL
  if(dim(N_vec_4search)[1] > 20){
    Nrows <- dim(N_vec_4search)[1]
    increment <- floor(Nrows / 5)
    N_rough_search <- N_vec_4search[c(1, increment * c(2,3,4), Nrows), ]
  }

  row_start <- 1

  ## Rough search OR full Search when small enough ##
  if(!is.null(N_rough_search)){
  for(i in 1:dim(N_rough_search)[1]){

    n1 = N_rough_search[i, 1]
    n2 = N_rough_search[i, 2]
    n3 = N_rough_search[i, 3]

    an_power <- purrr::quietly(pwr2ppl::anova1f_3)(m1 = means[1], m2 = means[2], m3 = means[3],
                          s1 = sds[1], s2 = sds[2], s3 = sds[3],
                          n1 = n1, n2 = n2, n3 = n3, alpha = alphalevel)$result

    if(an_power$Power >= powerlevel){ # if power is sufficient
      if(i == i){# AND if this is first value then last calc is first calculation
        last_calc <- c(n1, n2, n3, an_power$Power, i)
      }
      break # end for loop.
    }

    # this will be from previous calculation because it will break if met
    last_calc <- c(n1, n2, n3, an_power$Power, i)
  }
    row_start <- which(last_calc[1] == N_vec_4search[, 1])
    } # end of rough search


  ## rest of search start at the last n values or at 1 if we did not do prior loop.
  N_search <- matrix(N_vec_4search[row_start:dim(N_vec_4search)[1], ], ncol = 3)

  for(i in 1:dim(N_search)[1]){

    n1 <- N_search[i, 1]
    n2 <- N_search[i, 2]
    n3 <- N_search[i, 3]

    if(n1 + n2 + n3 >= Nmax){ # more ppl than max
      Note <- "Larger Sample Size Required"
      if(i == 1){
        # use previous at max
        Final <- data.frame(n1 = an_power$n1, n2 = an_power$n2, n3 = an_power$n3,
                          Power = an_power$Power, Note = Note)
        break
      }
      if(i != 1){
        # use final value from this iteration
       Final <- data.frame(n1 = an_power2$n1, n2 = an_power2$n2, n3 = an_power2$n3,
                          Power = an_power2$Power, Note = Note)
        break
      }
    }

    an_power2 <- purrr::quietly(pwr2ppl::anova1f_3)(m1 = means[1], m2 = means[2], m3 = means[3],
                          s1 = sds[1], s2 = sds[2], s3 = sds[3],
                          n1 = n1, n2 = n2, n3 = n3, alpha = alphalevel)$result

    if(an_power2$Power >= powerlevel){
      # stop loop if power is sufficient
      Note <- ifelse(n1 + n2 + n3 > Nmax, "Larger Sample Size Required", " ")
      Final <- data.frame(n1 = n1, n2 = n2, n3 = n3, Power = an_power2$Power,
                          Note = Note)
      break
    }

  }
  return(Final)
}
