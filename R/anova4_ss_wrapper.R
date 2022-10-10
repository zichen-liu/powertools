#' One way ANOVA sample size calculation for omnibus F test with 4 group levels
#'
#' This function searches through different group sample sizes until the specified power level is achieved.
#' It allows for unequal standard deviations between groups and unequal group allocation.
#'
#' @param means A vector of group means.
#' @param sds A vector of group standard deviations. This function allows unequal standard deviations between groups.
#' @param alphalevel The significance level.
#' @param powerlevel The level of power to achieve.
#' @param groupweights A vector of allocation to each group in terms of decimals. This needs to sum to 1.
#' @param Nmin The minimum sample size N to search through. If not supplied the function starts at 30.
#' @param Nmax The maximum sample size N to search through. If not supplied the function ends at 150.
#' @param Nchange The initial sample size step to use in the rough search. If not supplied the function searches by increments of 30 people.
#'
#' @return Returns a data frame with sample sizes for each group indicated in n1, n2, n3, and n4, actual power, and any notes from the computation.
#' @export
#'
#' @examples
#' # explores range
#' anova4_ss_wrapper(means = c(5, 10, 12, 15), sds = c(10, 10, 10, 10), alphalevel = 0.05,
#'                   powerlevel = 0.8, groupweights = c(.3, .2, .2, .3),
#'                   Nmin = 60, Nmax = 120, Nchange = 20)

anova4_ss_wrapper <- function(means, sds, alphalevel, powerlevel, groupweights,
                              Nmin = NULL, Nmax = NULL, Nchange = NULL){
  # make sequence to explore by
  if(is.null(Nmin)){
    Nmin <- 30
  }
  if(is.null(Nmax)){
    Nmax <- 150
  }
  if(is.null(Nchange)){
    Nchange <- 30
  }

  N_vec_first <- seq(Nmin, Nmax, Nchange)

  if(N_vec_first[length(N_vec_first)] != Nmax){
    N_vec_first <- c(N_vec_first, Nmax)
  }

  for(i in N_vec_first){
    # start with smallest N
    n1 = max(ceiling(groupweights[1] * i), 5) # minimum of 5 ppl per group
    n2 = max(ceiling(groupweights[2] * i), 5)
    n3 = max(ceiling(groupweights[3] * i), 5)
    n4 = max(ceiling(groupweights[4] * i), 5)

    an_power <- purrr::quietly(pwr2ppl::anova1f_4)(m1 = means[1], m2 = means[2], m3 = means[3], m4 = means[4],
                          s1 = sds[1], s2 = sds[2], s3 = sds[3], s4 = sds[4],
                          n1 = n1, n2 = n2, n3 = n3, n4 = n4, alpha = alphalevel)$result

    if(an_power$Power >= powerlevel){
      if(i == N_vec_first[1]){ # if this is first value and power is sufficient
        last_calc <- c(n1, n2, n3, n4, an_power$Power, i)
      }
      break
    }
    # this will be from previous calculation because it will break if met
    last_calc <- c(n1, n2, n3, n4, an_power$Power, i)
  }

  # update the n's to be the last n's that weren't enough
    n1_last <- last_calc[1]; n2_last <- last_calc[2];
    n3_last <- last_calc[3]; n4_last <- last_calc[4];
  for(i in 1:Nchange){
    # now go up by 1 each time but keeping the ratio proper
    n1 <- ceiling(n1_last + i * groupweights[1])
    n2 <- ceiling(n2_last + i * groupweights[2])
    n3 <- ceiling(n3_last + i * groupweights[3])
    n4 <- ceiling(n4_last + i * groupweights[4])

    if(n1 + n2 + n3 + n4 >= Nmax){ # more ppl than max
      Note <- "Larger Sample Size Required"
      if(i == 1){
        # use previous at max
        Final <- data.frame(n1 = an_power$n1, n2 = an_power$n2, n3 = an_power$n3,
                            n4 = an_power$n4, Power = an_power$Power, Note = Note)
        break
      }
      if(i != 1){
        # use final value from this iteration
       Final <- data.frame(n1 = an_power2$n1, n2 = an_power2$n2, n3 = an_power2$n3,
                           n4 = an_power2$n4, Power = an_power2$Power, Note = Note)
        break
      }
    }

    an_power2 <- purrr::quietly(pwr2ppl::anova1f_4)(m1 = means[1], m2 = means[2], m3 = means[3], m4 = means[4],
                          s1 = sds[1], s2 = sds[2], s3 = sds[3], s4 = sds[4],
                          n1 = n1, n2 = n2, n3 = n3, n4 = n4, alpha = alphalevel)$result

    if(an_power2$Power >= powerlevel){
      # stop loop if power is sufficient
      Note <- ifelse(n1 + n2 + n3 + n4 > Nmax, "Larger Sample Size Required", " ")
      Final <- data.frame(n1 = n1, n2 = n2, n3 = n3, n4 = n4,
                          Power = an_power2$Power, Note = Note)
      break
    }

  }
  return(Final)
}
