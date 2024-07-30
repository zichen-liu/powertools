#' Additive Treatments Stepped Wedge Design Power Calculation
#'
#' Calculate power for detecting treatment effects in a two-treatment stepped wedge design with additive treatment effects.
#' We only consider complete designs. This function is based on Power Analysis for Stepped Wedge Trials with Two Treatments published in Stat Med in 2022 by Phillip Sundin and Catherine Crespi. This manuscript uses the normal distribution and the Wald Test to derive power.
#'
#' @param RhoW The intraclass correlation coefficient, which refers to the correlation between outcomes of two different individuals in the same cluster at the same time.
#' @param ModelChoice A string indicating either a repeated cross sectional design ("RCS"), nested exchangeable ("NEC"), Cohort ("Cohort"), or any one of these three designs in their standardized form. Only one of "RCS Standardized", "RCS", "Cohort Standardized", "Cohort", "NEC Standardized", or "NEC" can be specified. The standardized design refers to dividing by the standard deviation of y (the outcome) in the covariance specification.
#' @param n.individuals The number of individuals in each condition at one time period.
#' @param n.clusters The number of clusters
#' @param n.periods The number of time periods.
#' @param delta_1 The standardized effect size for treatment 1
#' @param delta_2 The standardized effect size for treatment 2
#' @param RhoA The correlation between two observations in the same individual in the same cluster but different time periods. RCS models assume RhoW = RhoA, and Cohort models assume RhoA >= RhoW.
#' @param IAC The individual auto-correlation, which is the proportion of the individual-level variance that is time-invariant. IAC = 0 in a RCS.
#' @param Sequence_Tx1 A vector representing the time points (periods) at which each cluster transitions from control to treatment 1 (ie the first occurence).
#' @param Sequence_Tx2 A vector representing the time points (periods) at which each cluster transitions from control to treatment 2 (ie the first occurence).
#' @param typeIerror The significance level.
#'
#' @return A list with power for detecting treatment effect for treatment 1 and 2
#' @export
#'
#' @examples
#' ## For example, the 3 cluster, 4 period SWD,
#' ##     0 1 1   1
#' ##     0 0 1+2 1+2
#' ##     0 0 2   1+2
#' ## where 0 denotes control condition, 1 denotes the condition with only
#' ## treatment 1, 2 denote the condition with only treatment 2, and 1+2 denotes the condition
#' ## where a cluster receives both treatment 1 and 2
#' ## would have Sequence1 <- c(2,3,4) and Sequence2 <- c(NA,3,3)
#' ## notice sequence 1 cluster 3 has time period 4 even though time period 3 had treatment 2,
#' ## so it is the first occurrence of trt 1 in this cluster.
#'
#'
#' # 12-cluster concurrent repeated cross sectional design
#' # Sequencing has 2 clusters transition to treatment 1 at time 2,
#' # 2 clusters transition to treatment 1 at time 3, and 2 clusters transition
#' # to treatment 1 at time 3. These clusters never receive treatment 2
#' Sequence_Tx1 <- c(2, 2, 3, 3, 4, 4, NA, NA, NA, NA, NA, NA)
#' # similar sequencing for treatment 2
#' Sequence_Tx2 <- c(NA, NA, NA, NA, NA, NA, 4, 4, 3, 3, 2, 2)
#' swd_2trt_additive_power(RhoW = 0.05, ModelChoice = "RCS", RhoA = 0.2, IAC = 0,
#'                n.individuals = 30, n.clusters = 12, n.periods = 4, delta_1 = 0.4, delta_2 = 0.1,
#'                Sequence_Tx1 = Sequence_Tx1, Sequence_Tx2 = Sequence_Tx2,
#'                typeIerror = 0.05)
swd_2trt_additive_power <- function(RhoW, ModelChoice, n.individuals, n.clusters, n.periods,
                           delta_1, delta_2, RhoA, IAC, Sequence_Tx1,
                           Sequence_Tx2, typeIerror){


  if(ModelChoice %in% c("RCS Standardized", "RCS", "Cohort Standardized",
                        "Cohort", "NEC Standardized", "NEC")){ }
  else{
    print("Please specify ModelChoice as one of 6 options: 'RCS Standardized', 'RCS', 'Cohort Standardized', 'Cohort', 'NEC Standardized', or 'NEC'")
    stop()
  }


  ## variance parameters for repeated cross-sectional design
  if (ModelChoice == "RCS")
    { Diagonal <-  (1 - RhoW) / n.individuals
      OffDiagonal <- (RhoW) }
  ## variance parameters for standardized repeated cross-sectional design
  if (ModelChoice == "RCS Standardized")
	  { Diagonal <- RhoW + (1 - RhoW) / n.individuals
	    OffDiagonal <- (RhoW) }


  ## variance parameters for cohort design
  ## Updated 2/21 based on Phil's response
  if (ModelChoice == "Cohort")
    { Diagonal <-  RhoW + IAC * (1 - RhoW) / n.individuals
      # Diagonal <-  (1 - IAC - RhoW +  IAC * RhoW ) / n.individuals
      OffDiagonal <- (1 - RhoW - IAC + (IAC * RhoW)) / n.individuals
      # OffDiagonal <- RhoW + IAC * (1 - RhoW)/ n.individuals
      }
  ## variance parameters for standardized cohort design
	if (ModelChoice == "Cohort Standardized")
	  { Diagonal <-  RhoW + (1 - RhoW) / n.individuals
	    OffDiagonal <- RhoW + IAC * (1 - RhoW) / n.individuals }



  ## variance parameters for nested exchangeable model
  if (ModelChoice == "NEC")
    { Diagonal <-  RhoW + (1 - RhoW) / n.individuals  - RhoA
      OffDiagonal <- RhoA}
  else if (ModelChoice == "NEC Standardized")
    { Diagonal <-  RhoW + (1 - RhoW) / n.individuals
      OffDiagonal <- RhoA}

  ## creating data frame for design matrix
  StudyDesign <- data.frame(site = c(1:n.clusters))

  StudyDesign$TimeA <- NA
  StudyDesign$TimeB <- NA

  ## assignment of sequences
  StudyDesign$TimeA <- Sequence_Tx1
  StudyDesign$TimeB <- Sequence_Tx2

  ## repeat each cluster for each time period
  StudyDesign <- StudyDesign[rep(seq_len(nrow(StudyDesign)), each = n.periods), ]
  StudyDesign$Time <- rep(1:n.periods)

  ## Creating treatment indicators
  StudyDesign$Tx_A <- 0
  StudyDesign$Tx_B <- 0
  StudyDesign$Tx_A[StudyDesign$Time >= StudyDesign$TimeA] <- 1
  StudyDesign$Tx_B[StudyDesign$Time >= StudyDesign$TimeB] <- 1
  StudyDesign$Tx_AB <- StudyDesign$Tx_A * StudyDesign$Tx_B


### Calculate precision matrix ###

  ## quantities in paper
  f <- n.clusters / (Diagonal + n.periods * OffDiagonal)
  g <- n.clusters * OffDiagonal / (Diagonal*(Diagonal + n.periods * OffDiagonal))

  SumofXijbyperiodbycluster <- sum(StudyDesign$Tx_A)
  y_1 <- SumofXijbyperiodbycluster / (Diagonal + n.periods * OffDiagonal)
  h_1 <- SumofXijbyperiodbycluster / (Diagonal*(Diagonal + n.periods * OffDiagonal))
  l_1 <- SumofXijbyperiodbycluster / Diagonal

  SumofYijbyperiodbycluster <- sum(StudyDesign$Tx_B)
  y_2 <- SumofYijbyperiodbycluster / (Diagonal + n.periods * OffDiagonal)
  h_2 <- SumofYijbyperiodbycluster / (Diagonal * (Diagonal + n.periods * OffDiagonal))
  l_2 <- SumofYijbyperiodbycluster / Diagonal


  z_1 <-  0
  AggregateData_TreatmentASquared_ByCluster <- stats::aggregate(Tx_A ~ site,
          data = StudyDesign, FUN = function(x) c(mn = (sum(x))^2))
  z_1 <- sum(AggregateData_TreatmentASquared_ByCluster$Tx_A) *
         OffDiagonal / (Diagonal * (Diagonal + n.periods * OffDiagonal))

  z_2 <-  0
  AggregateData_TreatmentBSquared_ByCluster <- stats::aggregate(Tx_B ~ site,
          data=StudyDesign,  FUN = function(x) c(mn = (sum(x))^2))
  z_2 <- sum(AggregateData_TreatmentBSquared_ByCluster$Tx_B) *
         OffDiagonal / (Diagonal * (Diagonal + n.periods * OffDiagonal))

  z_3 <-  0
  AggregateData_TreatmentABSquared_ByCluster <- stats::aggregate(Tx_AB ~ site,
          data=StudyDesign,   FUN = function(x) c(mn = (sum(x))^2))
  z_3 <- sum(AggregateData_TreatmentABSquared_ByCluster$Tx_AB) *
         OffDiagonal / (Diagonal * (Diagonal + n.periods * OffDiagonal))

  q_part1 <-  sum(StudyDesign$Tx_AB) / Diagonal
  q_part2 <- 0
  AggregateData_TreatmentA_ByCluster <- stats::aggregate(Tx_A ~ site,
          data=StudyDesign, FUN = function(x) c(mn = (sum(x))))
  AggregateData_TreatmentB_ByCluster <- stats::aggregate(Tx_B~ site,
          data=StudyDesign, FUN = function(x) c(mn = (sum(x))))
  AggregateData_TreatmentAB_ByCluster <- stats::aggregate(Tx_AB~ site,
          data=StudyDesign, FUN = function(x) c(mn = (sum(x))))

  q_part2 <-  as.numeric(t(AggregateData_TreatmentA_ByCluster$Tx_A) %*%
                           AggregateData_TreatmentB_ByCluster$Tx_B *
                           OffDiagonal / ((Diagonal*(Diagonal + (n.periods * OffDiagonal)))))

  q_1 <- q_part1 - q_part2


  AggregateData_TreatmentASquared_ByPeriod <- stats::aggregate(Tx_A ~ Time,
                                              data=StudyDesign, FUN = function(x) c(mn = (sum(x))^2/Diagonal^2))
  w_1 <- sum(AggregateData_TreatmentASquared_ByPeriod$Tx_A)

  AggregateData_TreatmentBSquared_ByPeriod <- stats::aggregate(Tx_B ~ Time,
                                              data=StudyDesign, FUN = function(x) c(mn = (sum(x))^2))
  w_2 <- sum(AggregateData_TreatmentBSquared_ByPeriod$Tx_B) / Diagonal^2


  AggregateData_TxA <- stats::aggregate(Tx_A ~ Time, data=StudyDesign, FUN = function(x) c(mn = (sum(x))))

  AggregateData_TxB <- stats::aggregate(Tx_B ~ Time, data=StudyDesign, FUN = function(x) c(mn = (sum(x))))

  w_xy <- t(AggregateData_TxA$Tx_A) %*% AggregateData_TxB$Tx_B * (1/Diagonal^2)

  A_11 <- ((l_1 - z_1) * (f * n.periods * (f + g * n.periods)) -
             y_1^2 * (f + g * n.periods) - f * (n.periods * w_1 - l_1^2)) /
                                           (f * n.periods * (f + g * n.periods))

  A_22 <- ((l_2 - z_2) * (f * n.periods * (f + g * n.periods)) -
             y_2^2 * (f + g * n.periods) - f * (n.periods * w_2 - l_2^2)) /
                                           (f * n.periods * (f + g * n.periods))

  A_12 <- q_1 - y_1 * y_2 / (f * n.periods) - (1 / (f + g * n.periods)) * (w_xy - l_1 * l_2 / n.periods)

  PrecisionMatrix <- matrix(nrow = 2, ncol = 2)
  PrecisionMatrix[1, 1] <- A_11
  PrecisionMatrix[1, 2] <- PrecisionMatrix[2, 1] <- A_12
  PrecisionMatrix[2, 2] <- A_22


  ## get standard error of main effect of treatment 1
  SE_Tx1 <- sqrt((solve(PrecisionMatrix))[1, 1])
  TestStatistic_Tx1 <- delta_1 / SE_Tx1
  PowerTx1 <- stats::pnorm(TestStatistic_Tx1 - stats::qnorm(1 - typeIerror / 2)) +
              stats::pnorm(stats::qnorm(typeIerror / 2) - TestStatistic_Tx1)

  ## get standard error of main effect of treatment 2
  SE_Tx2 <- sqrt((solve(PrecisionMatrix))[2, 2])
  TestStatistic_Tx2 <- delta_2 / SE_Tx2
  PowerTx2 <- stats::pnorm(TestStatistic_Tx2 - stats::qnorm(1 - typeIerror / 2)) +
              stats::pnorm(stats::qnorm(typeIerror / 2) - TestStatistic_Tx2)

  return(list(PowerTx1 = PowerTx1, PowerTx2 = PowerTx2))
}
