#' Stepped Wedge Design Power Calculation for Interaction Term
#'
#' Calculate power for detecting interaction effects in a two-treatment stepped wedge design.
#' We only consider complete designs. This function is based on Power Analysis for Stepped Wedge Trials with Two Treatments published in Stat Med in 2022 by Phillip Sundin and Catherine Crespi. This manuscript uses the normal distribution and the Wald Test to derive power.
#'
#'
#' @param RhoW The intraclass correlation coefficient, which refers to the correlation between outcomes of two different individuals in the same cluster at the same time.
#' @param ModelChoice A string indicating either a repeated cross sectional design ("RCS"), nested exchangeable ("NEC"), Cohort ("Cohort"), or any one of these three designs in their standardized form. Only one of "RCS Standardized", "RCS", "Cohort Standardized", "Cohort", "NEC Standardized", or "NEC" can be specified. The standardized design refers to dividing by the standard deviation of y (the outcome) in the covariance specification.
#' @param n.individuals The number of individuals in each condition at one time period.
#' @param n.clusters The number of clusters
#' @param n.periods The number of time periods.
#' @param delta_1 The standardized effect size for treatment 1
#' @param delta_2 The standardized effect size for treatment 2
#' @param delta_3 The standardized effect size for the interaction effect.
#' @param RhoA The correlation between two observations in the same individual in the same cluster but different time periods. RCS models assume RhoW = RhoA, and Cohort models assume RhoA >= RhoW.
#' @param IAC The individual auto-correlation, which is the proportion of the individual-level variance that is time-invariant. IAC = 0 in a RCS.
#' @param typeIerror The significance level.
#' @param Sequence_Tx1 A vector representing the time points (periods) at which each cluster transitions from control to treatment 1 (ie the first occurence).
#' @param Sequence_Tx2 A vector representing the time points (periods) at which each cluster transitions from control to treatment 2 (ie the first occurence).
#'
#' @return A matrix with the standard error, statistic, and power for detecting effects for treatment 1, treatment 2, and interaction effects.
#' @export
#'
#' @examples
#' # 8 Cluster repeated cross sectional SWD with main and interaction effects assumed
#' # to be the same with 0.6 effect size.
#' swd_2trt_interaction_power(RhoW = 0.2, ModelChoice = "RCS", RhoA = 0.2, IAC = 0,
#'                 n.individuals = 15, n.clusters = 8, n.periods = 5, delta_1 = 0.6,
#'                 delta_2 = 0.6, delta_3 = 0.6, typeIerror = 0.05/3,
#'                 Sequence_Tx1 = c(2,3,2,3,4,5,NA,NA), Sequence_Tx2 = c(NA,NA,2,3,4,5,3,2))
swd_2trt_interaction_power <- function(RhoW, ModelChoice, n.individuals, n.clusters, n.periods,
                                       delta_1, delta_2, delta_3, RhoA, IAC,
                                       typeIerror, Sequence_Tx1, Sequence_Tx2){

   if(ModelChoice %in% c("RCS Standardized", "RCS", "Cohort Standardized",
                        "Cohort", "NEC Standardized", "NEC")){ }
    else{
      print("Please specify ModelChoice as one of 6 options: 'RCS Standardized', 'RCS', 'Cohort Standardized', 'Cohort', 'NEC Standardized', or 'NEC'")
      stop()
    }


  Sequence_Tx1 <- unlist(Sequence_Tx1)
  Sequence_Tx2 <- unlist(Sequence_Tx2)

  ## variance parameters for repeated cross-sectional design
  if (ModelChoice == "RCS")
  { Diagonal <-  (1-RhoW) / n.individuals
    OffDiagonal <- (RhoW) }
  ## variance parameters for standardized repeated cross-sectional design
  if (ModelChoice == "RCS Standardized")
	  { Diagonal <- RhoW + (1 - RhoW) / n.individuals
	    OffDiagonal <- (RhoW) }

  ## variance parameters for cohort design
  else if (ModelChoice == "Cohort")
  { Diagonal <- RhoW + IAC * (1 - RhoW) / n.individuals
    OffDiagonal <- (1 - RhoW - IAC + (IAC * RhoW)) / n.individuals }
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

  StudyDesign <- data.frame(site = c(1:n.clusters))

  StudyDesign$TimeA <- Sequence_Tx1
  StudyDesign$TimeB <- Sequence_Tx2

  ## repeat each cluster for each time period
  StudyDesign <- StudyDesign[rep(seq_len(nrow(StudyDesign)), each = n.periods), ]
  StudyDesign$Time <- rep(1:n.periods)

  StudyDesign$Tx_A <- 0
  StudyDesign$Tx_B <- 0
  StudyDesign$Tx_A[StudyDesign$Time >= StudyDesign$TimeA] <- 1
  StudyDesign$Tx_B[StudyDesign$Time >= StudyDesign$TimeB] <- 1
  StudyDesign$Tx_AB <- StudyDesign$Tx_A * StudyDesign$Tx_B

  rownames(StudyDesign) <- 1:nrow(StudyDesign)

  #####################
  ##########Calculate precision matrix
  #####################

  ## quantities in paper
  f <- n.clusters /   (Diagonal + n.periods * OffDiagonal)
  g <- n.clusters*OffDiagonal / (Diagonal*(Diagonal + n.periods * OffDiagonal))

  SumofXijbyperiodbycluster <- sum(StudyDesign$Tx_A)
  y_1 <- SumofXijbyperiodbycluster / (Diagonal + n.periods * OffDiagonal)
  h_1 <- SumofXijbyperiodbycluster / (Diagonal*(Diagonal + n.periods * OffDiagonal))
  l_1 <-  SumofXijbyperiodbycluster / Diagonal

  SumofYijbyperiodbycluster <- sum(StudyDesign$Tx_B)
  y_2 <- SumofYijbyperiodbycluster / (Diagonal + n.periods * OffDiagonal)
  h_2 <- SumofYijbyperiodbycluster / (Diagonal*(Diagonal + n.periods * OffDiagonal))
  l_2 <- SumofYijbyperiodbycluster / Diagonal

  SumofXYijbyperiodbycluster <- sum(StudyDesign$Tx_AB)
  y_3 <- SumofXYijbyperiodbycluster / (Diagonal + n.periods * OffDiagonal)
  h_3 <- SumofXijbyperiodbycluster / (Diagonal*(Diagonal + n.periods * OffDiagonal))
  l_3 <-  SumofXYijbyperiodbycluster / Diagonal

  z_1 <-  0
  AggregateData_TreatmentASquared_ByCluster <- stats::aggregate(Tx_A ~ site,
          data=StudyDesign, FUN = function(x) c(mn = (sum(x))^2))
  z_1 <- sum(AggregateData_TreatmentASquared_ByCluster$Tx_A) *
  OffDiagonal / (Diagonal*(Diagonal + n.periods*OffDiagonal))

  z_2 <-  0
  AggregateData_TreatmentBSquared_ByCluster <- stats::aggregate(Tx_B ~ site,
          data=StudyDesign, FUN = function(x) c(mn = (sum(x))^2))
  z_2 <- sum(AggregateData_TreatmentBSquared_ByCluster$Tx_B) *
  OffDiagonal / (Diagonal*(Diagonal + n.periods*OffDiagonal))

  z_3 <-  0
  AggregateData_TreatmentABSquared_ByCluster <- stats::aggregate(Tx_AB ~ site,
          data=StudyDesign, FUN = function(x) c(mn = (sum(x))^2))
  z_3 <- sum(AggregateData_TreatmentABSquared_ByCluster$Tx_AB) *
  OffDiagonal / (Diagonal*(Diagonal + n.periods*OffDiagonal))

  q_part1 <-  sum(StudyDesign$Tx_AB) / Diagonal

  AggregateData_TreatmentA_ByCluster <- stats::aggregate(Tx_A ~ site,
          data=StudyDesign,
      FUN = function(x) c(mn = (sum(x))))
  AggregateData_TreatmentB_ByCluster <- stats::aggregate(Tx_B~ site,
          data=StudyDesign,
      FUN = function(x) c(mn = (sum(x))))
  AggregateData_TreatmentAB_ByCluster <- stats::aggregate(Tx_AB~ site,
          data=StudyDesign,
      FUN = function(x) c(mn = (sum(x))))

  q_part2 <-  as.numeric(t(AggregateData_TreatmentA_ByCluster$Tx_A) %*% AggregateData_TreatmentB_ByCluster$Tx_B *
  OffDiagonal / ((Diagonal*(Diagonal + (n.periods * OffDiagonal)))))

  q_1 <- q_part1 - q_part2

  q2_part2 <- as.numeric(t(AggregateData_TreatmentA_ByCluster$Tx_A) %*% AggregateData_TreatmentAB_ByCluster$Tx_AB *
                     OffDiagonal/(Diagonal*(Diagonal+n.periods*OffDiagonal)))

  q3_part2 <- as.numeric(t(AggregateData_TreatmentB_ByCluster$Tx_B) %*% AggregateData_TreatmentAB_ByCluster$Tx_AB *
                     OffDiagonal/(Diagonal*(Diagonal+n.periods*OffDiagonal)))

  q2 <- l_3 - q2_part2
  q3 <- l_3 - q3_part2

  AggregateData_TreatmentASquared_ByPeriod <- stats::aggregate(Tx_A ~ Time,
  data=StudyDesign,
      FUN = function(x) c(mn = (sum(x))^2/Diagonal^2))
  w_1 <- sum(AggregateData_TreatmentASquared_ByPeriod$Tx_A)

  AggregateData_TreatmentBSquared_ByPeriod <- stats::aggregate(Tx_B ~ Time,
  data=StudyDesign,
      FUN = function(x) c(mn = (sum(x))^2))
  w_2 <- sum(AggregateData_TreatmentBSquared_ByPeriod$Tx_B) / Diagonal^2

  AggregateData_TreatmentABSquared_ByPeriod <- stats::aggregate(Tx_AB ~ Time,
  data=StudyDesign,
      FUN = function(x) c(mn = (sum(x))^2))
  w_3 <- sum(AggregateData_TreatmentABSquared_ByPeriod$Tx_AB) / Diagonal^2

  AggregateData_TxA <- stats::aggregate(Tx_A ~ Time,
  data=StudyDesign,  FUN = function(x) c(mn = (sum(x))))

  AggregateData_TxB <- stats::aggregate(Tx_B ~ Time,
      data=StudyDesign, FUN = function(x) c(mn = (sum(x))))

  AggregateData_TxAB <- stats::aggregate(Tx_AB ~ Time,
  data=StudyDesign, FUN = function(x) c(mn = (sum(x))))

  w_xy <- t(AggregateData_TxA$Tx_A) %*% AggregateData_TxB$Tx_B * (1/Diagonal^2)

  w_x_xy <- t(AggregateData_TxA$Tx_A) %*% AggregateData_TxAB$Tx_AB *  (1/Diagonal^2)

  w_y_xy <- t(AggregateData_TxB$Tx_B) %*% AggregateData_TxAB$Tx_AB *  (1/Diagonal^2)

  A_11 <-  ((l_1 - z_1)*(f*n.periods*(f+g*n.periods)) -
    y_1^2*(f+g*n.periods) - f*(n.periods*w_1 - l_1^2)) / (f*n.periods*(f+g*n.periods))

  A_22 <- ((l_2 - z_2)*(f*n.periods*(f+g*n.periods)) -
    y_2^2*(f+g*n.periods) - f*(n.periods*w_2 - l_2^2)) /  (f*n.periods*(f+g*n.periods))

  A_33 <- ((l_3 - z_3)*(f*n.periods*(f+g*n.periods)) -
    y_3^2*(f+g*n.periods) - f*(n.periods*w_3 - l_3^2)) /  (f*n.periods*(f+g*n.periods))


  A_12 <- q_1 - y_1*y_2/(f*n.periods) - (1/(f+g*n.periods))*(
  w_xy - l_1*l_2/n.periods)
  A_13 <- q2 - y_1*y_3/(f*n.periods) -
  (1/(f+g*n.periods))*(w_x_xy - l_1*l_3/n.periods)
  A_23 <- q3 - y_2*y_3/(f*n.periods) -
  (1/(f+g*n.periods))*(w_y_xy - l_2*l_3/n.periods)


  PrecisionMatrix <- matrix(nrow=3,ncol=3)
  PrecisionMatrix[1,1] <- A_11
  PrecisionMatrix[1,2] <-  PrecisionMatrix[2,1] <- A_12
  PrecisionMatrix[2,2] <- A_22
  PrecisionMatrix[1,3] <- PrecisionMatrix[3,1] <- A_13
  PrecisionMatrix[3,3] <- A_33
  PrecisionMatrix[2,3] <- PrecisionMatrix[3,2] <- A_23

  InversePrecision <- solve(PrecisionMatrix)

  ResultsTable <- matrix(0,nrow=3,ncol=3)
  colnames(ResultsTable) <- c("Main Effect of Treatment 1",
                    "Main Effect of Treatment 2",
                    "Interaction Effect")
  rownames(ResultsTable) <- c("Standard Error (Sqrt of Variances)",
                    "Statistic (Estimated Value / SE)", "Power")

  ## Reading in the Standard Errors
  ResultsTable[1,1] <- round(sqrt(InversePrecision[1,1]),4)
  ResultsTable[1,2] <- round(sqrt(InversePrecision[2,2]),4)
  ResultsTable[1,3] <- round(sqrt(InversePrecision[3,3]),4)

  ResultsTable[2,1] <- delta_1 / ResultsTable[1,1]
  ResultsTable[2,2] <- delta_2 / ResultsTable[1,2]
  ResultsTable[2,3] <- delta_3 / ResultsTable[1,3]

  ############calculating power for treatment A
  ResultsTable[3,1] <- stats::pnorm(ResultsTable[2,1]  - stats::qnorm(1-typeIerror/2),)+
  stats::pnorm(stats::qnorm(typeIerror/2) - ResultsTable[2,1] )

  ############calculating power for treatment B
  ResultsTable[3,2] <- stats::pnorm(ResultsTable[2,2] - stats::qnorm(1-typeIerror/2))+
  stats::pnorm(stats::qnorm(typeIerror/2) - ResultsTable[2,2])

  ############calculating power for treatment AB
  ResultsTable[3,3] <- stats::pnorm(ResultsTable[2,3] - stats::qnorm(1-typeIerror/2),)+
  stats::pnorm(stats::qnorm(typeIerror/2) - ResultsTable[2,3])

  return(ResultsTable)
}
