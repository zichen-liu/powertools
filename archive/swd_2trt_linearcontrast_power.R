#' Stepped Wedge Design Power Calculation for Linear Contrasts of Treatment Effects
#'
#' Calculate the power for detecting linear contrasts of two treatment effects in a two-treatment stepped wedge design (SWD) with additive treatment effects.
#' We only consider complete designs. This function is based on Power Analysis for Stepped Wedge Trials with Two Treatments published in Stat Med in 2022 by Phillip Sundin and Catherine Crespi. This manuscript uses the normal distribution and the Wald Test to derive power.
#'
#' @param RhoW The intraclass correlation coefficient, which refers to the correlation between outcomes of two different individuals in the same cluster at the same time.
#' @param ModelChoice A string indicating either a repeated cross sectional design ("RCS"), nested exchangeable ("NEC"), Cohort ("Cohort"), or any one of these three designs in their standardized form. Only one of "RCS Standardized", "RCS", "Cohort Standardized", "Cohort", "NEC Standardized", or "NEC" can be specified. The standardized design refers to dividing by the standard deviation of y (the outcome) in the covariance specification.
#' @param n.individuals The number of individuals in each condition at one time period.
#' @param n.clusters The number of clusters
#' @param n.periods The number of time periods.
#' @param delta_12 The standardized difference between the two treatment effects to be powered on.
#' @param RhoA The correlation between two observations in the same individual in the same cluster but different time periods. RCS models assume RhoW = RhoA, and Cohort models assume RhoA >= RhoW.
#' @param IAC The individual auto-correlation, which is the proportion of the individual-level variance that is time-invariant. IAC = 0 in a RCS.
#' @param Sequence_Tx1 A vector representing the time points (periods) at which each cluster transitions from control to treatment 1 (ie the first occurence).
#' @param Sequence_Tx2 A vector representing the time points (periods) at which each cluster transitions from control to treatment 2 (ie the first occurence).
#' @param typeIerror The significance level.
#'
#' @return Power for detecting the difference of treatment effects
#' @export
#'
#' @examples
#' # Late factorial SWD design with 12 clusters
#' swd_2trt_linearcontrast_power(RhoW = 0.2, ModelChoice="RCS", IAC = 0,
#'                       RhoA = 0.2, n.individuals = 15, delta_12 = 0.4,
#'                       n.clusters = 12, n.periods = 4,
#'                       Sequence_Tx1 = c(2,2,3,3,4,4,4,4,4,4,4,4),
#'                       Sequence_Tx2 = c(4,4,4,4,4,4,2,2,3,3,4,4),
#'                       typeIerror = 0.05)
#'
#' # Example where there are different potential designs to be considered #
#' #creating data frame for power as a function of rho_w
#' PowerTable <- data.frame(RhoW = c(rep(seq(0,0.4,by=0.01),times=3)))
#' PowerTable$DesignChoice <- c(rep('"Late" Factorial Design, 12 clusters',times=41),
#'                              rep('"Early" Factorial Design, 10 clusters',times=41),
#'                              rep("Concurrent Design, 12 clusters",times=41))
#'
#' PowerTable$n.clusters <- c(rep(12,times=41), rep(10,times=41), rep(12,times=41))
#'
#' PowerTable$Sequence1 <- c(rep(list(c(2,2,3,3,4,4,4,4,4,4,4,4)),times=41),
#'                           rep(list(c(2,2,2,3,4,4,3,4,4,4)),times=41),
#'                           rep(list(c(2,2,3,3,4,4,NA,NA,NA,NA,NA,NA)),times=41))
#'
#' PowerTable$Sequence2 <- c(rep(list(c(4,4,4,4,4,4,2,2,3,3,4,4)),times=41),
#'                           rep(list(c(4,4,4,3,4,4,3,2,2,2)),times=41),
#'                           rep(list(c(NA,NA,NA,NA,NA,NA,4,4,3,3,2,2)),times=41))
#'
#' PowerTable$Power <- mapply(swd_2trt_linearcontrast_power,
#'                         RhoW = PowerTable$RhoW,
#'                         ModelChoice="RCS",
#'                         IAC = 0, RhoA = PowerTable$RhoW,
#'                         n.individuals = 15, delta_12 = 0.4,
#'                         n.clusters = PowerTable$n.clusters,
#'                         n.periods = 4, typeIerror = 0.05,
#'                         Sequence_Tx1 = PowerTable$Sequence1,
#'                         Sequence_Tx2 = PowerTable$Sequence2)
#'
#' ##removing unnecessary columns of clusters
#' PowerTable <- PowerTable[ ,-which(names(PowerTable) %in%
#'             c("Sequence1","Sequence2","n.clusters"))]
#' PowerTable
swd_2trt_linearcontrast_power <- function(RhoW, ModelChoice, n.individuals, n.clusters,
                                          n.periods, delta_12, RhoA, IAC,
                                          Sequence_Tx1, Sequence_Tx2, typeIerror){

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

  StudyDesign <- data.frame(site = c(1:n.clusters))

  StudyDesign$TimeA <- Sequence_Tx1
  StudyDesign$TimeB <- Sequence_Tx2


##repeat each cluster for each time period
StudyDesign <- StudyDesign[rep(seq_len(nrow(StudyDesign)), each = n.periods), ]
StudyDesign$Time <- rep(1:n.periods)

#Creating treatment indicators
StudyDesign$Tx_A <- 0
StudyDesign$Tx_B <- 0
StudyDesign$Tx_A[StudyDesign$Time >= StudyDesign$TimeA] <- 1
StudyDesign$Tx_B[StudyDesign$Time >= StudyDesign$TimeB] <- 1
StudyDesign$Tx_AB <- StudyDesign$Tx_A * StudyDesign$Tx_B

#####################
##########Calculate precision matrix
#####################

##quantities in paper
f <- n.clusters /
  (Diagonal + n.periods * OffDiagonal)
g <- n.clusters*OffDiagonal /
  (Diagonal*(Diagonal + n.periods * OffDiagonal))

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
          data=StudyDesign,
      FUN = function(x) c(mn = (sum(x))^2))
z_1 <- sum(AggregateData_TreatmentASquared_ByCluster$Tx_A) *
  OffDiagonal / (Diagonal*(
  Diagonal + n.periods*OffDiagonal
))

z_2 <-  0
AggregateData_TreatmentBSquared_ByCluster <- stats::aggregate(Tx_B ~ site,
          data=StudyDesign,
      FUN = function(x) c(mn = (sum(x))^2))
z_2 <- sum(AggregateData_TreatmentBSquared_ByCluster$Tx_B) *
  OffDiagonal / (Diagonal*(
  Diagonal + n.periods*OffDiagonal
))

z_3 <-  0
AggregateData_TreatmentABSquared_ByCluster <- stats::aggregate(Tx_AB ~ site,
          data=StudyDesign,
      FUN = function(x) c(mn = (sum(x))^2))
z_3 <- sum(AggregateData_TreatmentABSquared_ByCluster$Tx_AB) *
  OffDiagonal / (Diagonal*(
  Diagonal + n.periods*OffDiagonal
))

q_part1 <- 0
q_part1 <-  sum(StudyDesign$Tx_AB) /Diagonal

q_part2 <- 0
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

w_1 <- 0
AggregateData_TreatmentASquared_ByPeriod <- stats::aggregate(Tx_A ~ Time,
  data=StudyDesign,
      FUN = function(x) c(mn = (sum(x))^2/Diagonal^2))
w_1 <- sum(AggregateData_TreatmentASquared_ByPeriod$Tx_A)

AggregateData_TreatmentBSquared_ByPeriod <- stats::aggregate(Tx_B ~ Time,
  data=StudyDesign,
      FUN = function(x) c(mn = (sum(x))^2))
w_2 <- sum(AggregateData_TreatmentBSquared_ByPeriod$Tx_B) / Diagonal^2

AggregateData_TxA <- stats::aggregate(Tx_A ~ Time,
          data=StudyDesign,
      FUN = function(x) c(mn = (sum(x))))

AggregateData_TxB <- stats::aggregate(Tx_B ~ Time,
          data=StudyDesign,
      FUN = function(x) c(mn = (sum(x))))


w_xy <- t(AggregateData_TxA$Tx_A) %*% AggregateData_TxB$Tx_B *
  (1/Diagonal^2)


A_11 <-  ((l_1 - z_1)*(f*n.periods*(f+g*n.periods)) -
    y_1^2*(f+g*n.periods) - f*(n.periods*w_1 - l_1^2)) /
  (f*n.periods*(f+g*n.periods))

A_22 <- ((l_2 - z_2)*(f*n.periods*(f+g*n.periods)) -
    y_2^2*(f+g*n.periods) - f*(n.periods*w_2 - l_2^2)) /
  (f*n.periods*(f+g*n.periods))


A_12 <- q_1 - y_1*y_2/(f*n.periods) - (1/(f+g*n.periods))*(
  w_xy - l_1*l_2/n.periods)

PrecisionMatrix <- matrix(nrow=2,ncol=2)
PrecisionMatrix[1,1] <- A_11
PrecisionMatrix[1,2] <- PrecisionMatrix[2,1] <- A_12
PrecisionMatrix[2,2] <- A_22


InversePrecision2x2 <- (solve(PrecisionMatrix))
SE_TxAMinusTxB <- sqrt(InversePrecision2x2[1,1] +
  InversePrecision2x2[2,2] - 2 * InversePrecision2x2[2,1])

TestStat_TxAminusTxB <- delta_12 / SE_TxAMinusTxB

PowerTxAMinusTxB <- stats::pnorm(TestStat_TxAminusTxB - stats::qnorm(1-typeIerror/2))+
 stats::pnorm(stats::qnorm(typeIerror/2) - TestStat_TxAminusTxB)

return(PowerTxAMinusTxB)
}
