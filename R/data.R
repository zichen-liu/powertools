# #' Multisite Dataframe
# #'
# #' A simulated dataset for multisite trials with and without interactions.
# #'
# #' @format ## `multisite.data`
# #'
# #' A data frame with 200 rows and 9 columns:
# #' \describe{
# #'   \item{id}{Individual observation number which goes from id = 1 to 200 (N)}
# #'   \item{j}{Site number which runs from j = 1 to 10 (J)}
# #'   \item{i}{Within site observation number. This runs from i = 1 to 20 (m)}
# #'   \item{xij}{Group level indicator variable coded so that within each site xij sums to 0.}
# #'   \item{Yij}{Outcome variable with a site by group interaction.}
# #'   \item{Wij}{Outcome variable without a site by group interaction.}
# #' }
# #'
# "multisite.data"

#' Delta Sign Table
#'
#' A table documenting the sign of delta for tests comparing two parameters.
#'
#' @format ## `delta.sign`
#' A data frame with 3 rows and 3 columns:
#' \describe{
#'   \item{columns}{Whether the test is for noninferiority of superirority by a margin}
#'   \item{rows}{Whether a higher or lower parameter value is better}
#'   \item{xij}{The sign of delta given the row & column}
#' }
"delta.sign"
