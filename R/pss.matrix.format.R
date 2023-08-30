pss.matrix.format <- function (matrix = NULL, nspaces = 18) {
  mrows <- c()
  for (i in 1:nrow(matrix))
    mrows <- c(mrows, paste(matrix[i,], collapse = ', '))
  div <- paste0("\n", strrep(" ", nspaces))
  paste(mrows, collapse = div)
}


