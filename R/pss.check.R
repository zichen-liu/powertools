pss.check <- function(param = NULL, name = NULL,
                      type = c("unit", "pos", "bool")) {

  if (type == "unit") {
    if (!is.numeric(param))
      stop(paste0(name, " should be a numeric value"))
    if (param > 1 || param < 0)
      stop(paste0(name, " should be between 0 and 1"))
  }

  else if (type == "pos") {
    if (!is.numeric(param))
      stop(paste0(name, " should be a numeric value"))
    if (param <= 0)
      stop(paste0(name, " should be positive"))
  }

  else if (type == "bool") {
    if (!isTRUE(param) & !isFALSE(param))
      stop(paste0(name, " should be either TRUE or FALSE"))
  }
}
