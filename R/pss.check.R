pss.check <- function(param = NULL,
                      type = c("req", "num", "unit", "pos", "int", "bool", "vals", "min"),
                      valslist = NULL, min = NULL) {

  name <- deparse(substitute(param))

  # required variable
  if (type == "req") {
    if (is.null(param))
      stop(paste(name, "must be specified"))
  }

  # the following checks are done assuming the variable isn't null
  else {
    if (!is.null(param)){

      # numeric variables
      if (type == "num") {
        if (!is.numeric(param))
          stop(paste(name, "should be a numeric value"))
      }

      # variables between 0 and 1
      else if (type == "unit") {
        if (!is.numeric(param))
          stop(paste(name, "should be a numeric value"))
        if (param >= 1 || param <= 0)
          stop(paste(name, "should be between 0 and 1"))
      }

      # positive values only
      else if (type == "pos") {
        if (!is.numeric(param))
          stop(paste(name, "should be a numeric value"))
        if (param <= 0)
          stop(paste(name, "should be a positive numeric value"))
      }

      # integers only
      else if (type == "int") {
        if (!is.numeric(param))
          stop(paste(name, "should be a numeric value"))
        if (param %% 1 != 0)
          stop(paste(name, "should be a positive integer"))
        if (param <= 0)
          stop(paste(name, "should be positive"))
      }

      # true / false only
      else if (type == "bool") {
        if (!isTRUE(param) & !isFALSE(param))
          stop(paste(name, "should be either TRUE or FALSE"))
      }

      # only certain values allowed
      else if (type == "vals") {
        if (!(param %in% valslist))
          stop(paste(name, "should be one of:", paste(valslist, collapse = ", ")))
      }

      # value has a minimum
      else if (type == "min") {
        if (!is.numeric(param))
          stop(paste(name, "should be a numeric value"))
        if (param < min)
          stop(paste(name, "should be greater than", min))
      }
    }
  }
}
