check.many <- function(paramlist = NULL,
                      type = c("oneof")) {

  names <- lapply(substitute(paramlist), deparse)
  names <- paste(names[-1], collapse = ", ")

  # variables where only one must be missing
  if (type == "oneof") {
    if (sum(sapply(paramlist, is.null)) != 1)
      stop(paste("exactly one of the following arguments must be NULL:", names))
  }
}
