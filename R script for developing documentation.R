library(devtools)
library(roxygen2)
library(Rdpack)
devtools::check()

devtools::document()
devtools::build_manual()

