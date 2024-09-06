library(devtools)
library(roxygen2)
library(Rdpack)
install_github("zichen-liu/powertools")
library(powertools)

devtools::check()

devtools::document()
devtools::build_manual()

