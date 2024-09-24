

# Load the required packages
library(devtools)
library(roxygen2)


devtools::install_github("zichen-liu/powertools")

library(powertools)


# Check the package
devtools::check()

# If there are errors, fix them as suggested in the output, then regenerate documentation
document()

# build the manual
build_manual()









