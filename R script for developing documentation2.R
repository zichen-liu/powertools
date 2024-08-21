

# Load the required packages
library(devtools)
library(roxygen2)

# Check the package
devtools::check()

# If there are errors, fix them as suggested in the output, then regenerate documentation
document()

# build the manual
build_manual()




relrisk(n1 = 40, n.ratio = 1, p1 = 0.5, p2 = 0.3, power = NULL)

epi.sscohortc(N = NA, irexp1 = 0.3, irexp0 = 0.5, pexp = NA, n = 80,
              power = NA, r = 1, design = 1, sided.test = 2,
               conf.level = 0.95)

