
<!-- README.md is generated from README.Rmd. Please edit that file -->

# powertools

<!-- badges: start -->
<!-- badges: end -->

The goal of powertools is to provide Power and Sample Size Tools in R for researchers. This package accompanies the text "Power and Sample Size in R" by Catherine M. Crespi. 

This package loads all R packages that are used in the textbook and provides new functions for sample size and power calculations that did not previously exist.

## Installation

You can install the development version of PSSTools like so:

``` r
# install.packages("devtools")
devtools::install_github("zichen-liu/powertools")
```

## Example

This is a basic example which calculates the sample size needed to adequately power a study hoping to show that a new experimental therapy is promising based on the proportion of participants expected to respond to treatment. The treatment is considered promising if at least 20% of participants respond well, and the researchers believe the true response proportion is 30%. 

This example follows Example 5.1 in the textbook. 

``` r
library(powertools)
oneprop_ss(p0 = 0.2, pA = 0.3, alpha = 0.05, power = 0.8,
           method = "conditional", one.or.two.sided = "one")
#> [1] 109
```

`oneprop_ss` returns n, the required sample size for a one proportion test. 
