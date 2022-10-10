
<!-- README.md is generated from README.Rmd. Please edit that file -->

# PSSTools

<!-- badges: start -->
<!-- badges: end -->

The goal of PSSTools is to provide Power and Sample Size Tools in R for
researchers.

## Installation

You can install the development version of PSSTools like so:

``` r
# install.packages("devtools")
devtools::install_github("kristenmcgreevy/PSSTools")
```

## Example

This is a basic example which shows you how to solve a common problem:

``` r
library(PSSTools)
(x <- "alfa,bravo,charlie,delta")
#> [1] "alfa,bravo,charlie,delta"
str_split_one(x, pattern = ",")
#> [1] "alfa"    "bravo"   "charlie" "delta"
```

Use `str_split_one()` when the input is known to be a single string. For
safety, it will error if its input has length greater than one.

`str_split_one()` is built on `stringr::str_split()`, so you can use its
`n` argument and stringr’s general interface for describing the
`pattern` to be matched.

``` r
str_split_one(x, pattern = ",", n = 2)
#> [1] "alfa"                "bravo,charlie,delta"

y <- "192.168.0.1"
str_split_one(y, pattern = stringr::fixed("."))
#> [1] "192" "168" "0"   "1"
```

What is special about using `README.Rmd` instead of just `README.md`?
You can include R chunks like so:

You’ll still need to render `README.Rmd` regularly, to keep `README.md`
up-to-date. `devtools::build_readme()` is handy for this. You could also
use GitHub Actions to re-render `README.Rmd` every time you push. An
example workflow can be found here:
<https://github.com/r-lib/actions/tree/v1/examples>.

In that case, don’t forget to commit and push the resulting figure
files, so they display on GitHub and CRAN.
