# Hello, world!
#
# This is an example function named 'hello'
# which prints 'Hello, world!'.
#
# You can learn more about package authoring with RStudio at:
#
#   http://r-pkgs.had.co.nz/
#
# Some useful keyboard shortcuts for package authoring:
#
#   Install Package:           'Cmd + Shift + B'
#   Check Package:             'Cmd + Shift + E'
#   Test Package:              'Cmd + Shift + T'

hello <- function() {
  print("Hello, world!")
}


### Beginning of R Package PSStools #
### https://r-pkgs.org/Whole-game.html
library(devtools)

packageVersion("devtools")
# create_package("~/Box/Mine/Crespi_GSR/regexcite")
# create_package("/Users/kristenmae/Box/regexcite")
# neither of those were working properly, so I just made a new one in my Volumes/KMM_LAB2


use_git() # once I restart this project (by closing and reopening the .Rproj file in /regexcite),
# you can see the history of commands (from this and previous sessions) in the history panel! thats pretty cool


use_r("strsplit1")

# strsplit1 <- function(x, split) {
#  strsplit(x, split = split)[[1]]
# }

rm(strsplit1)

# call this new function by loading it this way:
load_all()

(x <- "alfa,bravo,charlie,delta")
# so even though I don't see the function in the environment pane, it is still loaded!
strsplit1(x, split = ",")

# make sure this is false and we are using the function loaded
exists("strsplit1", where = globalenv(), inherits = FALSE)

# use this to check the R package
check()

use_mit_license()

# trigger the conversion of roxygen code in the functions by using
document("~/Box/Mine/Crespi_GSR/")

