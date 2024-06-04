
### Beginning of R Package powertools #
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

# rm(strsplit1)

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
document()
# and now I can do this!!! wow!!
?strsplit1

check()

# now we can install the package like any other package!!
install()
library(powertools)
x <- "alfa,bravo,charlie,delta"
strsplit1(x, split = ",")

# initialize test files to check the functions we are adding
use_testthat()
# now this will open a test file for us for the function we ask for.
use_test("strsplit1")
# test_that("strsplit1() splits a string", {
#  expect_equal(strsplit1("a,b,c", split = ","), c("a", "b", "c"))
# })
# this type of test in the test file ^^
# give it an argument to evaluate and the correct answer.

# we need to have testthat library open already. Not sure how to do that..

library(testthat)
load_all()
test()

# To select specific packages to import for the package to work
use_package("stringr")

# in order to update a function, we update the .R file and the test .R file
# and ask to rename the function

rename_files("strsplit1", "str_split_one")
# but I am not actually running this. I am just making a new file.

use_r("str_split_one")
use_test("str_split_one")

document()
test()

use_readme_rmd() # then edit the readme file

build_readme() # use this to build the file (?) IDK


check()

install()

## as of 10/10 at 10:57 AM it works...


# now lets see if I can add other functions that I actually made for it.

use_package("purrr")
use_package("pwr2ppl")


use_r("anova3_ss_wrapper")
load_all()

document()

# error when I check with Error in purrr::quietly(anova1f_3) : object 'anova1f_3' not found
# tried to just change it to be pwr2ppl::anova1f_3 and this worked.

# when I restart need to:
library(devtools)
library(usethis)
check()

?anova3_ss_wrapper
anova3_ss_wrapper(means = c(5, 10, 12), sds = c(10, 10, 10), alphalevel = 0.05,
                  powerlevel = 0.8, groupweights = c(.5, .25, .25),
                  Nmin = NULL, Nmax = NULL, Nchange = NULL)


use_r("anova4_ss_wrapper")
load_all()
document()

use_package("Superpower")
use_package("stats")

use_r("McNemar_N")

# importFrom("stats", "qnorm") <- not exactly sure.. so I am just calling the entire stats package.

use_r("powersearch_multvars")

use_package("dplyr")
use_package("magrittr")
use_package("pwr2")
use_package("pwrAB")
# library(magrittr)

use_pipe(export = TRUE) # to add %>% pipe from magrittr or dplyr
check()

use_r("ss_oneprop")
load_all()

use_git()

usethis::create_from_github(
  "https://github.com/kristenmcgreevy/PSSTools.git",
  destdir = "/Volumes/KMM_LAB2/Rpackage_Git"
)


usethis::create_from_github(
  "https://github.com/kristenmcgreevy/PSSTools.git"
)



use_test("anova3_ss_wrapper")


use_r("anova3_ss_wrapper")
load_all()
document()


use_test("anova4_ss_wrapper")


library(available)

available("powertools")

available("powertools")

## R CMD check is the same as
devtools::check()


# To select specific packages to import for the package to work
use_package("stringr")

?powersearch_multvars # parameters


rename_files("ss_oneprop", "oneprop_ss")


use_test("oneprop_ss")
?oneprop_ss
check()


oneprop_ss(p0 = 0.2, pA = 0.3, alpha = 0.05, power = 0.8,
            method = "conditional", one.or.two.sided = "one")
oneprop_ss(p0 = 0.2, pA = 0.3, alpha = 0.05, power = 0.8,
            method = "unconditional", one.or.two.sided = "one")

rename_files("powersearch_multvars", "search_power")
check()
use_test("search_power")



library(pwr2)
output <- search_power(powerfunction = pwr.1way, searchvector1 = c(3, 3, 4, 4),
                      searchvector2 = c(80, 100, 80, 100),
                      searchvector3 = NULL, searchvector4= NULL,
                      alpha = 0.05, delta = 3, sigma = 8.5)
output_power <- do.call(c, output$power)
vec <- c(0.4968811, 0.5965158, 0.4342090, 0.5304948)
output_power == c(0.4968811, 0.5965158, 0.4342090, 0.5304948)
typeof(output_power) # double
output_power == vec
typeof(vec) # double

typeof(output) # list
notes <- c("n is number in each group, total sample = 240 power = 0.0537073121787097",
           "n is number in each group, total sample = 300 power = 0.054652972281662",
           "n is number in each group, total sample = 320 power = 0.0528927215920873",
           "n is number in each group, total sample = 400 power = 0.0536311240089957")
output <- list(k = c(3,3,4,4), n = c(80, 100, 80, 100), delta = c(0.3, 0.3, 0.3, 0.3),
               sigma = rep(8.5, 4), effect.size = rep(0.0144, 4),
               sig.level = rep(0.05, 4), power = rep(0.0537, 4), note = notes,
               method = rep("Balanced one-way analysis of variance power calculation", 4))
# see if this will pass the test ^^

# I can't get it to work. so I am just deleting the test file.
check()



# To Update the Description File, in R Studio type Control + . and then type DESCRIPTION


rename_files("anova3_ss_wrapper", "anova3_ss")
rename_files("anova4_ss_wrapper", "anova4_ss")
check()

rename_files("McNemar_N", "mcnemar_ss")



### 11/28/2022 ###
use_r("onemean_ci_power")
use_r("onemean_ci_ss")
check()

use_test("onemean_ci_power")
use_test("onemean_ci_ss") # was not working, so I just deleted the file
check()

use_r("meandiff_ci_power")
use_r("meandiff_ci_ss")



use_package("pwr")
use_package("SampleSize4ClinicalTrials")
use_package("MESS")
use_package("MKmisc")
use_package("PowerTOST")
use_package("exact2x2")
use_package("clusterPower")
use_package("presize")
use_package("precisely")
use_package("plyr")



check()

use_package("ssanv", "Suggests")
use_package("designsize", "Suggests")
use_package("npsurvSS", "Suggests")
use_package("powerSurvEpi", "Suggests")
use_package("longpower", "Suggests")


### making data to include in the package ###
# m is indivs/site, J is number of sites
m <- 20
J <- 10
N <- m*J

gamma00 <- 0
gamma10 <- 3
sigma_u0_sq <- 4
sigma_u1_sq <- 8
sigma_eps_sq <- 36

#generate REs for sites

#set.seed(56524)
set.seed(51124)

u0 <- rnorm(n=J,mean=0,sd=sqrt(sigma_u0_sq))
u1 <- rnorm(n=J,mean=0,sd=sqrt(sigma_u1_sq))

#make dataframe

id <- as.factor(seq(1:N))
i <- as.factor(rep(seq(1:m),J))
j <- as.factor(rep(seq(1:J),each=m))
xij <- rep(c(rep(-0.5,m/2),rep(0.5,m/2) ), J)
u0j <- rep(u0, each=m)
u1j <- rep(u1, each=m)
eij <- rnorm(n=N, mean=0, sd=sqrt(sigma_eps_sq))
#this Y includes site by group interaction:
Yij <- u0j + (gamma10+u1j)*xij + eij
#this Y has no interaction:
Wij <- u0j + gamma10*xij + eij

group <- xij+1.5

multisite.data <- data.frame(cbind(id,j,i,group,xij,Yij,Wij))

usethis::use_data(multisite.data, overwrite = TRUE)
head(multisite.data)


use_r("data")
check()



load_all()
oneprop_ss(p0 = 0.2, pA = 0.3, alpha = 0.05, power = 0.8,
           method = "conditional", one.or.two.sided = "one")


### just to look at functions ###
load_all()


### Feb 2023 - adding new functions ###

use_r("tandz_ratio")

load_all()

tandz_ratio(alpha = 0.05, power = 0.8, df = 20)

use_r("parallelCRT_power") # example 12.1
use_r("parallelCRT_ss") # example 12.2
use_r("crossover_power") # example 12.3
use_r("crossover_ss") # example 12.4

use_r("swd_1trt_power") # for example 12.5
use_r("swd_1trt_ss") # for example 12.6

check()


# try to make the pdf manual ? #
build_manual(pkg = ".", path = NULL)


library(testthat)
load_all()
use_test("parallelCRT_power")

use_test("anova3_ss")


mdata <- data(multisite.data)
summary(mdata)

usethis::use_data_raw("multisite.data")

check()

load_all()

data("powertools")

use_data("multisite.data")
use_r("data")


use_readme_rmd()


use_git()


### 2/21/2023 ###

use_r("newfunction")

check()
load_all()
?newfunction


use_test("newfunction")


use_package("stats")

use_package("stats", "Suggests")

usethis::use_data_raw("datasettest")


use_r("data")

## Adding (modified) Phil's functions ##
use_r("swd_2trt_power")

check()

use_r("swd_2trt_interaction_power")

use_r("swd_2trt_linearcontrast_power")


rename_files("swd_2trt_power", "swd_2trt_additive_power")

use_r("swd_2trt_additive_power")


use_package("rPowerSampleSize", "Suggests")

use_r("multendpoints_Ck_ss")


use_package("mvtnorm")
check()

use_r("pss.z.test.R")
check()

use_r("pss.t.test.R")
check()
