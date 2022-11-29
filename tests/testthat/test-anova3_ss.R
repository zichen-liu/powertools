test_that("anova3_ss() outputs proper sample size and power", {
  expect_equal(anova3_ss(means = c(4,5,6), sds = c(1,2,5), alphalevel = 0.05,
                        powerlevel = 0.8, groupratios = c(1,1,3)),
               data.frame(n1 = 49, n2 = 49, n3 = 147, Power = 0.8,
                          Note = " "))
})
