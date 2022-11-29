test_that("anova4_ss() outputs proper sample size and power", {
  expect_equal(anova4_ss(means = c(5, 10, 12, 15), sds = c(10, 10, 10, 10), alphalevel = 0.05,
                                 powerlevel = 0.8, groupratios = c(2, 1, 1, 2), Nmin = 60, Nmax = 120),
        data.frame(n1 = 24, n2 = 12, n3 = 12, n4=24, Power = 0.83,
                          Note = " "))
})
