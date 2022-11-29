test_that("oneprop_ss() returns correct n", {
  expect_equal(oneprop_ss(p0 = 0.2, pA = 0.3, alpha = 0.05, power = 0.8,
            method = "conditional", one.or.two.sided = "one"), 109)
})
