#  File tests/testthat/test-summary-networkDynamic.R in package tergm, part of the
#  Statnet suite of packages for network analysis, https://statnet.org .
#
#  This software is distributed under the GPL-3 license.  It is free,
#  open source, and has the attribution requirements (GPL Section 7) at
#  https://statnet.org/attribution .
#
#  Copyright 2008-2024 Statnet Commons
################################################################################

test_that("summary_formula.networkDynamic behaves reasonably", {
  nw <- network.initialize(10, directed = FALSE)
  nwd <- simulate(nw ~ edges, coef = c(0), time.slices = 10, dynamic = TRUE)
  s1 <- summary(nwd ~ edges, at = c(1,3,6,8))
  s2 <- summary(nwd ~ edges + triangle, at = c(1,3,6,8))
  s3 <- summary(nwd ~ edges, at = c(3))
  s4 <- summary(nwd ~ edges + triangle, at = c(2))
  expect_true(is(s1, "matrix"))
  expect_true(is(s2, "matrix"))
  expect_true(is(s3, "matrix"))
  expect_true(is(s4, "matrix"))
  expect_identical(colnames(s1), c("edges"))
  expect_identical(colnames(s2), c("edges", "triangle"))
  expect_identical(colnames(s3), c("edges"))
  expect_identical(colnames(s4), c("edges", "triangle"))
  expect_identical(rownames(s1), NULL)
  expect_identical(rownames(s2), NULL)
  expect_identical(rownames(s3), NULL)
  expect_identical(rownames(s4), NULL)
})
