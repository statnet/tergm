#  File tests/testthat/test-basis.R in package tergm, part of the
#  Statnet suite of packages for network analysis, https://statnet.org .
#
#  This software is distributed under the GPL-3 license.  It is free,
#  open source, and has the attribution requirements (GPL Section 7) at
#  https://statnet.org/attribution .
#
#  Copyright 2008-2022 Statnet Commons
################################################################################

test_that("basis works as expected for CMLE", {
  data(samplk)

  net1 <- list(samplk1, samplk2)
  net2 <- list(samplk2, samplk3)

  set.seed(0)
  sampfit1 <- tergm(net1 ~ Form(~edges + mutual + transitiveties + cyclicalties) +
                    Diss(~edges + mutual + transitiveties + cyclicalties),
                    eval.loglik = TRUE,
                    estimate = "CMLE")

  set.seed(0)
  sampfit2 <- tergm(~Form(~edges + mutual + transitiveties + cyclicalties) +
                    Diss(~edges + mutual + transitiveties + cyclicalties),
                    eval.loglik = TRUE,
                    basis = net1,
                    estimate = "CMLE")

  set.seed(0)
  sampfit3 <- tergm(net2 ~ Form(~edges + mutual + transitiveties + cyclicalties) +
                    Diss(~edges + mutual + transitiveties + cyclicalties),
                    eval.loglik = TRUE,
                    basis = net1,
                    estimate = "CMLE")

  expect_equal(coef(sampfit1), coef(sampfit2))
  expect_equal(coef(sampfit2), coef(sampfit3))
  expect_equal(logLik(sampfit1), logLik(sampfit2))
  expect_equal(logLik(sampfit2), logLik(sampfit3))
})
