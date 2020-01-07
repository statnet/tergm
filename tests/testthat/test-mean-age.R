#  File tests/testthat/test-mean-age.R in package tergm, part of the Statnet suite
#  of packages for network analysis, https://statnet.org .
#
#  This software is distributed under the GPL-3 license.  It is free,
#  open source, and has the attribution requirements (GPL Section 7) at
#  https://statnet.org/attribution
#
#  Copyright 2003-2020 Statnet Commons
#######################################################################

context("test-mean-age.R")

test_that("mean.age behaves correctly with summary and godfather", {
  
  el <- matrix(c(1, 2,
                 1, 3,
                 2, 3,
                 5, 6,
                 7, 8,
                 7, 10), ncol = 2, byrow = TRUE)
  
  init.lasttoggle.times <- c(0,1,2,3,3,5)
  
  nw <- network(el, dir = FALSE)
  
  nw %n% "time" <- 5
  
  nw %n% "lasttoggle" <- cbind(el, init.lasttoggle.times)
  
  expect_equal(unname(summary(nw ~ mean.age)), mean(5 - init.lasttoggle.times + 1))
  
  toggles <- matrix(c(6, 2, 3,
                      7, 1, 4,
                      7, 1, 5,
                      8, 7, 8,
                      8, 7, 10,
                      8, 9, 10,
                      9, 1, 4,
                      9, 5, 6,
                      10, 5, 6,
                      10, 5, 6,
                      11, 1, 2,
                      11, 1, 2,
                      12, 5, 8,
                      12, 5, 8,
                      12, 5, 8,
                      13, 1, 3,
                      13, 1, 3,
                      13, 1, 3), ncol = 3, byrow = TRUE)
  
  rv <- tergm.godfather(nw ~ edges + mean.age, toggles=toggles)
  rv <- as.matrix(rv)
  
  lt.times.6 <- init.lasttoggle.times[-3]
  lt.times.7 <- c(lt.times.6, 7, 7)
  lt.times.8 <- c(lt.times.7[-c(4,5)], 8)
  lt.times.9 <- c(lt.times.8[-c(3,4)])
  lt.times.10 <- lt.times.9
  lt.times.11 <- lt.times.10
  lt.times.12 <- c(lt.times.11, 12)
  lt.times.13 <- lt.times.12[-2]
  
  
  expect_equal(mean(6 - lt.times.6 + 1), unname(rv[1,2]))
  expect_equal(mean(7 - lt.times.7 + 1), unname(rv[2,2]))
  expect_equal(mean(8 - lt.times.8 + 1), unname(rv[3,2]))
  expect_equal(mean(9 - lt.times.9 + 1), unname(rv[4,2]))
  expect_equal(mean(10 - lt.times.10 + 1), unname(rv[5,2]))
  expect_equal(mean(11 - lt.times.11 + 1), unname(rv[6,2]))
  expect_equal(mean(12 - lt.times.12 + 1), unname(rv[7,2]))
  expect_equal(mean(13 - lt.times.13 + 1), unname(rv[8,2]))
  
})


test_that("mean.(log)age behaves correctly with summary and godfather", {
  
  el <- matrix(c(1, 2,
                 1, 3,
                 2, 3,
                 5, 6,
                 7, 8,
                 7, 10), ncol = 2, byrow = TRUE)
  
  init.lasttoggle.times <- c(0,1,2,3,3,5)
  
  nw <- network(el, dir = FALSE)
  
  nw %n% "time" <- 5
  
  nw %n% "lasttoggle" <- cbind(el, init.lasttoggle.times)
  
  expect_equal(unname(summary(nw ~ mean.age(log=TRUE))), mean(log(5 - init.lasttoggle.times + 1)))
  
  toggles <- matrix(c(6, 2, 3,
                      7, 1, 4,
                      7, 1, 5,
                      8, 7, 8,
                      8, 7, 10,
                      8, 9, 10,
                      9, 1, 4,
                      9, 5, 6,
                      10, 5, 6,
                      10, 5, 6,
                      11, 1, 2,
                      11, 1, 2,
                      12, 5, 8,
                      12, 5, 8,
                      12, 5, 8,
                      13, 1, 3,
                      13, 1, 3,
                      13, 1, 3), ncol = 3, byrow = TRUE)
  
  rv <- tergm.godfather(nw ~ edges + mean.age(log=TRUE), toggles=toggles)
  rv <- as.matrix(rv)
  
  lt.times.6 <- init.lasttoggle.times[-3]
  lt.times.7 <- c(lt.times.6, 7, 7)
  lt.times.8 <- c(lt.times.7[-c(4,5)], 8)
  lt.times.9 <- c(lt.times.8[-c(3,4)])
  lt.times.10 <- lt.times.9
  lt.times.11 <- lt.times.10
  lt.times.12 <- c(lt.times.11, 12)
  lt.times.13 <- lt.times.12[-2]
  
  
  expect_equal(mean(log(6 - lt.times.6 + 1)), unname(rv[1,2]))
  expect_equal(mean(log(7 - lt.times.7 + 1)), unname(rv[2,2]))
  expect_equal(mean(log(8 - lt.times.8 + 1)), unname(rv[3,2]))
  expect_equal(mean(log(9 - lt.times.9 + 1)), unname(rv[4,2]))
  expect_equal(mean(log(10 - lt.times.10 + 1)), unname(rv[5,2]))
  expect_equal(mean(log(11 - lt.times.11 + 1)), unname(rv[6,2]))
  expect_equal(mean(log(12 - lt.times.12 + 1)), unname(rv[7,2]))
  expect_equal(mean(log(13 - lt.times.13 + 1)), unname(rv[8,2]))
  
})

