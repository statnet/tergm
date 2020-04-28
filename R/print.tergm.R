#  File R/print.stergm.R in package tergm, part of the Statnet suite
#  of packages for network analysis, https://statnet.org .
#
#  This software is distributed under the GPL-3 license.  It is free,
#  open source, and has the attribution requirements (GPL Section 7) at
#  https://statnet.org/attribution
#
#  Copyright 2008-2019 Statnet Commons
#######################################################################


#' @describeIn tergm Print the parameter estimates.
#' 
#' @param x A \code{\link{tergm}} object.
#' @param digits Significant digits for coefficients
#' @export
print.tergm <- function (x, digits = max(3, getOption("digits") - 3), ...) {
  cat("Coefficients:\n")
  print.default(format(x$coef, digits = digits), print.gap = 2, quote = FALSE)
}

