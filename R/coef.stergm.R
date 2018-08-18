#  File R/coef.stergm.R in package tergm, part of the Statnet suite
#  of packages for network analysis, http://statnet.org .
#
#  This software is distributed under the GPL-3 license.  It is free,
#  open source, and has the attribution requirements (GPL Section 7) at
#  http://statnet.org/attribution
#
#  Copyright 2008-2018 Statnet Commons
#######################################################################


#' @describeIn stergm Extract parameter estimates.
#' 
#' 
#' @param object A \code{\link{stergm}} fit.
#' @return `coef` and `coefficients` methods return parameter
#'   estimates extracted from \code{object} in the form of a list with
#'   two elements: \code{formation}, a vector of formation
#'   coefficients and \code{dissolution}, a vector of dissolution
#'   coefficients.
#' @keywords regression models
#' @importFrom stats coef
#' @export
coef.stergm <- function(object, ...){list(formation=object$formation.fit$coef,
                                          dissolution=object$dissolution.fit$coef)}

#' @describeIn stergm An \emph{alias} for the `coef` method.
#' @importFrom stats coefficients
#' @export
coefficients.stergm <- coef.stergm
