#  File R/is.durational.R in package tergm, part of the Statnet suite
#  of packages for network analysis, https://statnet.org .
#
#  This software is distributed under the GPL-3 license.  It is free,
#  open source, and has the attribution requirements (GPL Section 7) at
#  https://statnet.org/attribution
#
#  Copyright 2008-2020 Statnet Commons
#######################################################################
###############################################################################
# These functions are used to detect whether a ERGM formula/model/etcs are 
# durational dependent or not, based on the (T)ERGM term used.
# To make an (T)ERGM term durational dependent, simply add an "duration" object
# to terms in InitErgmTerm.duration.R
###############################################################################


#' Testing for duration dependent models
#' 
#' These functions test whether an ERGM is duration dependent or not.
#' 
#' @param object An ERGM formula, \code{\link{ergm_model}} object, or
#'                 \code{\link{ergm_state}} object.
#' @param \dots Unused at this time.
#' @return \code{TRUE} if the ERGM terms in the model are duration dependent; 
#'           \code{FALSE} otherwise.
#' @keywords model
#' @export
is.durational<-function(object,...) UseMethod("is.durational")

#' @rdname is.durational
#' @description The method for `NULL` always returns `FALSE` by
#'   convention.
#' @export
is.durational.NULL <- function(object, ...) FALSE # By convention.

#' @describeIn is.durational Test if the \code{\link{ergm_model}} has duration-dependent terms, which call for \code{\link{lasttoggle}} data structures.
#' @export
is.durational.ergm_model <- function(object, ...){
#' @import purrr
  map(object$terms, "duration") %>% unlist %>% NVL(FALSE) %>% max
}

#' @describeIn is.durational Test if the \code{\link{ergm_state}} has duration-dependent terms, which call for \code{\link{lasttoggle}} data structures.
#' @export
is.durational.ergm_state <- function(object, ...){
  is.durational(object$model)
}

#' @rdname is.durational
#' @param response,basis See [ergm()].
#' @export
is.durational.formula<-function(object,response=NULL,basis=ergm.getnetwork(object),...){
  # If basis is not null, replace network in formula by basis.
  # In either case, let nw be network object from formula.
  
  nw <- basis
  ergm_preprocess_response(nw,response)
  m<-ergm_model(object, nw, dynamic=TRUE, ...)
  is.durational(m)
}
