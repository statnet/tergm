#  File R/is.lasttoggle.R in package tergm, part of the
#  Statnet suite of packages for network analysis, https://statnet.org .
#
#  This software is distributed under the GPL-3 license.  It is free,
#  open source, and has the attribution requirements (GPL Section 7) at
#  https://statnet.org/attribution .
#
#  Copyright 2008-2023 Statnet Commons
################################################################################
###############################################################################
# is.lasttoggle function tests whether a nw object, with the inherited model 
# info requires duration information. 
###############################################################################

#' @rdname lasttoggle
#' @name lasttoggle
#' @title Lasttoggle
#' 
#' @description A data structure used by \code{tergm} for tracking of limited information
#'                about dyad edge histories.
#' 
#' @details The \code{tergm} package handles durational information attached to
#'            [`network`] objects by way of the \code{time} and 
#'            \code{lasttoggle} network attributes.  The \code{lasttoggle} data
#'            structure is a 3-column matrix; the first two columns are tails 
#'            and heads (respectively) of dyads, and the third column is the last
#'            time at which the dyad was toggled.  The default last toggle time
#'            is \code{-INT_MAX/2}.  Last toggle times for non-edges are 
#'            periodically cleared in the C code.  The \code{time} network 
#'            attribute is simply an integer, and together with the 
#'            \code{lasttoggle} data it determines the age of an extant
#'            tie as \code{time + 1} minus the last toggle time for that dyad.
#'            The default value for \code{time} is 0. 
NULL

is.lasttoggle <- function(nw, formation=NULL, dissolution=NULL, monitor=NULL, targets=NULL) {  
  if(!is.null(formation))
    formation <- nonsimp_update.formula(formation, nw~., from.new="nw")
  
  if(!is.null(dissolution))  
    dissolution <- nonsimp_update.formula(dissolution, nw~., from.new="nw")

  if(is(monitor, "formula"))
    monitor <- nonsimp_update.formula(monitor, nw~., from.new="nw")
     
  if(is(targets, "formula"))
    targets <- nonsimp_update.formula(targets, nw~., from.new="nw")
  
  as.integer(is.durational(formation) || is.durational(dissolution) || is.durational(monitor) || is.durational(targets))
}
