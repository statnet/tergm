#  File R/is.lasttoggle.R in package tergm, part of the Statnet suite
#  of packages for network analysis, https://statnet.org .
#
#  This software is distributed under the GPL-3 license.  It is free,
#  open source, and has the attribution requirements (GPL Section 7) at
#  https://statnet.org/attribution
#
#  Copyright 2008-2019 Statnet Commons
#######################################################################
###############################################################################
# is.lasttoggle function tests whether a nw object, with the inherited model 
# info requires duration information. 
###############################################################################

#' @rdname lasttoggle
#' @name lasttoggle
#' @title Lasttoggle placeholder RD
#' 
#' This needs to be written (or links to it removed) before the 4.0 release.
NULL

is.lasttoggle <- function(nw, formation=NULL, dissolution=NULL, monitor=NULL, targets=NULL) {  
  if(!is.null(formation))
    formation<-nonsimp_update.formula(formation,nw~., from.new="nw")
  
  if(!is.null(dissolution))  
    dissolution<-nonsimp_update.formula(dissolution,nw~., from.new="nw")

  if(is(monitor, "formula"))
    monitor <- nonsimp_update.formula(monitor, nw~., from.new="nw")
     
  if(is(targets, "formula"))
    targets <- nonsimp_update.formula(targets, nw~., from.new="nw")
  
  as.integer(is.durational(formation) || is.durational(dissolution) || is.durational(monitor) || is.durational(targets))
}
