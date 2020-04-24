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

#' @describeIn lasttoggle
#' @name lasttoggle
#' @title Lasttoggle placeholder RD
#' @descrption This needs to be written (or links to it removed) before the 4.0 release.
NULL

is.lasttoggle <- function(nw,formation=NULL,dissolution=NULL,monitor=NULL,target=NULL){
  
  if(!is.null(formation))
    formation<-nonsimp_update.formula(formation,nw~., from.new="nw")
  
  if(!is.null(dissolution))  
    dissolution<-nonsimp_update.formula(dissolution,nw~., from.new="nw")
  
  if(!is.null(monitor)){    
    if(is.character(monitor)){
      monitor <- switch(monitor,
          formation = formation,
          dissolution = dissolution,
          all = append_rhs.formula(~nw, unique(lapply(c(list_rhs.formula(formation),list_rhs.formula(dissolution)), unset.offset.call)))
      )
    }
    
    if(!is.null(monitor)) 
      monitor <- nonsimp_update.formula(monitor,nw~., from.new="nw")
  }
  
  
  if(!is.null(target)){
    if(is.character(targets)){
      targets <- switch(targets,
          formation = formation,
          dissolution = dissolution)}
      
      targets <- nonsimp_update.formula(targets,nw~., from.new="nw")
    }
  
  
  
    duration.dependent <- if(is.durational(formation) || is.durational(dissolution)|| is.durational(monitor))
        {1} else {0}
    
    duration.dependent
    
  }
