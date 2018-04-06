#  File R/InitErgmProposal.DynMLE.R in package tergm, part of the Statnet suite
#  of packages for network analysis, http://statnet.org .
#
#  This software is distributed under the GPL-3 license.  It is free,
#  open source, and has the attribution requirements (GPL Section 7) at
#  http://statnet.org/attribution
#
#  Copyright 2008-2017 Statnet Commons
#######################################################################
#===========================================================================
# The <InitErgmProposal> file contains the following 24 functions for
# initializing the proposal object; each is prepended with 'InitErgmProposal.'
#        <dissolutionMLE>
#        <formationNonObservedMLE>
#        <dissolutionNonObservedMLE>
#        <formationMLE>       
#============================================================================


##########################################2##############################
# Each of the <InitErgmProposal.X> functions initializes and returns a
# proposal list; when appropriate, proposal types are checked against
# covariates and network types for 1 of 2 side effects: to print warning
# messages or to halt execution (only <InitErgmProposal.nobetweengroupties> can
# halt execution)
#
# --PARAMETERS--
#   arguments: is ignored by all but <InitErgmProposal.nobetweengroupties>,
#              where 'arguments' is used to get the nodal attributes
#              via <get.node.attr>
#   nw       : the network given by the model
#   model    : the model for 'nw', as returned by <ergm_model>
#
# --RETURNED--
#   proposal: a list containing:
#        name   : the name of the proposal
#        inputs : a vector to be passed to the proposal
#        package: is "tergm"
#
############################################################################

InitErgmProposal.formationMLE <- function(arguments, nw) {
  proposal <- list(name = "FormationMLE", inputs=to_ergm_Cdouble(arguments$constraints$atleast$nw))
  proposal
}

InitErgmProposal.formationMLETNT <- function(arguments, nw) {
  proposal <- list(name = "FormationMLETNT", inputs=to_ergm_Cdouble(arguments$constraints$atleast$nw))
  proposal
}

InitErgmProposal.dissolutionMLE <- function(arguments, nw) {
  proposal <- list(name = "randomtoggleList", inputs=to_ergm_Cdouble(arguments$constraints$atmost$nw), pkgname="ergm")
  proposal
}

InitErgmProposal.dissolutionMLETNT <- function(arguments, nw) {
  proposal <- list(name = "DissolutionMLETNT", inputs=to_ergm_Cdouble(arguments$constraints$atmost$nw))
  proposal
}

InitErgmProposal.formationNonObservedMLE <- function(arguments, nw) {
  ## Precalculate toggleable dyads: dyads which
  ## * are unobserved in y[t]
  ## * are non-ties in y[t-1]
  
  y0<-arguments$costraints$atleast$nw
  y.miss<-is.na(nw)

  ## Given the list of toggleable dyads, no formation-specific proposal function is needed:
  proposal <- list(name = "randomtoggleList", inputs=to_ergm_Cdouble(y.miss-y0), pkgname="ergm")
  proposal
}

InitErgmProposal.formationNonObservedMLETNT <- function(arguments, nw) {
  ## Precalculate toggleable dyads: dyads which
  ## * are unobserved in y[t]
  ## * are non-ties in y[t-1]
  
  y0<-arguments$costraints$atleast$nw
  y.miss<-is.na(nw)

  ## Given the list of toggleable dyads, no formation-specific proposal function is needed:
  proposal <- list(name = "listTNT", inputs=to_ergm_Cdouble(y.miss-y0), pkgname="ergm")
  proposal
}

InitErgmProposal.dissolutionNonObservedMLE <- function(arguments, nw) {
  ## Precalculate toggleable dyads: dyads which
  ## * are unobserved in y[t]
  ## * are ties in y[t-1]

  y0<-arguments$constraints$atmost$nw
  y.miss<-is.na(nw)

  ## Given the list of toggleable dyads, no formation-specific proposal function is needed:
  proposal <- list(name = "randomtoggleList", inputs=to_ergm_Cdouble(y.miss & y0), pkgname="ergm")
  proposal
}

InitErgmProposal.dissolutionNonObservedMLETNT <- function(arguments, nw) {
  ## Precalculate toggleable dyads: dyads which
  ## * are unobserved in y[t]
  ## * are ties in y[t-1]

  y0<-arguments$constraints$atmost$nw
  y.miss<-is.na(nw)

  ## Given the list of toggleable dyads, no formation-specific proposal function is needed:
  proposal <- list(name = "listTNT", inputs=to_ergm_Cdouble(y.miss & y0), pkgname="ergm")
  proposal
}

