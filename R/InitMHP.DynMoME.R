#  File R/InitMHP.DynMoME.R in package tergm, part of the Statnet suite
#  of packages for network analysis, https://statnet.org .
#
#  This software is distributed under the GPL-3 license.  It is free,
#  open source, and has the attribution requirements (GPL Section 7) at
#  https://statnet.org/attribution
#
#  Copyright 2008-2019 Statnet Commons
#######################################################################
#===================================================================
# This file contains the 5 following MHP initializers, each
# prepended with 'InitErgmProposal.'  All of these functions may also be
# found in the <InitErgmProposal> file.
#      <formation>       <formationTNT>
#      <dissolution>
#===================================================================

InitErgmProposal.discordTNT <- function(arguments, nw, model) {
  proposal <- list(name = "discordTNT", inputs=NULL, auxiliaries = ~.lasttoggle)
  proposal
}

InitErgmProposal.discordStratTNT <- function(arguments, nw, model) {
  stp <- ergm:::InitErgmProposal.StratTNT(arguments, nw)
  proposal <- list(name = "discordStratTNT", inputs=stp$inputs, auxiliaries = ~.lasttoggle)
  proposal
}

InitErgmProposal.discordBDTNT <- function(arguments, nw, model) {
  stp <- ergm:::InitErgmProposal.BDTNT(arguments, nw)
  proposal <- list(name = "discordBDTNT", inputs=stp$inputs, auxiliaries = ~.lasttoggle)
  proposal
}

InitErgmProposal.discordBDStratTNT <- function(arguments, nw, model) {
  stp <- ergm:::InitErgmProposal.BDStratTNT(arguments, nw)
  proposal <- list(name = "discordBDStratTNT", inputs=stp$inputs, auxiliaries = ~.lasttoggle)
  proposal
}

InitErgmProposal.formation <- function(arguments, nw, model) {
  proposal <- list(name = "Formation", inputs=NULL)
  proposal
}

InitErgmProposal.formationTNT <- function(arguments, nw, model) {
  proposal <- list(name = "FormationTNT", inputs=NULL)
  proposal
}

InitErgmProposal.dissolution <- function(arguments, nw, model) {
  proposal <- list(name = "Dissolution", inputs=NULL)
  proposal
}

InitErgmProposal.dissolutionTNT <- function(arguments, nw, model) {
  proposal <- list(name = "DissolutionTNT", inputs=NULL)
  proposal
}
