#  File R/zzz.R in package tergm, part of the Statnet suite
#  of packages for network analysis, https://statnet.org .
#
#  This software is distributed under the GPL-3 license.  It is free,
#  open source, and has the attribution requirements (GPL Section 7) at
#  https://statnet.org/attribution
#
#  Copyright 2008-2020 Statnet Commons
#######################################################################
#' @import statnet.common
.onAttach <- function(lib, pkg){
  sm <- statnetStartupMessage("tergm", c("statnet"), FALSE)
  if(!is.null(sm)) packageStartupMessage(sm)
}

#' @import ergm
.onLoad <- function(lib, pkg){
  .RegisterProposals()
}

.RegisterProposals <- function(){
  ergm_proposal_table("c", "Bernoulli", "|.dyads&TNT&discord",  1, "discordTNT", "staticDiscordTNT")
  ergm_proposal_table("t", "Bernoulli", "|discord&TNT",  1, "discordTNT", "discordTNT")
  ergm_proposal_table("t", "Bernoulli", "|Strat|discord&TNT",  0, "discordStratTNT", "discordStratTNT")
  ergm_proposal_table("t", "Bernoulli", "|bdmax|blocks|discord&TNT",  0, "discordBDTNT", "discordBDTNT")
  ergm_proposal_table("t", "Bernoulli", "|bdmax|blocks|Strat|discord&TNT",  0, "discordBDStratTNT", "discordBDStratTNT")
}
