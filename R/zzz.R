#  File R/zzz.R in package tergm, part of the Statnet suite
#  of packages for network analysis, https://statnet.org .
#
#  This software is distributed under the GPL-3 license.  It is free,
#  open source, and has the attribution requirements (GPL Section 7) at
#  https://statnet.org/attribution
#
#  Copyright 2008-2019 Statnet Commons
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
  ergm_proposal_table("c", "Bernoulli", "",  0, "discordTNT", "discordTNT")
  ergm_proposal_table("c", "Bernoulli", "",  0, "discordStratTNT", "discordStratTNT")
  ergm_proposal_table("c", "Bernoulli", "",  0, "discordBDTNT", "discordBDTNT")
  ergm_proposal_table("c", "Bernoulli", "",  0, "discordBDStratTNT", "discordBDStratTNT")

  ergm_proposal_table("c", "Bernoulli", "atleast",  0, "random", "formationMLE")
  ergm_proposal_table("c", "Bernoulli", "atleast+bd",  0, "random", "formationMLE")
  ergm_proposal_table("c", "Bernoulli", "atleast",  1, "TNT", "formationMLETNT")
  ergm_proposal_table("c", "Bernoulli", "atleast+bd",  1, "TNT", "formationMLETNT")

  ergm_proposal_table("c", "Bernoulli", "atmost",  0, "random", "dissolutionMLE")
  ergm_proposal_table("c", "Bernoulli", "atmost+bd",  0, "random", "dissolutionMLE")
  ergm_proposal_table("c", "Bernoulli", "atmost",  1, "TNT", "dissolutionMLETNT")
  ergm_proposal_table("c", "Bernoulli", "atmost+bd",  1, "TNT", "dissolutionMLETNT")

  ergm_proposal_table("c", "Bernoulli", "atleast+observed",  0, "random", "formationNonObservedMLE")
  ergm_proposal_table("c", "Bernoulli", "atleast+bd+observed",  0, "random", "formationNonObservedMLE")
  ergm_proposal_table("c", "Bernoulli", "atleast+observed",  1, "TNT", "formationNonObservedMLETNT")
  ergm_proposal_table("c", "Bernoulli", "atleast+bd+observed",  1, "TNT", "formationNonObservedMLETNT")
  ergm_proposal_table("c", "Bernoulli", "atmost+observed",  0, "random", "dissolutionNonObservedMLE")
  ergm_proposal_table("c", "Bernoulli", "atmost+bd+observed",  0, "random", "dissolutionNonObservedMLE")
  ergm_proposal_table("c", "Bernoulli", "atmost+observed",  1, "TNT", "dissolutionNonObservedMLETNT")
  ergm_proposal_table("c", "Bernoulli", "atmost+bd+observed",  1, "TNT", "dissolutionNonObservedMLETNT")

  ergm_proposal_table("c", "Bernoulli", "atleast+blockdiag",  0, "random", "formationMLEblockdiag")
  ergm_proposal_table("c", "Bernoulli", "atleast+bd+blockdiag",  0, "random", "formationMLEblockdiag")
  ergm_proposal_table("c", "Bernoulli", "atleast+blockdiag",  1, "TNT", "formationMLEblockdiagTNT")
  ergm_proposal_table("c", "Bernoulli", "atleast+bd+blockdiag",  1, "TNT", "formationMLEblockdiagTNT")
  ergm_proposal_table("c", "Bernoulli", "atmost+blockdiag",  0, "random", "dissolutionMLEblockdiag")
  ergm_proposal_table("c", "Bernoulli", "atmost+bd+blockdiag",  0, "random", "dissolutionMLEblockdiag")
  ergm_proposal_table("c", "Bernoulli", "atmost+blockdiag",  1, "TNT", "dissolutionMLEblockdiagTNT")
  ergm_proposal_table("c", "Bernoulli", "atmost+bd+blockdiag",  1, "TNT", "dissolutionMLEblockdiagTNT")

  ergm_proposal_table("c", "Bernoulli", "atleast+blockdiag+observed",  0, "random", "formationNonObservedMLEblockdiag")
  ergm_proposal_table("c", "Bernoulli", "atleast+bd+blockdiag+observed",  0, "random", "formationNonObservedMLEblockdiag")
  ergm_proposal_table("c", "Bernoulli", "atleast+blockdiag+observed",  1, "TNT", "formationNonObservedMLEblockdiagTNT")
  ergm_proposal_table("c", "Bernoulli", "atleast+bd+blockdiag+observed",  1, "TNT", "formationNonObservedMLEblockdiagTNT")
  ergm_proposal_table("c", "Bernoulli", "atmost+blockdiag+observed",  0, "random", "dissolutionNonObservedMLEblockdiag")
  ergm_proposal_table("c", "Bernoulli", "atmost+bd+blockdiag+observed",  0, "random", "dissolutionNonObservedMLEblockdiag")
  ergm_proposal_table("c", "Bernoulli", "atmost+blockdiag+observed",  1, "TNT", "dissolutionNonObservedMLEblockdiagTNT")
  ergm_proposal_table("c", "Bernoulli", "atmost+bd+blockdiag+observed",  1, "TNT", "dissolutionNonObservedMLEblockdiagTNT")
  
  ergm_proposal_table("f", "Bernoulli", "",  0, "random", "formation")
  ergm_proposal_table("f", "Bernoulli", "bd",  0, "random", "formation")
  ergm_proposal_table("f", "Bernoulli", "",  1, "TNT", "formationTNT")
  ergm_proposal_table("f", "Bernoulli", "bd",  1, "TNT", "formationTNT")
  ergm_proposal_table("d", "Bernoulli", "",  0, "random", "dissolution")
  ergm_proposal_table("d", "Bernoulli", "bd",  0, "random", "dissolution")
  ergm_proposal_table("d", "Bernoulli", "",  1, "TNT", "dissolutionTNT")
  ergm_proposal_table("d", "Bernoulli", "bd",  1, "TNT", "dissolutionTNT")
}
