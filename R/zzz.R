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
.onAttach <- function(libname, pkgname){
  sm <- statnetStartupMessage("tergm", c("statnet"), FALSE)
  if(!is.null(sm)) packageStartupMessage(sm)
}

#' @import ergm
.onLoad <- function(libname, pkgname){
  eval(COLLATE_ALL_MY_CONTROLS_EXPR)

  .RegisterProposals()
}

# TODO: Figure out some automatic way to keep this in sync with statnet.common.
#' @name snctrl
#'
#' @title Statnet Control
#'
#' @description A utility to facilitate argument completion of control lists, reexported from `statnet.common`.
#'
#' @section Currently recognised control parameters:
#' This list is updated as packages are loaded and unloaded.
#'
#' \Sexpr[results=rd,stage=render]{statnet.common::snctrl_names()}
#'
#' @seealso [statnet.common::snctrl()]
#' @docType import
NULL
#' @export
snctrl <- statnet.common::snctrl
## BEGIN boilerplate: should be kept in sync with statnet.common.


eval(UPDATE_MY_SCTRL_EXPR)

.RegisterProposals <- function(){
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
