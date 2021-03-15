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
  # . is used as a placeholder by stantet.common::NVL3().
  utils::globalVariables(".")
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
  ergm_proposal_table("c", "Bernoulli", "|.dyads&sparse&discord",  1, "discordTNT", "staticDiscordTNT")
  ergm_proposal_table("t", "Bernoulli", "|discord&sparse",  1, "discordTNT", "discordTNT")
  ergm_proposal_table("t", "Bernoulli", "|strat|discord&sparse",  0, "discordStratTNT", "discordStratTNT")
  ergm_proposal_table("t", "Bernoulli", "|bdmax|blocks|discord&sparse",  0, "discordBDTNT", "discordBDTNT")
  ergm_proposal_table("t", "Bernoulli", "|bdmax|blocks|strat|discord&sparse",  0, "discordBDStratTNT", "discordBDStratTNT")
}
