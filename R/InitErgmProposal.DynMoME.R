#  File R/InitErgmProposal.DynMoME.R in package tergm, part of the
#  Statnet suite of packages for network analysis, https://statnet.org .
#
#  This software is distributed under the GPL-3 license.  It is free,
#  open source, and has the attribution requirements (GPL Section 7) at
#  https://statnet.org/attribution .
#
#  Copyright 2008-2023 Statnet Commons
################################################################################

#' @templateVar name discordTNT
#' @aliases InitErgmProposal.discordTNT
#' @title Temperal TNT proposal
#' @description A temporal version of \code{\link[ergm:TNT-ergmProposal]{TNT}}, with approximately
#'   \code{discordance_fraction} of proposed toggles being made on the set of discordant dyads,
#'   approximately \code{1 - discordance_fraction} of proposed toggles being TNT proposals from
#'   network.  The value of \code{discordance_fraction} can be set by the user as a proposal argument,
#'   defaults to `0.5`.
#' @template ergmProposal-general
NULL

InitErgmProposal.discordTNT <- function(arguments, nw, ...) {
  discordance_fraction <- NVL(arguments$discordance_fraction, 1/2)
  if(!is.numeric(discordance_fraction) || length(discordance_fraction) != 1 || discordance_fraction < 0 || discordance_fraction >= 1) {
    ergm_Init_abort("Argument ", sQuote("discordance_fraction"), " to ", sQuote("discordTNT"), " must be a number in [0,1).")
  }
  
  proposal <- list(name = "discordTNT", inputs=NULL, auxiliaries = trim_env(~.lasttoggle), discordance_fraction = discordance_fraction)
  proposal
}

#' @templateVar name discordBDStratTNT
#' @aliases InitErgmProposal.discordBDStratTNT
#' @title Temperal TNT proposal with degree bounds
#' @description A temporal version of \code{\link[ergm:ergm-proposals]{BDStratTNT}}.  Within each
#'   mixing type, approximately 50\% of proposed toggles are made on
#'   dyads, and approximately 50\% of proposed toggles are
#'   proposals from the network, all subject to the bounded degree
#'   and mixing type constraints.  The degree bound constraint is imposed
#'   the instantaneous network state
#'   rather than the temporal operator networks).
#'
#'   arguments are the same as for \code{\link[ergm:BDStratTNT-ergmProposal]{BDStratTNT}},
#'   and should be passed in via the \code{\link[ergm:bd-ergmConstraint]{bd}} and
#'   \code{\link[ergm:blocks-ergmConstraint]{blocks}} constraints and
#'   \code{\link[ergm:strat-ergmHint]{strat}} hint.
#' @template ergmProposal-general
NULL

InitErgmProposal.discordBDStratTNT <- function(arguments, nw, ...) {
  # Work around CRAN's ::: warning.
  InitErgmProposal.BDStratTNT <- eval(locate_prefixed_function("BDStratTNT", "InitErgmProposal", "Metropolis-Hastings proposal"))
  proposal <- InitErgmProposal.BDStratTNT(arguments, nw, ...)
  proposal$name <- "discordBDStratTNT"
  proposal$auxiliaries <- trim_env(~.lasttoggle)

  discordance_fraction <- NVL(arguments$discordance_fraction, 1/2)
  if(!is.numeric(discordance_fraction) || length(discordance_fraction) != 1 || discordance_fraction < 0 || discordance_fraction >= 1) {
    ergm_Init_abort("Argument ", sQuote("discordance_fraction"), " to ", sQuote("discordBDStratTNT"), " must be a number in [0,1).")
  }
  proposal$discordance_fraction <- as.numeric(discordance_fraction)
  
  proposal
}
