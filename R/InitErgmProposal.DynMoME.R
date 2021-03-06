#  File R/InitErgmProposal.DynMoME.R in package tergm, part of the
#  Statnet suite of packages for network analysis, https://statnet.org .
#
#  This software is distributed under the GPL-3 license.  It is free,
#  open source, and has the attribution requirements (GPL Section 7) at
#  https://statnet.org/attribution .
#
#  Copyright 2008-2021 Statnet Commons
################################################################################
#===================================================================
# This file contains the 5 following MHP initializers, each
# prepended with 'InitErgmProposal.'  All of these functions may also be
# found in the <InitErgmProposal> file.
#      <formation>       <formationTNT>
#      <dissolution>
#===================================================================

InitErgmProposal.discordTNT <- function(arguments, nw, ...) {
  discordance_fraction <- NVL(arguments$discordance_fraction, 1/2)
  if(!is.numeric(discordance_fraction) || length(discordance_fraction) != 1 || discordance_fraction < 0 || discordance_fraction >= 1) {
    ergm_Init_abort("Argument ", sQuote("discordance_fraction"), " to ", sQuote("discordTNT"), " must be a number in [0,1).")
  }
  
  proposal <- list(name = "discordTNT", inputs=NULL, auxiliaries = ~.lasttoggle, discordance_fraction = discordance_fraction)
  proposal
}

InitErgmProposal.discordBDStratTNT <- function(arguments, nw, ...) {
  # Work around CRAN's ::: warning.
  InitErgmProposal.BDStratTNT <- eval(locate_prefixed_function("BDStratTNT", "InitErgmProposal", "Metropolis-Hastings proposal"))
  proposal <- InitErgmProposal.BDStratTNT(arguments, nw, ...)
  proposal$name <- "discordBDStratTNT"
  proposal$auxiliaries <- ~.lasttoggle

  discordance_fraction <- NVL(arguments$discordance_fraction, 1/2)
  if(!is.numeric(discordance_fraction) || length(discordance_fraction) != 1 || discordance_fraction < 0 || discordance_fraction >= 1) {
    ergm_Init_abort("Argument ", sQuote("discordance_fraction"), " to ", sQuote("discordBDStratTNT"), " must be a number in [0,1).")
  }
  proposal$discordance_fraction <- as.numeric(discordance_fraction)
  
  proposal
}
