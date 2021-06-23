#  File R/InitErgmProposal.DynMLE.R in package tergm, part of the
#  Statnet suite of packages for network analysis, https://statnet.org .
#
#  This software is distributed under the GPL-3 license.  It is free,
#  open source, and has the attribution requirements (GPL Section 7) at
#  https://statnet.org/attribution .
#
#  Copyright 2008-2021 Statnet Commons
################################################################################
InitErgmProposal.staticDiscordTNT <- function(arguments, nw, model) {
  dissolvable <- as.rlebdm(arguments$constraints$discord$nw)
  formable <- !dissolvable
  dissolvable <- ergm_dyadgen_select(arguments, nw, dissolvable)
  formable <- ergm_dyadgen_select(arguments, nw, formable)
  bd <- ergm_bd_init(arguments, nw)
  list(name = "staticDiscordTNT", formable = formable, dissolvable = dissolvable, inputs = 0.5, bd = bd)
}
