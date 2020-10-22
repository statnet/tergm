#  File R/InitErgmConstraint.hints.R in package tergm, part of the Statnet suite
#  of packages for network analysis, https://statnet.org .
#
#  This software is distributed under the GPL-3 license.  It is free,
#  open source, and has the attribution requirements (GPL Section 7) at
#  https://statnet.org/attribution
#
#  Copyright 2008-2020 Statnet Commons
#######################################################################
InitErgmConstraint.discordTNT<-function(lhs.nw, ref, ...){
  nw <- if(is.character(ref)) lhs.nw %n% ref else lhs.nw

  if(length(list(...)))
     ergm_Init_abort(paste("discordTNT hint only takes one arguments at this time."))
   list(dependence = FALSE, priority=10, nw=nw)
}
