InitErgmConstraint.discordTNT<-function(lhs.nw, ref, ...){
  nw <- if(is.character(ref)) lhs.nw %n% ref else lhs.nw

  if(length(list(...)))
     ergm_Init_abort(paste("discordTNT hint only takes one arguments at this time."))
   list(dependence = FALSE, priority=10, nw=nw)
}
  
InitErgmConstraint.discord <- function(lhs.nw, ...) {
   if(length(list(...)))
     ergm_Init_abort(paste("discord hint does not take arguments at this time."))
   list(dependence = FALSE, priority=10)
}
