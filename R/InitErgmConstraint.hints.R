InitErgmConstraint.discordTNT<-function(lhs.nw, ref, ...){
  nw <- if(is.character(ref)) lhs.nw %n% ref else lhs.nw

  if(length(list(...)))
     ergm_Init_abort(paste("discordTNT hint only takes one arguments at this time."))
   list(dependence = FALSE, priority=10, impliedby=c("edges", "degrees", "edges", "idegrees", "odegrees", "b1degrees", "b2degrees", "idegreedist", "odegreedist", "degreedist", "b1degreedist", "b2degreedist"), nw=nw)
}
