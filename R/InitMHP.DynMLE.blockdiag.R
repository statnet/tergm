#===========================================================================
# The <InitMHP> file contains the following 24 functions for
# initializing the MHproposal object; each is prepended with 'InitMHP.'
#        <dissolutionMLE>
#        <formationNonObservedMLE>
#        <dissolutionNonObservedMLE>
#        <formationMLE>       
#============================================================================


##########################################2##############################
# Each of the <InitMHP.X> functions initializes and returns a
# MHproposal list; when appropriate, proposal types are checked against
# covariates and network types for 1 of 2 side effects: to print warning
# messages or to halt execution (only <InitMHP.nobetweengroupties> can
# halt execution)
#
# --PARAMETERS--
#   arguments: is ignored by all but <InitMHP.nobetweengroupties>,
#              where 'arguments' is used to get the nodal attributes
#              via <get.node.attr>
#   nw       : the network given by the model
#   model    : the model for 'nw', as returned by <ergm.getmodel>
#
# --RETURNED--
#   MHproposal: a list containing:
#        name   : the name of the proposal
#        inputs : a vector to be passed to the proposal
#        package: is "tergm"
#
############################################################################

InitMHP.formationMLEblockdiag <- function(arguments, nw) {
  if(is.bipartite(nw)) stop("Block-diagonal sampling is not implemented for bipartite networks at this time.")
  # rle() returns contigous runs of values.
  a <- rle(nw %v% arguments$constraints$blockdiag$attrname)
  # If we have more runs than unique values, the blocks must not be all contiguous.
  if(length(a$lengths)!=length(unique(a$values))) stop("Current implementation of block-diagonal sampling requires that the blocks be contiguous.")

  nd <- sum(a$lengths*(a$lengths-1)/(if(is.directed(nw)) 1 else 2))
  b <- cumsum(c(0,a$lengths)) # upper bounds of blocks
  w <- cumsum(a$lengths*(a$lengths-1)) # cumulative block weights ~ # dyads in the block
  w <- w/max(w)
  # Note that this automagically takes care of singleton blocks by giving them weight 0.
  
  MHproposal <- list(name = "FormationMLEblockdiag", inputs=c(ergm.Cprepare.el(arguments$constraints$atleast$nw),nd,length(b)-1,b,w))
  MHproposal
}

InitMHP.formationMLEblockdiagTNT <- function(arguments, nw) {
  if(is.bipartite(nw)) stop("Block-diagonal sampling is not implemented for bipartite networks at this time.")

  el <- as.edgelist(nw)
  a <- nw %v% arguments$constraints$blockdiag$attrname
  
  if(any(a[el[,1]]!=a[el[,2]])) stop("Block-diagonal TNT sampler implementation does not support sampling networks with off-block-diagonal ties at this time.")
  
  # rle() returns contigous runs of values.
  a <- rle(a)
  # If we have more runs than unique values, the blocks must not be all contiguous.
  if(length(a$lengths)!=length(unique(a$values))) stop("Current implementation of block-diagonal sampling requires that the blocks be contiguous.")

  nd <- sum(a$lengths*(a$lengths-1)/(if(is.directed(nw)) 1 else 2))
  b <- cumsum(c(0,a$lengths)) # upper bounds of blocks
  w <- cumsum(a$lengths*(a$lengths-1)) # cumulative block weights ~ # dyads in the block
  w <- w/max(w)
  # Note that this automagically takes care of singleton blocks by giving them weight 0.

  MHproposal <- list(name = "FormationMLEblockdiagTNT", inputs=c(ergm.Cprepare.el(arguments$constraints$atleast$nw),nd,length(b)-1,b,w))
  MHproposal
}

InitMHP.formationNonObservedMLEblockdiag <- function(arguments, nw) {
  ## Precalculate toggleable dyads: dyads which
  ## * are unobserved in y[t]
  ## * are non-ties in y[t-1]
  
  y0<-arguments$costraints$atleast$nw
  y.miss<-is.na(nw)

  if(is.bipartite(nw)) stop("Block-diagonal sampling is not implemented for bipartite networks at this time.")
  # rle() returns contigous runs of values.
  a <- nw %v% arguments$constraints$blockdiag$attrname
  # If we have more runs than unique values, the blocks must not be all contiguous.
  if(length(rle(a)$lengths)!=length(unique(rle(a)$values))) stop("Current implementation of block-diagonal sampling requires that the blocks be contiguous.")
  el <- as.edgelist(y.miss-y0)
  el <- el[a[el[,1]]==a[el[,2]],,drop=FALSE]

  ## Given the list of toggleable dyads, no formation-specific proposal function is needed:
  MHproposal <- list(name = "randomtoggleList", inputs=ergm.Cprepare.el(el), pkgname="ergm")
  MHproposal
}

InitMHP.dissolutionNonObservedMLEblockdiag <- function(arguments, nw) {
  ## Precalculate toggleable dyads: dyads which
  ## * are unobserved in y[t]
  ## * are ties in y[t-1]

  y0<-arguments$constraints$atmost$nw
  y.miss<-is.na(nw)

  if(is.bipartite(nw)) stop("Block-diagonal sampling is not implemented for bipartite networks at this time.")
  # rle() returns contigous runs of values.
  a <- nw %v% arguments$constraints$blockdiag$attrname
  # If we have more runs than unique values, the blocks must not be all contiguous.
  if(length(rle(a)$lengths)!=length(unique(rle(a)$values))) stop("Current implementation of block-diagonal sampling requires that the blocks be contiguous.")
  el <- as.edgelist(y.miss & y0)
  el <- el[a[el[,1]]==a[el[,2]],,drop=FALSE]

  
  ## Given the list of toggleable dyads, no formation-specific proposal function is needed:
  MHproposal <- list(name = "randomtoggleList", inputs=ergm.Cprepare.el(el), pkgname="ergm")
  MHproposal
}

InitMHP.dissolutionMLEblockdiag <- function(arguments, nw) {
  
  if(is.bipartite(nw)) stop("Block-diagonal sampling is not implemented for bipartite networks at this time.")

  y0<-arguments$constraints$atmost$nw
  
  # rle() returns contigous runs of values.
  a <- nw %v% arguments$constraints$blockdiag$attrname
  # If we have more runs than unique values, the blocks must not be all contiguous.
  if(length(rle(a)$lengths)!=length(unique(rle(a)$values))) stop("Current implementation of block-diagonal sampling requires that the blocks be contiguous.")
  el <- as.edgelist(y0)
  el <- el[a[el[,1]]==a[el[,2]],,drop=FALSE]

  
  MHproposal <- list(name = "randomtoggleList", inputs=ergm.Cprepare.el(el), pkgname="ergm")
  MHproposal
}

InitMHP.dissolutionMLEblockdiagTNT <- function(arguments, nw) {
  
  if(is.bipartite(nw)) stop("Block-diagonal sampling is not implemented for bipartite networks at this time.")

  y0<-arguments$constraints$atmost$nw
  
  # rle() returns contigous runs of values.
  a <- nw %v% arguments$constraints$blockdiag$attrname
  # If we have more runs than unique values, the blocks must not be all contiguous.
  if(length(rle(a)$lengths)!=length(unique(rle(a)$values))) stop("Current implementation of block-diagonal sampling requires that the blocks be contiguous.")
  el <- as.edgelist(y0)
  el <- el[a[el[,1]]==a[el[,2]],,drop=FALSE]

  # el enforces the block constraint.
  MHproposal <- list(name = "DissolutionMLETNT", inputs=ergm.Cprepare.el(el))
  MHproposal
}
