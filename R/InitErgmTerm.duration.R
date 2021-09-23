#  File R/InitErgmTerm.duration.R in package tergm, part of the
#  Statnet suite of packages for network analysis, https://statnet.org .
#
#  This software is distributed under the GPL-3 license.  It is free,
#  open source, and has the attribution requirements (GPL Section 7) at
#  https://statnet.org/attribution .
#
#  Copyright 2008-2021 Statnet Commons
################################################################################

#' @name edges.ageinterval-ergmTerm
#' @title Number of edges with age falling into a specified range
#' @description Number of edges with age falling into a specified range
#' @details This term counts the number of edges in the network for
#'   which the time elapsed since formation is greater than or equal to
#'   `from` but strictly less than `to` . In other words, it
#'   is in the semiopen interval `[from, to)` .
#'
#' @usage
#' # binary: edges.ageinterval(from, to=+Inf)
#' @param from,to parameters to specify the lower bound and strict upper bounds. Can be scalars, vectors of the same length, or one of them must have length one, in which case it is recycled.
#'
#' @template ergmTerm-general
#'
#' @concept durational
InitErgmTerm.edges.ageinterval<-function(nw, arglist, ...) {
  stopifnot_dynamic(nw, ...)
  a <- check.ErgmTerm(nw, arglist,
                      varnames = c("from","to"),
                      vartypes = c("numeric","numeric"),
                      defaultvalues = list(NULL, Inf),
                      required = c(TRUE, FALSE))
     
  from<-a$from
  to<-a$to

  if(length(to)==1 && length(from)>1) to <- rep(to, length(from))
  else if(length(from)==1 && length(to)>1) from <- rep(from, length(to))
  else if(length(from)!=length(to)) stop("The arguments of term edges.ageinterval must have arguments either of the same length, or one of them must have length 1.")
else if(any(from>=to)) stop("Term edges.ageinterval must have from<to.")

  if(any(from<1)) stop("An extant edge cannot have an \"age\" of less than 1.")
  list(name="edges_ageinterval",
       coef.names = paste("edges.age",from,"to",to,sep=""),
       inputs=c(rbind(from, ifelse(to==Inf, 0, to))),
       duration=TRUE,
       dependence=FALSE,
       auxiliaries = ~.lasttoggle)
}

#' @name edge.ages-ergmTerm
#' @title Sum of ages of extant ties
#' @description Sum of ages of extant ties
#' @details This term adds one statistic equaling sum, over all ties
#'   present in the network, of the amount of time elapsed since
#'   formation.
#'
#'   Unlike [`mean.age`][mean.age-ergmTerm] , this statistic is well-defined on
#'   an empty network. However, if used as a target, it appears to
#'   produce highly biased dissolution parameter estimates if the goal
#'   is to get an intended average duration.
#'
#' @usage
#' # binary: edge.ages
#'
#' @template ergmTerm-general
#'
#' @concept durational
InitErgmTerm.edge.ages<-function(nw, arglist, ...) {
  stopifnot_dynamic(nw, ...)
  a <- check.ErgmTerm(nw, arglist)
  
  list(name="edge_ages",
       coef.names = "edge.ages",
       duration=TRUE,
       dependence=FALSE,
       auxiliaries = ~.lasttoggle)
}

#' @name mean.age-ergmTerm
#' @title Average age of an extant tie
#' @description Average age of an extant tie
#' @details This term adds one statistic equaling the average, over all ties
#'   present in the network, of the amount of time elapsed since
#'   formation.
#'
#' @usage
#' # binary: mean.age(emptyval=0, log=FALSE)
#' @param emptyval can be used to specify the value returned if the network is empty. This is, technically, an arbitrary value, but it should
#'   not have a substantial effect unless a non-negligible fraction of
#'   networks at the parameter configuration of interest is empty.
#' @param log logical specifying if mean log age should be returned instead of mean age
#'
#' @template ergmTerm-general
#'
#' @concept durational
InitErgmTerm.mean.age<-function(nw, arglist, ...) {
  stopifnot_dynamic(nw, ...)
  a <- check.ErgmTerm(nw, arglist,
                      varnames = c("emptyval","log"),
                      vartypes = c("numeric","logical"),
                      defaultvalues = list(0,FALSE),
                      required = c(FALSE,FALSE))
  
  
  list(name="mean_age",
       coef.names = if(a$log) "mean.log.age" else "mean.age",
       inputs = c(a$emptyval,if(a$log) 1 else 0),
       emptynwstats = a$emptyval,
       duration=TRUE,
       dependence=FALSE,
       auxiliaries = ~.lasttoggle)
}

#' @name nodefactor.mean.age-ergmTerm
#' @title Average ages of extant half-ties incident on nodes of specified attribute levels
#' @description Average ages of extant half-ties incident on nodes of specified attribute levels
#' @details This term adds one statistic for each level of `attr` ,
#'   equaling the average, over all half-ties incident on nodes of that level,
#'   of the amount of time elapsed since formation.
#'
#' @usage
#' # binary: nodefactor.mean.age(attr, levels=NULL, emptyval=0, log=FALSE)
#' @template ergmTerm-attr
#' @param levels controls what levels are included. Note that the default
#'   `levels` value for `nodefactor.mean.age` retains all levels, unlike the default
#'   for `nodefactor` , which omits the first level.
#' @param emptyval can be used to specify the value returned if the network is empty. A different value may be
#'   specified for each level of `attr`. The length of `emptyval` should either be 1 (in which case that value
#'   is used for every level of `attr` ) or should be equal to the number of retained levels of `attr` , in
#'   which case the `i` th value in `emptyval` is used for the `i` th retained level of `attr`. This is,
#'   technically, an arbitrary value, but it should
#'   not have a substantial effect unless a non-negligible fraction of
#'   networks at the parameter configuration of interest is empty.
#' @param log logical specifying if mean log age should be returned instead of mean age
#'
#' @template ergmTerm-general
#'
#' @concept durational
InitErgmTerm.nodefactor.mean.age <- function(nw, arglist, ...) {
  stopifnot_dynamic(nw, ...)
  a <- check.ErgmTerm(nw, arglist,
                      varnames = c("attr", "levels", "emptyval", "log"),
                      vartypes = c(ERGM_VATTR_SPEC, ERGM_LEVELS_SPEC, "numeric", "logical"),
                      defaultvalues = list(NULL, NULL, 0, FALSE),
                      required = c(TRUE, FALSE, FALSE, FALSE))
  attrarg <- a$attr                        
  levels <- a$levels  

  nodecov <- ergm_get_vattr(attrarg, nw)
  attrname <- attr(nodecov, "name")
  u <- ergm_attr_levels(levels, nodecov, nw, levels = sort(unique(nodecov)))  
  
  if(!is.logical(a$log) || length(a$log) > 1) {
    ergm_Init_abort("log argument must be a length 1 logical")
  }  

  cn <- paste(if(a$log) "nodefactor.mean.log.age" else "nodefactor.mean.age", paste(attrname,collapse="."), u, sep=".")
  
  if(!is.numeric(a$emptyval)) {
    ergm_Init_abort("emptyval must be numeric")
  }
  
  if(length(a$emptyval) == 1) {
    a$emptyval <- rep(a$emptyval, length(cn))
  } else if(length(a$emptyval) != length(cn)) {
    ergm_Init_abort("emptyval must be either length 1 or have length equal to the number of statistics")
  }

  #   Recode to numeric
  nodecov <- c(0, match(nodecov, u, nomatch = 0)) - 1

  ### Construct the list to return
  list(name="nodefactor_mean_age",
       coef.names = cn,
       dependence=FALSE,
       duration=TRUE,
       auxiliaries = ~.lasttoggle,
       inputs = NULL, # passed by name
       emptynwstats = as.double(a$emptyval),
       nodecov = as.integer(nodecov),
       log = as.integer(a$log))  
}

#' @name nodemix.mean.age-ergmTerm
#' @title Average ages of extant ties of specified mixing types
#' @description Average ages of extant ties of specified mixing types
#' @details This term adds one statistic for each mixing type of `attr` ,
#'   equaling the average, over all ties of that mixing type,
#'   of the amount of time elapsed since formation.
#'
#' @usage
#' # binary: nodemix.mean.age(attr, b1levels=NULL, b2levels=NULL, levels=NULL,
#' #                          levels2=NULL, emptyval=0, log=FALSE)
#' @template ergmTerm-attr
#' @param b1levels,b2levels,levels,level2 control what statistics are included in the model and the order in which they appear. `levels2` apply to all networks; `levels` applies to unipartite networks; `b1levels` and `b2levels` apply to bipartite networks (see Specifying Vertex attributes and Levels (`?nodal_attributes`) for details)
#' @param emptyval can be used to specify the value returned if the network is empty. A different value may be
#'   specified for each mixing type of `attr`. The length of `emptyval` should either be 1 (in which case that value
#'   is used for every mixing type of `attr` ) or should be equal to the number of retained mixing types of `attr` , in
#'   which case the `i` th value in `emptyval` is used for the `i` th retained mixing type of `attr`. This is,
#'   technically, an arbitrary value, but it should
#'   not have a substantial effect unless a non-negligible fraction of
#'   networks at the parameter configuration of interest is empty.
#' @param log logical specifying if mean log age should be returned instead of mean age
#'
#' @template ergmTerm-general
#'
#' @concept durational
InitErgmTerm.nodemix.mean.age <- function(nw, arglist, ...) {
  stopifnot_dynamic(nw, ...)
  a <- check.ErgmTerm(nw, arglist,
                      varnames = c("attr", "b1levels", "b2levels", "levels", "levels2", "emptyval", "log"),
                      vartypes = c(ERGM_VATTR_SPEC, ERGM_LEVELS_SPEC, ERGM_LEVELS_SPEC, ERGM_LEVELS_SPEC, ERGM_LEVELS_SPEC, "numeric", "logical"),
                      defaultvalues = list(NULL, NULL, NULL, NULL, NULL, 0, FALSE),
                      required = c(TRUE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE))
                      
  attrarg <- a$attr
  b1levels <- a$b1levels
  b2levels <- a$b2levels  

  nodecov <- ergm_get_vattr(attrarg, nw)
  attrname <- attr(nodecov, "name")
  
  if(!is.logical(a$log) || length(a$log) > 1) {
    ergm_Init_abort("log argument must be a length 1 logical")
  }  
  
  if (is.bipartite(nw)) {
    #  So undirected network storage but directed mixing
    
    b1nodecov <- ergm_get_vattr(attrarg, nw, bip = "b1")
    b2nodecov <- ergm_get_vattr(attrarg, nw, bip = "b2")
    
    b1namescov <- ergm_attr_levels(b1levels, b1nodecov, nw, sort(unique(b1nodecov)))
    b2namescov <- ergm_attr_levels(b2levels, b2nodecov, nw, sort(unique(b2nodecov)))
    
    nr <- length(b1namescov)
    nc <- length(b2namescov)
    
    levels2.list <- transpose(expand.grid(row = b1namescov, col = b2namescov, stringsAsFactors=FALSE))
    indices2.grid <- expand.grid(row = 1:nr, col = nr + 1:nc)
   
    levels2.sel <- ergm_attr_levels(a$levels2, list(row = b1nodecov, col = b2nodecov), nw, levels2.list)
    
    rows2keep <- match(levels2.sel,levels2.list, NA)
    rows2keep <- rows2keep[!is.na(rows2keep)]
  
    u <- indices2.grid[rows2keep,]
  
    # Recode to numeric
    b1nodecov <- match(b1nodecov,b1namescov,nomatch=length(b1namescov)+1)
    b2nodecov <- match(b2nodecov,b2namescov,nomatch=length(b2namescov)+1)
  
    namescov <- c(b1namescov, b2namescov)
    
    nodecov <- c(b1nodecov, b2nodecov)
    
    cn <- paste(if(a$log) "nodemix.mean.log.age" else "nodemix.mean.age", paste(attrname,collapse="."), apply(matrix(namescov[as.matrix(u)],ncol=2),
                                       1,paste,collapse="."), sep=".")
                                       
    ## the +1 for nrow and ncol are needed by the way we code non-included b1 and b2 levels above
    indmat <- matrix(0L, nrow = nr + 1, ncol = nc + 1)
    u[,2L] <- u[,2L] - nr
    indmat[as.matrix(u)] <- seq_len(NROW(u))
    indmat <- indmat - 1L
  } else { # So one mode, but could be directed or undirected
    u <- ergm_attr_levels(a$levels, nodecov, nw, sort(unique(nodecov)))
    namescov <- u 
    
    nr <- length(u)
    nc <- length(u)

    levels2.list <- transpose(expand.grid(row = u, col = u, stringsAsFactors=FALSE))
    indices2.grid <- expand.grid(row = 1:nr, col = 1:nc)
    uun <- as.vector(outer(u,u,paste,sep="."))
    
    if (!is.directed(nw)) {
        rowleqcol <- indices2.grid$row <= indices2.grid$col
        levels2.list <- levels2.list[rowleqcol]
        indices2.grid <- indices2.grid[rowleqcol,]
        uun <- uun[rowleqcol]
    }    
   
    levels2.sel <- ergm_attr_levels(a$levels2, list(row = nodecov, col = nodecov), nw, levels2.list)
    
    rows2keep <- match(levels2.sel,levels2.list, NA)
    rows2keep <- rows2keep[!is.na(rows2keep)]
  
    u <- indices2.grid[rows2keep,]
    uun <- uun[rows2keep]

    nodecov <- match(nodecov,namescov,nomatch=length(namescov)+1)
    
    cn <- paste(if(a$log) "nodemix.mean.log.age" else "nodemix.mean.age", paste(attrname,collapse="."), uun, sep=".")

    ## the +1 for nrow and ncol are needed by the way we code non-included b1 and b2 levels above
    indmat <- matrix(0L, nrow = nr + 1, ncol = nc + 1)
    indmat[as.matrix(u)] <- seq_len(NROW(u))
    if(!is.directed(nw)) indmat <- indmat + t(indmat) - diag(diag(indmat))
    indmat <- indmat - 1L
  }
  
  if(!is.numeric(a$emptyval)) {
    ergm_Init_abort("emptyval must be numeric")
  }
  
  if(length(a$emptyval) == 1) {
    a$emptyval <- rep(a$emptyval, length(cn))
  } else if(length(a$emptyval) != length(cn)) {
    ergm_Init_abort("emptyval must be either length 1 or have length equal to the number of statistics")
  }
    
  ### Construct the list to return
  list(name = "nodemix_mean_age", 
       coef.names = cn, # will need to incorporate log in name
       dependence = FALSE,
       emptynwstats = as.double(a$emptyval),
       duration = TRUE,
       auxiliaries = ~.lasttoggle,
       inputs = NULL, # passed by name below; we will also use emptynwstats in the changestat functions
       log = as.integer(a$log),
       nr = as.integer(nr + 1),
       nc = as.integer(nc + 1),
       indmat = as.integer(t(indmat)),
       nodecov = as.integer(c(0L, nodecov) - 1L)) # two shifts to make the C code cleaner
}

#' @name edgecov.mean.age-ergmTerm
#' @title Weighted average age of an extant tie
#' @description Weighted average age of an extant tie
#' @details This term adds one statistic equaling the average, over all ties
#'   present in the network, of the amount of time elapsed since
#'   formation, weighted by a (nonnegative) dyadic covariate.
#'
#'   The behavior when there are negative weights is undefined.
#'
#' @usage
#' # binary: edgecov.mean.age(x, attrname=NULL, emptyval=0)
#' @template ergmTerm-x-attrname
#' @param emptyval can be used to specify the value returned if the network is empty (or all extant edges have been weighted 0). This is, technically, an arbitrary value, but it should
#'   not have a substantial effect unless a non-negligible fraction of
#'   networks at the parameter configuration of interest is empty
#'   and/or if only a few dyads have nonzero weights.
#'
#' @template ergmTerm-general
#'
#' @concept durational
InitErgmTerm.edgecov.mean.age<-function(nw, arglist, ...) {
  stopifnot_dynamic(nw, ...)
  a <- check.ErgmTerm(nw, arglist, 
                      varnames = c("x", "attrname", "emptyval", "log"),
                      vartypes = c("matrix,network,character", "character", "numeric", "logical"),
                      defaultvalues = list(NULL, NULL, 0, FALSE),
                      required = c(TRUE, FALSE, FALSE, FALSE))
  
  ### Check the network and arguments to make sure they are appropriate.
  ### Process the arguments
  if(is.network(a$x))
    xm<-as.matrix(a$x,matrix.type="adjacency",a$attrname)
  else if(is.character(a$x))
    xm<-get.network.attribute(nw,a$x)
  else
    xm<-as.matrix(a$x)
  
  ### Construct the list to return
  if(!is.null(a$attrname)) {
    # Note: the sys.call business grabs the name of the x object from the 
    # user's call.  Not elegant, but it works as long as the user doesn't
    # pass anything complicated.
    cn<-paste("mean",if(a$log) ".log" else "", ".age", as.character(a$attrname), sep = ".")
  } else {
    cn<-paste("mean",if(a$log) ".log" else "", ".age", as.character(sys.call(0)[[3]][2]), sep = ".")
  }
  inputs <- c(a$emptyval, if(a$log) 1 else 0, NCOL(xm), as.double(xm))
  attr(inputs, "ParamsBeforeCov") <- 3
  list(name="edgecov_mean_age", coef.names = cn, inputs = inputs, duration=TRUE, dependence=FALSE, emptynwstats = a$emptyval, auxiliaries = ~.lasttoggle)
}

#' @name edgecov.ages-ergmTerm
#' @title Weighted sum of ages of extant ties
#' @description Weighted sum of ages of extant ties
#' @details This term adds one statistic equaling sum, over all ties
#'   present in the network, of the amount of time elapsed since
#'   formation, multiplied by a dyadic covariate.
#'
#'   "Weights" can be negative.
#'
#'   Unlike [`edgecov.mean.age`][edgecov.mean.age-ergmTerm] , this statistic is well-defined on
#'   an empty network. However, if used as a target, it appears to
#'   produce highly biased dissolution parameter estimates if the goal
#'   is to get an intended average duration.
#'
#' @usage
#' # binary: edgecov.ages(x, attrname=NULL)
#' @template ergmTerm-x-attrname
#'
#' @template ergmTerm-general
#'
#' @concept durational
InitErgmTerm.edgecov.ages<-function(nw, arglist, ...) {
  stopifnot_dynamic(nw, ...)
  a <- check.ErgmTerm(nw, arglist, 
                      varnames = c("x", "attrname"),
                      vartypes = c("matrix,network,character", "character"),
                      defaultvalues = list(NULL, NULL),
                      required = c(TRUE, FALSE))
  
  ### Check the network and arguments to make sure they are appropriate.
  ### Process the arguments
  if(is.network(a$x))
    xm<-as.matrix(a$x,matrix.type="adjacency",a$attrname)
  else if(is.character(a$x))
    xm<-get.network.attribute(nw,a$x)
  else
    xm<-as.matrix(a$x)
  
  ### Construct the list to return
  if(!is.null(a$attrname)) {
    # Note: the sys.call business grabs the name of the x object from the 
    # user's call.  Not elegant, but it works as long as the user doesn't
    # pass anything complicated.
    cn<-paste("edgecov.ages", as.character(a$attrname), sep = ".")
  } else {
    cn<-paste("edgecov.ages", as.character(sys.call(0)[[3]][2]), sep = ".")
  }
  inputs <- c(NCOL(xm), as.double(xm))
  attr(inputs, "ParamsBeforeCov") <- 1
  list(name="edgecov_ages", coef.names = cn, inputs = inputs, duration=TRUE, dependence=FALSE, auxiliaries = ~.lasttoggle)
}

################################################################################

#' @name degree.mean.age-ergmTerm
#' @title Average age of ties incident on nodes having a given degree
#' @description Average age of ties incident on nodes having a given degree
#' @details This term adds one
#'   network statistic to the model for each element in `d` ; the \eqn{i} th
#'   such statistic equals the average, among all ties incident on nodes
#'   with degree exactly `d[i]` , of the amount of time elapsed
#'   since the tie's formation. The optional argument
#'   `byarg` specifies a vertex attribute (see
#'   Specifying Vertex Attributes and Levels
#'   for details). If specified, then separate degree
#'   statistics are calculated for nodes having each separate
#'   value of the attribute.
#'
#' @usage
#' # binary: degree.mean.age(d, byarg=NULL, emptyval=0)
#' @param d a vector of distinct integers
#' @template ergmTerm-byattr
#' @template ergmTerm-emptyval
#'
#' @template ergmTerm-general
#'
#' @concept durational
InitErgmTerm.degree.mean.age<-function(nw, arglist, ...) {
  stopifnot_dynamic(nw, ...)
  a <- check.ErgmTerm(nw, arglist, directed=FALSE,
                      varnames = c("d", "byarg", "emptyval", "log"),
                      vartypes = c("numeric", ERGM_VATTR_SPEC, "numeric", "logical"),
                      defaultvalues = list(NULL, NULL, 0, FALSE),
                      required = c(TRUE, FALSE, FALSE, FALSE))
  d<-a$d; byarg <- a$byarg
  if(!is.null(byarg)) {
    nodecov <- ergm_get_vattr(byarg, nw)
    u<-sort(unique(nodecov))
    if(any(is.na(nodecov))){u<-c(u,NA)}
    nodecov <- match(nodecov,u) # Recode to numeric
    if (length(u)==1)
         stop ("Attribute given to degree.mean.age() has only one value", call.=FALSE)
  }
  if(!is.null(byarg)) {
    # Combine degree and u into 2xk matrix, where k=length(d)*length(u)
    lu <- length(u)
    du <- rbind(rep(d,lu), rep(1:lu, rep(length(d), lu)))
    if (any(du[1,]==0)) {
      stop("Age of ties incident on nodes with degree 0 is meaningless.")
    }
  }
  
  if(is.null(byarg)) {
    if(length(d)==0){return(NULL)}
    coef.names <- paste("degree",d,".mean",if(a$log) ".log" else "", ".age",sep="")
    name <- "degree_mean_age"
    inputs <- c(d)
  } else {
    if(ncol(du)==0) {return(NULL)}
    #  No covariates here, so "ParamsBeforeCov" unnecessary
    # See comment in d_degree_by_attr function
    coef.names <- paste("deg", du[1,], ".", byarg,u[du[2,]],".mean",if(a$log) ".log" else "", ".age", sep="")
    name <- "degree_by_attr_mean_age"
    inputs <- c(as.vector(du), nodecov)
  }
  list(name=name,coef.names=coef.names, inputs=c(a$emptyval, if(a$log) 1 else 0, inputs),
       emptynwstats=rep(a$emptyval,length(coef.names)), duration=TRUE, dependence=TRUE, auxiliaries = ~.lasttoggle)
}

################################################################################

#' @name degrange.mean.age-ergmTerm
#' @title Average age of ties incident on nodes having degree in a given range
#' @description Average age of ties incident on nodes having degree in a given range
#' @details This term adds one
#'   network statistic to the model for each element of `from` (or `to` ); the \eqn{i} th
#'   such statistic equals the average, among all ties incident on nodes
#'   with degree greater than or equal to
#'   `from[i]` but strictly less than `to[i]` , of the amount of time elapsed
#'   since the tie's formation. The optional argument
#'
#' @usage
#' # binary: degrange.mean.age(from, to=+Inf, byarg=NULL, emptyval=0)
#' @param from,to vectors of distinct
#'   integers or `+Inf` , for `to` . If one of the vectors has
#'   length 1, it is recycled to the length of the other. Otherwise, they
#'   must have the same length.
#' @template ergmTerm-byattr
#' @template ergmTerm-emptyval
#' @template ergmTerm-general
#'
#' @concept durational
InitErgmTerm.degrange.mean.age<-function(nw, arglist, ...) {
  stopifnot_dynamic(nw, ...)
  a <- check.ErgmTerm(nw, arglist, directed=FALSE,
                      varnames = c("from", "to", "byarg", "emptyval", "log"),
                      vartypes = c("numeric", "numeric", ERGM_VATTR_SPEC, "numeric", "logical"),
                      defaultvalues = list(NULL, Inf, NULL, 0, FALSE),
                      required = c(TRUE, FALSE, FALSE, FALSE, FALSE))
  from<-a$from; to<-a$to; byarg <- a$byarg
  to <- ifelse(to==Inf, network.size(nw)+1, to)

  if(length(to)==1 && length(from)>1) to <- rep(to, length(from))
  else if(length(from)==1 && length(to)>1) from <- rep(from, length(to))
  else if(length(from)!=length(to)) stop("The arguments of term edges.ageinterval must have arguments either of the same length, or one of them must have length 1.")
  else if(any(from>=to)) stop("Term degrange.mean.age must have from<to.")
  
  if(!is.null(byarg)) {
    nodecov <- ergm_get_vattr(byarg, nw)
    u<-sort(unique(nodecov))
    if(any(is.na(nodecov))){u<-c(u,NA)}
    nodecov <- match(nodecov,u) # Recode to numeric
    if (length(u)==1)
         stop ("Attribute given to degrange.mean.age() has only one value", call.=FALSE)
  }
  if(!is.null(byarg)) {
    # Combine degree and u into 3xk matrix, where k=length(from)*length(u)
    lu <- length(u)
    du <- rbind(rep(from,lu), rep(to,lu), rep(1:lu, rep(length(from), lu)))
    if (any(du[1,]==0)) {
      stop("Age of ties incident on nodes with degree 0 is meaningless.")
    }
  }
  
  if(is.null(byarg)) {
    if(length(from)==0){return(NULL)}
    coef.names <- ifelse(to>=network.size(nw)+1,
                         paste("deg",from,"+",".mean",if(a$log) ".log" else "", ".age",sep=""),
                         paste("deg",from,"to",to,".mean",if(a$log) ".log" else "", ".age",sep=""))
    name <- "degrange_mean_age"
    inputs <- c(rbind(from,to))
  } else {
    if(ncol(du)==0) {return(NULL)}
    #  No covariates here, so "ParamsBeforeCov" unnecessary
    # See comment in d_degrange_by_attr function
    coef.names <- ifelse(du[2,]==network.size(nw)+1,
                         paste("deg", du[1,], "+.", byarg, u[du[3,]],".mean",if(a$log) ".log" else "", ".age", sep=""),
                         paste("deg", du[1,], "to", du[2,], ".", byarg, u[du[3,]],".mean",if(a$log) ".log" else "", ".age", sep=""))
    name <- "degrange_by_attr_mean_age"
    inputs <- c(as.vector(du), nodecov)
  }
  list(name=name,coef.names=coef.names, inputs=c(a$emptyval, if(a$log) 1 else 0, inputs),
       emptynwstats=rep(a$emptyval,length(coef.names)), duration=TRUE, dependence=TRUE, auxiliaries = ~.lasttoggle)
}


