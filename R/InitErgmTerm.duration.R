#  File R/InitErgmTerm.duration.R in package tergm, part of the Statnet suite
#  of packages for network analysis, https://statnet.org .
#
#  This software is distributed under the GPL-3 license.  It is free,
#  open source, and has the attribution requirements (GPL Section 7) at
#  https://statnet.org/attribution
#
#  Copyright 2008-2020 Statnet Commons
#######################################################################


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

InitErgmTerm.edge.ages<-function(nw, arglist, ...) {
  stopifnot_dynamic(nw, ...)
  a <- check.ErgmTerm(nw, arglist)
  
  list(name="edge_ages",
       coef.names = "edge.ages",
       duration=TRUE,
       dependence=FALSE,
       auxiliaries = ~.lasttoggle)
}

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


