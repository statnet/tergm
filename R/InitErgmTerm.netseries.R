#  File R/InitErgmTerm.netseries.R in package tergm, part of the Statnet suite
#  of packages for network analysis, https://statnet.org .
#
#  This software is distributed under the GPL-3 license.  It is free,
#  open source, and has the attribution requirements (GPL Section 7) at
#  https://statnet.org/attribution
#
#  Copyright 2008-2020 Statnet Commons
#######################################################################

.same_constraints <- function(nwl, nattr){
  map(nwl, get.network.attribute, nattr) %>% map(NVL, ~.) %>% map(empty_env) %>% all_identical
}

join_nets <- function(nwl, blockID, blockName){
  if(!.same_constraints(nwl, "constraints")) stop("Networks have differing constraint structures. This is not supported at this time.")
  if(!.same_constraints(nwl, "obs.constraints")) stop("Networks have differing observation processes. This is not supported at this time.")

  nw <- combine_networks(nwl, blockID.vattr=blockID, blockName.vattr=blockName, ignore.nattr = c(eval(formals(combine_networks)$ignore.nattr), "constraints", "obs.constraints", "ergm"), subnet.cache=TRUE)

  nw %n% "ergm" <- .combine_ergmlhs(nwl)

  nw %ergmlhs% "constraints" <-
      if(NVL(nwl[[1]] %ergmlhs% "constraints",base_env(~.))==base_env(~.))
        base_env(substitute(~blockdiag(blockID), list(blockID=blockID)))
      else
        append_rhs.formula(nwl[[1]]%ergmlhs%"constraints", list(call("blockdiag",".NetworkID")), TRUE)
  if(!is.null(nwl[[1]]%ergmlhs%"obs.constraints")) nw %ergmlhs% "obs.constraints" <- nwl[[1]]%ergmlhs%"obs.constraints"

  nw
}

#' A network series specification for modelling for conditional modeling.
#'
#' A function for specifying the LHS of a temporal network series ERGM.
#'
#' @param ... series specification, in one of three formats:
#' 
#'   1. A list of identically- dimensioned and directed networks.
#'
#'   1. Several networks as arguments.
#'
#'   1. A [`networkDynamic`] object and a numeric vector of time indices.
#'
#' @param order how many previous networks to store as an accessible
#'   covariate of the model.
#'
#' @param NA.impute How missing dyads in transitioned-from networks
#'   are be imputed when using conditional estimation. See argument
#'   \code{imputers} of \code{\link{impute.network.list}} for details.
#'
#' @return A network object with temporal metadata.
#'
#' @seealso [Help on model specification][ergm-terms] for specific terms.
#' 
#' @examples
#'
#' data(samplk)
#'
#' # Method 1: list of networks
#' monks <- NetSeries(list(samplk1,samplk2,samplk3))
#' ergm(monks ~ Form(~edges)+Diss(~edges))
#'
#' # Method 2: networks as arguments
#' monks <- NetSeries(samplk1,samplk2,samplk3)
#' ergm(monks ~ Form(~edges)+Diss(~edges))
#'
#' # Method 3: networkDynamic and time points:
#' ## TODO
#'
#' @export
NetSeries <- function(..., order=1, NA.impute=NULL){
  if(order>1) stop("Higher-order network models are not supported at this time.")
  
  args <- list(...)
  #' @importFrom methods is
  if(is(args[[1]], "networkDynamic")){
    nw <- args[[1]]
    times <- if(length(args)>=2) args[[2]]
             else{warning("the times argument was not provided to specify sampling time points for networkDynamic object. Modeling transition from time 0 to 1."); c(0,1)}
    # Grab only the needed vertices.
    subnw <- network.extract(nw, onset=min(times), terminus=max(times)+min(abs(diff(times)))/2)
    # Grab the vector of vertex activity indicators for each time
    # point, bind them into an n*T matrix. If any rows have
    # variability, we have a changing composition.
    if(any(apply(sapply(times, function(t) networkDynamic::is.active(subnw, at=t, v=1:network.size(subnw))), 1, var)>0)) warning("Active vertex set varies from time point to time point. Estimation may not work.")
    
    nwl <- lapply(times, function(t) networkDynamic::network.collapse(subnw, at=t, retain.all.vertices=TRUE))
  }else if(all(sapply(args, is, "network"))){
    nwl <- args
    times <- seq_along(nwl)
  }else if(is.list(args[[1]]) && all(sapply(args[[1]], is, "network"))){
    nwl <- args[[1]]
    times <- seq_along(nwl)
  }else stop("Unrecognized format for network series specification. See help for information.")

  nwl0 <- impute.network.list(nwl[-length(nwl)], NA.impute, nwl.append=nwl[length(nwl)])

  nwl0.NA <- sapply(nwl0, network.naedgecount)>0 # Update which networks have missing dyads.
  
  if(any(nwl0.NA)) stop("Transitioned-from network(s) at time(s) ", paste.and(times[-length(times)][nwl0.NA]), " has missing dyads that cannot be imputed using the selected imputation options. Fix or add imputation options via CMLE.NA.impute control parameter.")

  nwl <- lapply(seq_along(nwl)[-1], function(t){
    nwl[[t]] %n% ".PrevNets" <- list(nwl0[[t-1]])
    nwl[[t]]
  })

  # nwl now has all networks in the series but the first, which each network's previous network attached as a network attribute. In other words, it's a list of *transitions*.

  if(!all_identical(sapply(nwl,network.size)) || !all_identical(sapply(nwl,is.directed)) || !all_identical(sapply(nwl, `%n%`, "bipartite"))){
    stop("Networks in the network series passed must all be of the same size, directedness, and bipartite status.")
  }
  
  # Now, just combine them using the Networks() constructor.
  nw <- join_nets(nwl,".NetworkID",".NetworkName")
  # Add previous networks combined.
  .PrevNet <- join_nets(nwl0,".NetworkID",".NetworkName")
  nw %ergmlhs% "constraints" <- nonsimp_update.formula(nw%ergmlhs%"constraints", .~.+discord(.PrevNet), from.new=".PrevNet")

  nw
}


InitErgmTerm.Form <- function(nw, arglist,  ..., env=baseenv()) {
  a <- check.ErgmTerm(nw, arglist,
                      varnames = c("formula"),
                      vartypes = c("formula"),
                      defaultvalues = list(NULL),
                      required = c(TRUE))

  if(!is(nw, "combined_networks")) {
    return(InitErgmTerm.FormE(nw = nw, arglist = a["formula"], ...))
  }

  f <- a$formula
  ult(f) <- call("Form1",as.formula(call("~",ult(f)), env=env)) # e.g., a~b -> a~Form1(~b)
  environment(f) <- environment(a$formula)
  a$formula <- f

  # Just call N() operator.
  term <- as.call(c(list(as.name("Cross")),a))
  out <- call.ErgmTerm(term, env=env, nw=nw, ...)
  out
  # TODO: Fix coefficient names.
}

#' @importFrom utils modifyList
# One formation transition
InitErgmTerm.Form1 <- function(nw, arglist,  ...){
  a <- check.ErgmTerm(nw, arglist,
                      varnames = c("formula"),
                      vartypes = c("formula"),
                      defaultvalues = list(NULL),
                      required = c(TRUE))

  m <- ergm_model(a$formula, nw,...)
  
  c(list(name="on_union_net_Network", pkgname="ergm",
         auxiliaries = ~.union.net((nw%n%".PrevNets")[[1]], implementation="Network"),
         submodel = m),
    modifyList(wrap.ergm_model(m, nw, ergm_mk_std_op_namewrap("Form")),
               list(emptynwstats=summary(m, (nw%n%".PrevNets")[[1]])))
    )
}

InitErgmTerm.Diss <- function(nw, arglist,  ..., env=baseenv()) {
  a <- check.ErgmTerm(nw, arglist,
                      varnames = c("formula"),
                      vartypes = c("formula"),
                      defaultvalues = list(NULL),
                      required = c(TRUE))

  if(!is(nw, "combined_networks")) {
    return(InitErgmTerm.DissE(nw = nw, arglist = a["formula"], ...))
  }

  f <- a$formula
  ult(f) <- call("Diss1",as.formula(call("~",ult(f)), env=env)) # e.g., a~b -> a~Form1(~b)
  environment(f) <- environment(a$formula)
  a$formula <- f

  # Just call N() operator.
  term <- as.call(c(list(as.name("Cross")),a))
  out <- call.ErgmTerm(term, env=env, nw=nw, ...)
  out
  # TODO: Fix coefficient names.
}

# One dissolution transition
InitErgmTerm.Diss1 <- function(nw, arglist,  ...){
  a <- check.ErgmTerm(nw, arglist,
                      varnames = c("formula"),
                      vartypes = c("formula"),
                      defaultvalues = list(NULL),
                      required = c(TRUE))

  m <- ergm_model(a$formula, nw,...)
  
  c(list(name="on_intersect_net_Network", pkgname="ergm",
         auxiliaries = ~.intersect.net((nw%n%".PrevNets")[[1]], implementation="Network"),
         submodel = m),
    wrap.ergm_model(m, nw, ergm_mk_std_op_namewrap("Diss")))
}

# Auxiliary to extract crossectional networks.
InitErgmTerm..crossnets <- function(nw, arglist, ...){
  a <- check.ErgmTerm(nw, arglist,
                      varnames = c("attrname"),
                      vartypes = c("character"),
                      defaultvalues = list(NULL),
                      required = c(TRUE))

  list(name="_crossnets", coef.names=c(), iinputs=c(unlist(.block_vertexmap(nw, a$attrname))), dependence=FALSE)
}

# A term for a cross-sectional model for the network in a series.
#' @import purrr
InitErgmTerm.Cross <- function(nw, arglist,...){
  a <- check.ErgmTerm(nw, arglist,
                      varnames = c("formula"),
                      vartypes = c("formula"),
                      defaultvalues = list(NULL),
                      required = c(TRUE))

  auxiliaries <- base_env(~.crossnets(".NetworkID"))

  nwl <- .split_constr_network(nw, ".NetworkID", ".NetworkName")
  nwnames <- names(nwl)
  nn <- length(nwl)

  ms <- lapply(nwl, function(nw1){
    ergm_model(a$formula, nw1, ...)
  })

  nparams <- ms %>% map_int(nparam, canonical=FALSE)
  nstats <-  ms %>% map_int(nparam, canonical=TRUE)

  if(!all_identical(nparams) && !all_identical(ms %>% map(param_names, canonical=FALSE))) ergm_Init_abort("The transitions appear to have different sets of parameters. This may be a consequence of a change in network composition. The levels= argument may help.")
  nparam <- nparams[1]

  wms <- mapply(wrap.ergm_model, ms, nwl, SIMPLIFY=FALSE, USE.NAMES=FALSE)

  gss <- map(wms, "emptynwstats")
  gs <-
    if(all(map_lgl(gss, is.null))){ # Linear combination of 0s is 0.
      NULL
    }else{ # All numeric or NULL
      gss %>% map(NVL, 0) %>% Reduce(`+`, .)
    }

  c(list(name="OnCrossNets", submodels=ms, emptynwstats = gs, auxiliaries = auxiliaries),
    wms[[1]][c("coef.names","map","gradient","params","dependence","offset")])
}
