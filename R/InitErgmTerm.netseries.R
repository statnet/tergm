#  File R/InitErgmTerm.netseries.R in package tergm, part of the
#  Statnet suite of packages for network analysis, https://statnet.org .
#
#  This software is distributed under the GPL-3 license.  It is free,
#  open source, and has the attribution requirements (GPL Section 7) at
#  https://statnet.org/attribution .
#
#  Copyright 2008-2024 Statnet Commons
################################################################################

.same_constraints <- function(nwl, nattr){
  map(nwl, get.network.attribute, nattr) %>% map(NVL, ~.) %>% map(empty_env) %>% all_identical
}

join_nets <- function(nwl, blockID, blockName){
  if(!.same_constraints(nwl, "constraints")) stop("Networks have differing constraint structures. This is not supported at this time.")
  if(!.same_constraints(nwl, "obs.constraints")) stop("Networks have differing observation processes. This is not supported at this time.")

  nw <- combine_networks(nwl, blockID.vattr=blockID, blockName.vattr=blockName, ignore.nattr = c(eval(formals(combine_networks)$ignore.nattr), "constraints", "obs.constraints", "ergm"), subnet.cache=TRUE)

  nw %n% "ergm" <- combine_ergmlhs(nwl)

  nw %ergmlhs% "constraints" <-
      if(NVL(nwl[[1]] %ergmlhs% "constraints",base_env(~.))==base_env(~.))
        base_env(substitute(~blockdiag(blockID), list(blockID=blockID)))
      else
        append_rhs.formula(nwl[[1]]%ergmlhs%"constraints", list(call("blockdiag",".TimeID")), TRUE)
  if(!is.null(nwl[[1]]%ergmlhs%"obs.constraints")) nw %ergmlhs% "obs.constraints" <- nwl[[1]]%ergmlhs%"obs.constraints"

  nw
}

#' A network series specification for conditional modeling.
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
#'   \code{imputers} of [impute.network.list()] for details.
#'
#' @return A network object with temporal metadata.
#'
#' @note It is not recommended to modify the network returned by
#'   `NetSeries` except by adding and removing edges, and even that
#'   must be done with some care, to avoid putting it into an
#'   inconsistent state.
#'
#'   It is almost always better to modify the original networks and
#'   regenerate the series.
#'
#' @seealso [`ergmTerm`] for specific terms.
#' 
#' @examples
#'
#' data(samplk)
#'
#' # Method 1: list of networks
#' monks <- NetSeries(list(samplk1,samplk2,samplk3))
#' ergm(monks ~ Form(~edges)+Diss(~edges))
#' ergm(monks ~ Form(~edges)+Persist(~edges))
#'
#' # Method 2: networks as arguments
#' monks <- NetSeries(samplk1,samplk2,samplk3)
#' ergm(monks ~ Form(~edges)+Diss(~edges))
#' ergm(monks ~ Form(~edges)+Persist(~edges))
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
    nwl[[t]] %n% ".TimeID" <- t
    nwl[[t]] %n% ".Time" <- times[t]
    nwl[[t]] %n% ".TimeDelta" <- times[t] - times[t-1]
    nwl[[t]]
  })

  # nwl now has all networks in the series but the first, which each network's previous network attached as a network attribute. In other words, it's a list of *transitions*.

  if(!all_identical(sapply(nwl,network.size)) || !all_identical(sapply(nwl,is.directed)) || !all_identical(sapply(nwl, `%n%`, "bipartite"))){
    stop("Networks in the network series passed must all be of the same size, directedness, and bipartite status.")
  }
  
  # Now, just combine them using the Networks() constructor.
  nw <- join_nets(nwl,".TimeID",".Time")
  # Add previous networks combined.
  nw %n% ".PrevNet" <- join_nets(nwl0,".TimeID",".Time")
  nw %ergmlhs% "constraints" <- nonsimp_update.formula(nw%ergmlhs%"constraints", .~.+discord(".PrevNet"))

  structure(nw, class = c("tergm_NetSeries", class(nw)))
}

.call_N <- function(term, nw, arglist, ..., env=baseenv()){
  a <- check.ErgmTerm(nw, arglist,
                      varnames = c("formula", "lm", "subset", "weights", "contrasts", "offset", "label"),
                      vartypes = c("formula", "formula", "formula,logical,numeric,expression,call", "formula,logical,numeric,expression,call", "list", "formula,logical,numeric,expression,call", "character"),
                      defaultvalues = list(NULL, ~1, TRUE, 1, NULL, NULL, NULL),
                      required = c(TRUE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE))

  f <- a$formula
  ult(f) <- call(paste0(term,"1"), as.formula(call("~",ult(f)), env=env)) # e.g., a~b -> a~term1(~b)
  environment(f) <- environment(a$formula)
  a$formula <- f

  # Just call N() operator.
  call.ErgmTerm(as.call(c(list(as.name("N")),
                    c(a[c("formula", "lm", "subset", "weights", "contrasts", "offset")],
                      label=ergm_mk_std_op_namewrap(term,a$label),
                      .NetworkID=".TimeID", .NetworkName=".Time")
                    )),
                env=env, nw=nw, ...)
}

#' @templateVar name Form
#' @template ergmTerm-rdname
#' @usage NULL
#' @template ergmTerm-N-arguments
InitErgmTerm.Form <- function(nw, arglist,  ...){
  if(!is(nw, "tergm_NetSeries")) `InitErgmTerm.Form (dynamic)`(nw = nw, arglist = arglist, ...)
  else .call_N("Form", nw, arglist, ...)
}

#' @importFrom utils modifyList
# One formation transition
InitErgmTerm.Form1 <- function(nw, arglist,  ...){
  a <- check.ErgmTerm(nw, arglist,
                      varnames = c("formula"),
                      vartypes = c("formula"),
                      defaultvalues = list(NULL),
                      required = c(TRUE))

  m <- ergm_model(a$formula, nw, ..., offset.decorate=FALSE)
  ergm_no_ext.encode(m)

  c(list(name="on_union_net_Network", pkgname="ergm",
         auxiliaries = ~.union.net((nw%n%".PrevNets")[[1]], implementation="Network"),
         submodel = m),
    modifyList(wrap.ergm_model(m, nw, identity),
               list(emptynwstats=summary(m, (nw%n%".PrevNets")[[1]])))
    )
}

#' @templateVar name Persist
#' @template ergmTerm-rdname
#' @usage NULL
#' @template ergmTerm-N-arguments
InitErgmTerm.Persist <- function(nw, arglist,  ...) {
  if(!is(nw, "tergm_NetSeries")) `InitErgmTerm.Persist (dynamic)`(nw = nw, arglist = arglist,...)
  else .call_N("Persist", nw, arglist, ...)
}

# One dissolution transition
InitErgmTerm.Persist1 <- function(nw, arglist,  ...){
  a <- check.ErgmTerm(nw, arglist,
                      varnames = c("formula"),
                      vartypes = c("formula"),
                      defaultvalues = list(NULL),
                      required = c(TRUE))

  m <- ergm_model(a$formula, nw, ..., offset.decorate=FALSE)
  ergm_no_ext.encode(m)

  c(list(name="on_intersect_net_Network", pkgname="ergm",
         auxiliaries = ~.intersect.net((nw%n%".PrevNets")[[1]], implementation="Network"),
         submodel = m),
    wrap.ergm_model(m, nw, identity))
}

#' @templateVar name Diss
#' @template ergmTerm-rdname
#' @usage NULL
#' @template ergmTerm-N-arguments
InitErgmTerm.Diss <- function(nw, arglist,  ..., env=baseenv()) {
  if(!is(nw, "tergm_NetSeries")) `InitErgmTerm.Diss (dynamic)`(nw = nw, arglist = arglist, ...)
  else{
    .call_N("Diss", nw, arglist, ...)
  }
}

#
InitErgmTerm.Diss1 <- function(nw, arglist,  ..., env){
  a <- check.ErgmTerm(nw, arglist,
                      varnames = c("formula"),
                      vartypes = c("formula"),
                      defaultvalues = list(NULL),
                      required = c(TRUE))
  formula <- a$formula
  formula[2:3] <- c(-1, formula[[2]])
  term <- call("Persist1", substitute(~Sum(formula, I), list(formula=formula)))
  ergm_model(term_list(list(term), env=env), nw, ..., env=env, terms.only=TRUE)
}


#' @templateVar name Change
#' @template ergmTerm-rdname
#' @usage NULL
#' @template ergmTerm-N-arguments
InitErgmTerm.Change <- function(nw, arglist,  ...) {
  if(!is(nw, "tergm_NetSeries")) `InitErgmTerm.Change (dynamic)`(nw = nw, arglist = arglist, ...)
  else .call_N("Change", nw, arglist, ...)
}

# One difference transition
InitErgmTerm.Change1 <- function(nw, arglist,  ...){
  a <- check.ErgmTerm(nw, arglist,
                      varnames = c("formula"),
                      vartypes = c("formula"),
                      defaultvalues = list(NULL),
                      required = c(TRUE))

  m <- ergm_model(a$formula, nw, ..., offset.decorate=FALSE)
  ergm_no_ext.encode(m)

  c(list(name="on_discord_net_Network", pkgname="ergm",
         auxiliaries = ~.discord.net((nw%n%".PrevNets")[[1]], implementation="Network"),
         submodel = m),
    modifyList(wrap.ergm_model(m, nw, identity),
               list(emptynwstats=summary(m, (nw%n%".PrevNets")[[1]])))
    )
}

# A term for a cross-sectional model for the network in a series.

#' @templateVar name Cross
#' @title The Crossection Operator Term
#' @description The Crossection Operator Term
#' @details This term accepts a model formula
#'   and produces the corresponding model for the cross-sectional
#'   network. It is mainly useful for CMLE estimation, and has no effect (i.e.,
#'   `Cross(~TERM) == ~TERM` ) for EGMME and dynamic simulation.
#'
#' @usage
#' # binary: Cross(
#' #           formula,
#' #           lm = ~1,
#' #           subset = TRUE,
#' #           weights = 1,
#' #           contrasts = NULL,
#' #           offset = 0,
#' #           label = NULL
#' #         )
#' @template ergmTerm-formula
#' @template ergmTerm-N-arguments
#'
#' @template ergmTerm-general
#' @import purrr
#' @rawNamespace import(ergm.multi, except=c("snctrl"))
#'
#' @concept operator
#' @concept durational
InitErgmTerm.Cross <- function(nw, arglist, ..., env=baseenv()) {
  if (!is(nw, "tergm_NetSeries")) `InitErgmTerm.Cross (dynamic)`(nw = nw, arglist = arglist, ...)
  else .call_N("Cross", nw, arglist, ..., env=env)
}

# Pass through
InitErgmTerm.Cross1 <- function(nw, arglist,  ...){
  a <- check.ErgmTerm(nw, arglist,
                      varnames = c("formula"),
                      vartypes = c("formula"),
                      defaultvalues = list(NULL),
                      required = c(TRUE))

  ergm_model(a$formula, nw, ..., terms.only=TRUE)
}
