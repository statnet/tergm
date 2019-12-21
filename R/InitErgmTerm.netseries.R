######################################################################
#
# stergm CMLE estimation implemented as ergm terms, using auxiliary storage to
# keep track of the formation and dissolution networks.
#
# usage: ergm(y ~ formation(formula_f) + dissolution(formula_d))
# formula_f is the formation model formula, 
# formula_d is the dissolution model formula
# Optional parameter: markov=TRUE for the standard markov conditional tergm model
#                     markov=FALSE for unconditional sampling 
#                       (changing panel 2 affects transitions 1->2 and 2->3)
#
######################################################################

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
  #' @importFrom ergm.multi Networks
  Networks(nwl)
}


InitErgmTerm.Form <- function(nw, arglist, response=NULL,  ...) {
  a <- check.ErgmTerm(nw, arglist,
                      varnames = c("formula","lm","subset","weights","contrasts","offset","label"),
                      vartypes = c("formula","formula","formula,logical,numeric,expression,call","formula,logical,numeric,expression,call","list","formula,logical,numeric,expression,call","character"),
                      defaultvalues = list(NULL,~1,TRUE,1,NULL,NULL,NULL),
                      required = c(TRUE,FALSE,FALSE,FALSE,FALSE,FALSE,FALSE))

  f <- a$formula
  ult(f) <- call("Form1",call("~",ult(f))) # e.g., a~b -> a~Form1(~b)
  a$formula <- f

  # Just call N() operator. (TODO: implement an API in ergm to standardise this and obviate the need to use :::.)
  out <- ergm.multi:::InitErgmTerm.N(nw, a, response=response, ...)
  out$pkgname <- "ergm.multi"
  out
  # TODO: Fix coefficient names.
}

# One formation transition
InitErgmTerm.Form1 <- function(nw, arglist, response=NULL,  ...){
  a <- check.ErgmTerm(nw, arglist,
                      varnames = c("formula"),
                      vartypes = c("formula"),
                      defaultvalues = list(NULL),
                      required = c(TRUE))
  
  
  # get the network and formula
  f <- a$formula

  if(length(f)==2) f <- nonsimp_update.formula(f, nw~.)
  else nw <- ergm.getnetwork(f)
  
  m <- ergm_model(f, nw, response=response,...)
  
  gs <- summary(m, (nw%n%".PrevNets")[[1]])
  
  c(list(name="on_union_net_Network", pkgname="ergm",
         coef.names = paste0("Form",'(',param_names(m, canonical=TRUE),')'),
         auxiliaries = ~.union.net((nw%n%".PrevNets")[[1]], implementation="Network"),
         submodel = m,
         emptynwstats=gs,
         dependence=!is.dyad.independent(m)),
    passthrough.curved.ergm_model(m, function(x) paste0('Form(',x,')')))
}

InitErgmTerm.Diss <- function(nw, arglist, response=NULL,  ...) {
  a <- check.ErgmTerm(nw, arglist,
                      varnames = c("formula","lm","subset","weights","contrasts","offset","label"),
                      vartypes = c("formula","formula","formula,logical,numeric,expression,call","formula,logical,numeric,expression,call","list","formula,logical,numeric,expression,call","character"),
                      defaultvalues = list(NULL,~1,TRUE,1,NULL,NULL,NULL),
                      required = c(TRUE,FALSE,FALSE,FALSE,FALSE,FALSE,FALSE))

  f <- a$formula
  ult(f) <- call("Diss1",call("~",ult(f))) # e.g., a~b -> a~Form1(~b)
  a$formula <- f

  # Just call N() operator. (TODO: implement an API in ergm to standardise this and obviate the need to use :::.)
  out <- ergm.multi:::InitErgmTerm.N(nw, a, response=response, ...)
  out$pkgname <- "ergm.multi"
  out
  # TODO: Fix coefficient names.
}

# One dissolution transition
InitErgmTerm.Diss1 <- function(nw, arglist, response=NULL,  ...){
  a <- check.ErgmTerm(nw, arglist,
                      varnames = c("formula"),
                      vartypes = c("formula"),
                      defaultvalues = list(NULL),
                      required = c(TRUE))
  
  
  # get the network and formula
  f <- a$formula

  if(length(f)==2) f <- nonsimp_update.formula(f, nw~.)
  else nw <- ergm.getnetwork(f)
  
  m <- ergm_model(f, nw, response=response,...)
  
  gs <- summary(m)
  
  c(list(name="on_intersect_net_Network", pkgname="ergm",
         coef.names = paste0("Diss",'(',param_names(m, canonical=TRUE),')'),
         auxiliaries = ~.intersect.net((nw%n%".PrevNets")[[1]], implementation="Network"),
         submodel = m,
         emptynwstats=gs,
         dependence=!is.dyad.independent(m)),
    passthrough.curved.ergm_model(m, function(x) paste0('Diss(',x,')')))
}
