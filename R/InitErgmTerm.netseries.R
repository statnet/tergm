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

#' A network series representation.
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
NetSeries <- function(...){
  args <- list(...)
  if(all(sapply(args, is, "network"))){
    nwl <- args
  }else if(is.list(args[[1]]) && all(sapply(args[[1]], is, "network"))){
    nwl <- args[[1]]
  }else stop("Unrecognized format for network series specification. See help for information.")

  nwl <- lapply(seq_along(nwl)[-1], function(t){
    nwl[[t]] %n% ".previous" <- list(nwl[[t-1]])
    nwl[[t]]
  })

  # nwl now has all networks in the series but the first, which each network's previous network attached as an ergmlhs attribute. In other words, it's a list of *transitions*.

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
  inputs <- to_ergm_Cdouble(m)
  
  gs <- summary(m, (nw%n%".previous")[[1]])
  
  c(list(name="formation",
         coef.names = paste0("Form",'(',param_names(m, canonical=TRUE),')'),
         auxiliaries = ~.union.net((nw%n%".previous")[[1]], implementation="Network"),
         inputs=inputs,
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
  inputs <- to_ergm_Cdouble(m)
  
  gs <- summary(m)
  
  c(list(name="dissolution",
         coef.names = paste0("Diss",'(',param_names(m, canonical=TRUE),')'),
         auxiliaries = ~.intersect.net((nw%n%".previous")[[1]], implementation="Network"),
         inputs=inputs,
         emptynwstats=gs,
         dependence=!is.dyad.independent(m)),
    passthrough.curved.ergm_model(m, function(x) paste0('Diss(',x,')')))
}
