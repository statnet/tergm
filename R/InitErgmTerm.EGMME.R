#  File R/InitErgmTerm.EGMME.R in package tergm, part of the Statnet suite
#  of packages for network analysis, https://statnet.org .
#
#  This software is distributed under the GPL-3 license.  It is free,
#  open source, and has the attribution requirements (GPL Section 7) at
#  https://statnet.org/attribution
#
#  Copyright 2008-2020 Statnet Commons
#######################################################################

stopifnot_dynamic <- function(nw, ..., dynamic=FALSE, .netseries.OK=FALSE){
  if(!dynamic && ! "lasttoggle" %in% list.network.attributes(nw)){
    msg <- paste0("This term requires either dynamic data (",
                  if(.netseries.OK) "network series or ",
                  "network dynamic or last toggle information) ",
                  "or for dynamic mode to be set (typically by passing ",
                  sQuote("dynamic=TRUE)"), " to the top-level function.")
    ergm_Init_abort(msg)
  }
}

`InitErgmTerm.Form (dynamic)` <- function(nw, arglist,  ...) {
  stopifnot_dynamic(nw, .netseries.OK=TRUE, ...)
  a <- check.ErgmTerm(nw, arglist,
                      varnames = c("formula"),
                      vartypes = c("formula"),
                      defaultvalues = list(NULL),
                      required = c(TRUE))

  m <- ergm_model(a$formula, nw, ..., offset.decorate=FALSE)

  c(list(name="on_union_lt_net_Network",
         auxiliaries = ~.union.lt.net + .lasttoggle + .previous.lt.net,
         submodel = m,
         duration=TRUE),
    ergm_propagate_ext.encode(m),
    wrap.ergm_model(m, nw, ergm_mk_std_op_namewrap("Form")))
}

InitErgmTerm..union.lt.net<-function(nw, arglist, ...) {
  stopifnot_dynamic(nw, ...)
  a <- check.ErgmTerm(nw, arglist,
                      varnames = c(),
                      vartypes = c(),
                      defaultvalues = list(),
                      required = c())

  list(name="_union_lt_net_Network",
       coef.names=c(),
       auxiliaries = ~ .lasttoggle,
       duration=TRUE,
       dependence=FALSE)
}


`InitErgmTerm.Persist (dynamic)` <- function(nw, arglist,  ...) {
  stopifnot_dynamic(nw, .netseries.OK=TRUE, ...)
  a <- check.ErgmTerm(nw, arglist,
                      varnames = c("formula"),
                      vartypes = c("formula"),
                      defaultvalues = list(NULL),
                      required = c(TRUE))

  m <- ergm_model(a$formula, nw, ..., offset.decorate=FALSE)

  c(list(name="on_intersect_lt_net_Network",
         auxiliaries = ~.intersect.lt.net() + .lasttoggle,
         submodel = m,
         duration=TRUE),
    ergm_propagate_ext.encode(m),
    wrap.ergm_model(m, nw, ergm_mk_std_op_namewrap("Persist")))
}

`InitErgmTerm.Diss (dynamic)` <- function(nw, arglist,  ...) {
  stopifnot_dynamic(nw, .netseries.OK=TRUE, ...)
  a <- check.ErgmTerm(nw, arglist,
                      varnames = c("formula"),
                      vartypes = c("formula"),
                      defaultvalues = list(NULL),
                      required = c(TRUE))

  m <- ergm_model(a$formula, nw, ..., offset.decorate=FALSE)

  c(list(name="negate_on_intersect_lt_net_Network",
         auxiliaries = ~.intersect.lt.net() + .lasttoggle,
         submodel = m,
         duration=TRUE),
    ergm_propagate_ext.encode(m),
    wrap.ergm_model(m, nw, ergm_mk_std_op_namewrap("Diss")))
}

InitErgmTerm..intersect.lt.net<-function(nw, arglist, ...) {
  stopifnot_dynamic(nw, ...)
  a <- check.ErgmTerm(nw, arglist,
                      varnames = c(),
                      vartypes = c(),
                      defaultvalues = list(),
                      required = c())

  list(name="_intersect_lt_net_Network",
       coef.names=c(),
       auxiliaries = ~ .lasttoggle,
       duration=TRUE,
       dependence=FALSE)
}


`InitErgmTerm.Change (dynamic)` <- function(nw, arglist,  ...) {
  stopifnot_dynamic(nw, .netseries.OK=TRUE, ...)
  a <- check.ErgmTerm(nw, arglist,
                      varnames = c("formula"),
                      vartypes = c("formula"),
                      defaultvalues = list(NULL),
                      required = c(TRUE))

  m <- ergm_model(a$formula, nw, ..., offset.decorate=FALSE)

  c(list(name="on_discord_lt_net_Network",
         auxiliaries = ~.discord.lt.net() + .lasttoggle + .previous.lt.net,
         submodel = m,
         duration=TRUE),
    ergm_propagate_ext.encode(m),
    wrap.ergm_model(m, nw, ergm_mk_std_op_namewrap("Change")))
}


InitErgmTerm..discord.lt.net<-function(nw, arglist, ...) {
  stopifnot_dynamic(nw, ...)
  a <- check.ErgmTerm(nw, arglist,
                      varnames = c(),
                      vartypes = c(),
                      defaultvalues = list(),
                      required = c())

  list(name="_discord_lt_net_Network",
       coef.names=c(),
       auxiliaries = ~ .lasttoggle,
       duration=TRUE,
       dependence=FALSE)
}


InitErgmTerm..previous.lt.net<-function(nw, arglist, ...) {
  stopifnot_dynamic(nw, ...)
  a <- check.ErgmTerm(nw, arglist,
                      varnames = c(),
                      vartypes = c(),
                      defaultvalues = list(),
                      required = c())

  list(name="_previous_lt_net_Network",
       coef.names=c(),
       auxiliaries = ~ .lasttoggle,
       duration=TRUE,
       dependence=FALSE)
}

InitErgmTerm..lasttoggle <- function(nw, arglist, ...){
  stopifnot_dynamic(nw, ...)
  a <- check.ErgmTerm(nw, arglist,
                      varnames = c(),
                      vartypes = c(),
                      defaultvalues = list(),
                      required = c())

  list(name="_lasttoggle",
       coef.names=c(),
       duration=TRUE,
       dependence=FALSE,
       ext.encode = function(el, nw0) list(time=as.integer(NVL(nw0 %n% "time", 0)), lasttoggle=as.integer(nw0 %n% "lasttoggle")),
       ext.decode = function(ext.state, el, nw0){
         nw0 %n% "time" <- ext.state$time
         nw0 %n% "lasttoggle" <- matrix(ext.state$lasttoggle, ncol=3)
         list(el=el, nw0=nw0)
       }
       )
}

#' @title An Internal Operator Term for Subsetting Statistics
#'
#' @description The \code{.Statistics} term is used internally in \code{tergm}'s EGMME
#' to remove offset statistics from the targets model.  At this time it is not intended to
#' be used by the end user, but is documented here for completeness.  The behavior described
#' below may change without warning in the future, so do not rely on this term in your own code!
#'
#' The \code{.Statistics} term takes two arguments: \code{formula} and \code{statistics}.  The
#' \code{formula} argument is a model formula and the \code{statistics} argument is a
#' \code{levels}-type argument for selecting statistics to retain in the model generated by
#' \code{formula} (and the input network).  If the original model is linear, so is the derived
#' model (meaning thetas are dropped and/or rearranged just as the etas are), and if the original
#' model is curved, so is the derived model, which retains the full set of original thetas, even
#' if some of them do not influence any retained etas.
#' @noRd
InitErgmTerm..Statistics <- function(nw, arglist, ...) {
  a <- check.ErgmTerm(nw, arglist,
                      varnames = c("formula", "statistics"),
                      vartypes = c("formula", ERGM_LEVELS_SPEC),
                      defaultvalues = list(NULL, NULL),
                      required = c(TRUE, FALSE))

  m <- ergm_model(a$formula, nw, ..., offset.decorate=FALSE)

  ## extract levels and do a little safety check
  statistics <- ergm_attr_levels(a$statistics, seq_len(nparam(m, canonical = TRUE)), nw, seq_len(nparam(m, canonical = TRUE)))
  statistics <- as.integer(statistics[statistics >= 1 & statistics <= nparam(m, canonical = TRUE)])
  if(length(statistics) == 0) return(NULL)

  wem <- wrap.ergm_model(m, nw)

  if(is.curved(m)) {
    wem <- modifyList(wem,
                      list(coef.names = param_names(m, canonical = TRUE)[statistics],
                           map = function(x, n, ...) { ergm.eta(x, m$etamap)[statistics] },
                           gradient = function(x, n, ...) { ergm.etagrad(x, m$etamap)[,statistics,drop=FALSE] }))
  } else {
    for(name in c("minpar", "maxpar", "coef.names", "offset", "emptynwstats")) {
      if(!is.null(wem[[name]])) wem[[name]] <- wem[[name]][statistics]
    }
  }

  c(list(name = "_statistics_term",
         submodel = m,
         duration = is.durational(m),
         nstats = length(statistics),
         stats = statistics - 1L),
    ergm_propagate_ext.encode(m),
    wem)
}
