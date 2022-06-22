#  File R/InitErgmTerm.EGMME.R in package tergm, part of the
#  Statnet suite of packages for network analysis, https://statnet.org .
#
#  This software is distributed under the GPL-3 license.  It is free,
#  open source, and has the attribution requirements (GPL Section 7) at
#  https://statnet.org/attribution .
#
#  Copyright 2008-2022 Statnet Commons
################################################################################

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

#' @templateVar name Form
#' @title The Formation Operator Term
#' @description The Formation Operator Term
#' @details This term accepts a model formula
#'   and produces the corresponding model for the post-formation network:
#'   effectively a network containing both previous time step's ties and ties just formed,
#'   the union of the previous and current network. This is the equivalent of the
#'   old-style `formation` model.
#'
#' @usage
#' # binary: Form(formula)
#' @template ergmTerm-formula
#'
#' @template ergmTerm-general
#'
#' @concept operator
#' @concept durational
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
    modifyList(wrap.ergm_model(m, nw, ergm_mk_std_op_namewrap("Form")),
               list(emptynwstats=NULL)))
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

#' @templateVar name Persist
#' @title The Persistence Operator Term
#' @description The Persistence Operator Term
#' @details This term accepts a model formula
#'   and produces the corresponding model for the
#'   post-dissolution/persistence network: effectively the network containing
#'   ties that persisted since the last time step.
#'
#'   This is the equivalent of the old-style `dissolution` model. So
#'   a larger positive coefficient for `Persist()` operator means
#'   less dissolution. It
#'   produces the same results as the new `Diss()` operator, except the
#'   signs of the coefficients are negated.
#'
#' @usage
#' # binary: Persist(formula)
#' @template ergmTerm-formula
#'
#' @template ergmTerm-general
#'
#' @concept operator
#' @concept durational
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
    modifyList(wrap.ergm_model(m, nw, ergm_mk_std_op_namewrap("Persist")),
               list(emptynwstats=NULL)))
}

#' @templateVar name Diss
#' @title The Dissolution Operator Term
#' @description The Dissolution Operator Term
#' @details This term accepts a model formula
#'   and produces the corresponding model for the post-dissolution
#'   network (same as `Persist()` ), but with all statistics negated.
#'
#'   Note: This is not the equivalent of the old style `dissolution` model,
#'   because the signs of the coefficients are reversed. So a larger positive
#'   coefficient for `Diss()` operator means more dissolution.
#'
#' @usage
#' # binary: Diss(formula)
#' @template ergmTerm-formula
#'
#' @template ergmTerm-general
#'
#' @concept operator
#' @concept durational
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
    modifyList(wrap.ergm_model(m, nw, ergm_mk_std_op_namewrap("Diss")),
               list(emptynwstats=NULL)))
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

#' @templateVar name Change
#' @title The Change Operator Term
#' @description The Change Operator Term
#' @details This term accepts a model formula
#'   and produces the corresponding model for a network constructed
#'   by taking the dyads that have changed between time steps.
#'
#' @usage
#' # binary: Change(formula)
#' @template ergmTerm-formula
#'
#' @template ergmTerm-general
#'
#' @concept operator
#' @concept durational
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
    modifyList(wrap.ergm_model(m, nw, ergm_mk_std_op_namewrap("Change")),
               list(emptynwstats=NULL)))
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
#' @description The \code{.SubsetStatistics} term is used internally in \code{tergm}'s EGMME
#' to remove offset statistics from the targets model.  At this time it is not intended to
#' be used by the end user, but is documented here for completeness.  The behavior described
#' below may change without warning in the future, so do not rely on this term in your own code!
#'
#' The \code{.SubsetStatistics} term takes two arguments: \code{formula} and \code{statistics}.  The
#' \code{formula} argument is a model formula and the \code{statistics} argument is a
#' \code{levels}-type argument for selecting statistics to retain in the model generated by
#' \code{formula} (and the input network).  If the original model is linear, so is the derived
#' model (meaning thetas are dropped and/or rearranged just as the etas are), and if the original
#' model is curved, so is the derived model, which retains the full set of original thetas, even
#' if some of them do not influence any retained etas.
#' @noRd
InitErgmTerm..SubsetStatistics <- function(nw, arglist, ...) {
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

  wm <- wrap.ergm_model(m, nw)

  if(is.curved(m)) {
    wm <- replace(wm, c("coef.names", "emptynwstats", "map", "gradient"),
                  list(coef.names = param_names(m, canonical = TRUE)[statistics],
                       emptynwstats = wm$emptynwstats[statistics],
                       map = function(x, n, ...) { ergm.eta(x, m$etamap)[statistics] },
                       gradient = function(x, n, ...) { ergm.etagrad(x, m$etamap)[,statistics,drop=FALSE] }))
  } else {
    for(name in c("minpar", "maxpar", "coef.names", "offset", "emptynwstats")) {
      if(!is.null(wm[[name]])) wm[[name]] <- wm[[name]][statistics]
    }
  }

  c(list(name = "subset_stats",
         submodel = m,
         duration = is.durational(m),
         iinputs = statistics - 1L),
    ergm_propagate_ext.encode(m),
    wm)
}
