#  File R/InitErgmTerm.netseries.R in package tergm, part of the Statnet suite
#  of packages for network analysis, https://statnet.org .
#
#  This software is distributed under the GPL-3 license.  It is free, open
#  source, and has the attribution requirements (GPL Section 7) at
#  https://statnet.org/attribution .
#
#  Copyright 2008-2025 Statnet Commons
################################################################################

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
