#  File R/InitErgmTerm.netseries.R in package tergm, part of the Statnet suite
#  of packages for network analysis, https://statnet.org .
#
#  This software is distributed under the GPL-3 license.  It is free, open
#  source, and has the attribution requirements (GPL Section 7) at
#  https://statnet.org/attribution .
#
#  Copyright 2008-2025 Statnet Commons
################################################################################

.call_N <- function(term, nw, arglist, ..., env=baseenv(),
                    ownargs = list(varnames = "formula", vartypes = "formula", defaultvalues = list(NULL), required = TRUE)) {
  a <- check.ErgmTerm(nw, arglist,
                      varnames = c(ownargs$varnames, "lm", "subset", "weights", "contrasts", "offset", "label"),
                      vartypes = c(ownargs$vartypes, "formula", "formula,logical,numeric,expression,call", "formula,logical,numeric,expression,call", "list", "formula,logical,numeric,expression,call", "character"),
                      defaultvalues = c(ownargs$defaultvalues, list(~1, TRUE, 1, NULL, NULL, NULL)),
                      required = c(ownargs$required, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE))

  f <- a$formula
  # e.g., a~b -> a~term1(~b)
  ult(f) <- as.call(c(list(paste0(term,"1"),
                           as.formula(call("~",ult(f)), env=env)),
                      a[setdiff(ownargs$varnames, "formula")])) 
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
#' @concept temporal
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
#' @concept temporal
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
#' @concept temporal
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
#' @concept temporal
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


#' @templateVar name Lag
#' @title An edge covariate on the previous network's statistics
#' @description This term takes an ERGM formula that is evaluated on a
#'   previous time step's network; the change statistics (optionally
#'   transformed) are used as an edge covariate for the current
#'   network. That is, a statistic of the form \eqn{\sum_{i,j}
#'   y^{t}_{i,j} \Delta_{i,j} g(y^{t-1})} or a transformed variant
#'   \eqn{\sum_{i,j} y^{t}_{i,j} f(\Delta_{i,j} g(y^{t-1}),
#'   y^{t-1}_{i,j})}, where \eqn{Delta_{i,j} g(y) = g(y+(i,j)) -
#'   g(y-(i,j))}.
#'
#' @usage
#' # binary: Lag(
#' #           formula,
#' #           lag = 1,
#' #           transform = NULL,
#' #           lm = ~1,
#' #           subset = TRUE,
#' #           weights = 1,
#' #           contrasts = NULL,
#' #           offset = 0,
#' #           label = NULL
#' #         )
#' @template ergmTerm-formula
#' @param lag how many time steps to look back; at this time, only 1
#'   is implemented.
#' @param transform a `function(x, y)` given an array of covariates
#'   and the sociomatrix of the *previous* network to optionally
#'   transform `x`, `NULL` to leave unchanged, or a character string
#'   indicating a preset; current presets include \describe{
#'
#' \item{`"signed"`}{multiplies each change statistic by \eqn{(2
#' y^{t-1}_{i,j} - 1)}}
#'
#' \item{`"nonzero"`, "`positive`", "`negative`"}{indicator of whether
#' the change statistic has the respective value}
#' 
#' }
#' @template ergmTerm-N-arguments
#'
#' @template ergmTerm-general
#'
#' @references
#'
#' Almquist, Z. W. and Butts, C. T. (2014). Logistic Network
#' Regression for Scalable Analysis of Networks with Joint Edge/Vertex
#' Dynamics. *Sociological Methodology*, 44(1),
#' 273-321. \doi{10.1177/0081175013520159}
#'
#' @concept operator
#' @concept temporal
#' @concept dyad-independent
InitErgmTerm.Lag <- function(nw, arglist, ..., env=baseenv()) {
  if (!is(nw, "tergm_NetSeries")) ergm_Init_stop("This term does not support non-netseries input at this time.")
  ownargs <- list(varnames = c("formula", "lag", "transform"),
                  vartypes = c("formula", "numeric", "function,character"),
                  defaultvalues = list(NULL, 1, function(x, ...) x),
                  required = c(TRUE, FALSE, FALSE))
  .call_N("Lag", nw, arglist, ownargs = ownargs, ..., env=env)
}

InitErgmTerm.Lag1 <- function(nw, arglist, ...) {
  a <- check.ErgmTerm(nw, arglist,
                      varnames = c("formula", "lag", "transform"),
                      vartypes = c("formula", "numeric", "function,formula,character"),
                      defaultvalues = list(NULL, 1, NULL),
                      required = c(TRUE, FALSE, FALSE))

  NVL(a$transform) <- function(x, ...) x

  if (is.character(a$transform)) {
    presets <- c("signed", "nonzero", "positive")
    a$transform <- switch(match.arg(a$transform, presets),
                          signed = function(x, y, ...) c(2*y-1)*x,
                          nonzero = function(x, ...) x != 0,
                          positive = function(x, ...) x > 0,
                          negative = function(x, ...) x < 0,
                          ergm_Init_stop(sQuote(a$transform), " is not a valid preset; currently valid presets include ", paste.and(sQuote(presets)), "."))
  }

  if (a$lag != 1) ergm_Init_stop("Lags other than 1 are not currently implemented.")

  mple <- ergmMPLE(a$formula, basis = (nw%n%".PrevNets")[[1]],
                   expand.bipartite = TRUE, output = "array")
  y <- mple$response
  x <- mple$predictor

  pred <- a$transform(x = x, y = y)
  if (ncol(pred) == 0) ergm_Init_stop("Formula produced 0 statistics.")

  terms <- imap(dimnames(pred)[[3]],
                function(nm, i) term_list(call("edgecov", x = pred[, , i], attrname = nm), env = environment(a$formula))) |>
    do.call(c, args = _)

  m <- ergm_model(terms, nw, ..., terms.only = TRUE)
  param_names(m) <- map(list(param_names(m, canonical = FALSE), param_names(m, canonical = TRUE)),
                        substring, 9, 10000)
  m
}
