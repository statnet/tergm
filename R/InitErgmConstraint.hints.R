#  File R/InitErgmConstraint.hints.R in package tergm, part of the Statnet
#  suite of packages for network analysis, https://statnet.org .
#
#  This software is distributed under the GPL-3 license.  It is free, open
#  source, and has the attribution requirements (GPL Section 7) at
#  https://statnet.org/attribution .
#
#  Copyright 2008-2025 Statnet Commons
################################################################################

#' @templateVar name discord
#' @title Discordant dyads
#' @description Propose toggling discordant dyads with greater
#'   frequency (typically about 50 percent).  May be used in
#'   dynamic fitting and simulation.
#'
#' @usage
#' # discord
#'
#' @template ergmHint-general
#' @concept dyad-independent
InitErgmConstraint.discord <- function(nw, arglist, ...) {
  a <- check.ErgmTerm(nw, arglist,
                      varnames = "ref",
                      vartypes = "character",
                      defaultvalues = list(NULL),
                      required = c(FALSE))

  nw <- if(is.character(a$ref)) nw %n% a$ref else nw

  list(dependence = FALSE, priority = 10, nw = nw)
}
