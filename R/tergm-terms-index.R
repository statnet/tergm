#  File R/tergm-terms-index.R in package tergm, part of the Statnet suite
#  of packages for network analysis, https://statnet.org .
#
#  This software is distributed under the GPL-3 license.  It is free,
#  open source, and has the attribution requirements (GPL Section 7) at
#  https://statnet.org/attribution
#
#  Copyright 2008-2021 Statnet Commons
#######################################################################

#' Temporally-Sensitive Operator and Durational Terms used in Exponential
#' Family Random Graph Models
#'
#' @name ergmTerm
#' @aliases ergm-terms ergm.terms terms-ergm terms.ergm nodemix.mean.age InitErgmTerm.edgecov
#' @docType package
#' @description The terms described here are unique to temporal networks: each summarizes
#' some type of change or durational information.
#'
#' The operator terms include: `Form()`, `Persist()`, `Diss()`,
#' `Cross()` and `Change()`.  These are used to specify how the
#' \code{\link[ergm:ergm-terms]{ergm-terms}} in a formula are evaluated across a network
#' time-series.  Note, you cannot use an operator within another operator, so
#' `Cross(~Form(~edges))` is not a valid specification.
#'
#' The durational terms are distinguished either by their name,
#' `mean.age`, or their name extensions: `<name>.ages`,
#' `<name>.mean.age`, and `<name>.age.interval`.  In contrast to
#' their named \code{\link[ergm:ergm-terms]{ergm-terms}} counterparts, these durational
#' terms take into account the elapsed time since each (term-relevant) dyad in
#' the network was last toggled.
#'
#' As currently implemented, the package does not support use of durational
#' terms during estimation.  But durational terms may be used as targets,
#' monitors, or summary statistics.  The ability to use these terms in the
#' estimation of models is under development.
#'
#' @section Terms to represent network statistics included in the [`tergm`][tergm-package] package:
#' ## Term index
#' \if{latex}{\Sexpr[results=rd,stage=render]{ergm:::.formatIndexLatex(ergm:::.buildTermsDataframe("ergmTerm"))}}
#' \if{text}{\Sexpr[results=rd,stage=render]{ergm:::.formatIndexText(ergm:::.buildTermsDataframe("ergmTerm"))}}
#' \if{html}{\Sexpr[results=rd,stage=render]{ergm:::.formatIndexHtml(ergm:::.buildTermsDataframe("ergmTerm"))}}
#'
#' ## Frequently-used terms
#' \if{latex}{\Sexpr[results=rd,stage=render]{ergm:::.formatMatrixLatex(ergm:::.termMatrix("ergmTerm", categories=ergm:::FREQUENTLY_USED_TERM_CATEGORIES, only.include='frequently-used'))}}
#' \if{text}{\Sexpr[results=rd,stage=render]{ergm:::.formatMatrixText(ergm:::.termMatrix("ergmTerm", categories=ergm:::FREQUENTLY_USED_TERM_CATEGORIES, only.include='frequently-used'))}}
#' \if{html}{\Sexpr[results=rd,stage=render]{ergm:::.formatMatrixHtml(ergm:::.termMatrix("ergmTerm", categories=ergm:::FREQUENTLY_USED_TERM_CATEGORIES, only.include='frequently-used'))}}
#'
#' ## Operator terms
#' \if{latex}{\Sexpr[results=rd,stage=render]{ergm:::.formatMatrixLatex(ergm:::.termMatrix("ergmTerm", categories=ergm:::OPERATOR_CATEGORIES, only.include='operator'))}}
#' \if{text}{\Sexpr[results=rd,stage=render]{ergm:::.formatMatrixText(ergm:::.termMatrix("ergmTerm", categories=ergm:::OPERATOR_CATEGORIES, only.include='operator'))}}
#' \if{html}{\Sexpr[results=rd,stage=render]{ergm:::.formatMatrixHtml(ergm:::.termMatrix("ergmTerm", categories=ergm:::OPERATOR_CATEGORIES, only.include='operator'))}}
#'
#' ## All terms
#' \if{latex}{\Sexpr[results=rd,stage=render]{ergm:::.formatMatrixLatex(ergm:::.termMatrix("ergmTerm"))}}
#' \if{text}{\Sexpr[results=rd,stage=render]{ergm:::.formatMatrixText(ergm:::.termMatrix("ergmTerm"))}}
#' \if{html}{\Sexpr[results=rd,stage=render]{ergm:::.formatMatrixHtml(ergm:::.termMatrix("ergmTerm"))}}
#'
#' ## Terms by concept
#' \if{latex}{\Sexpr[results=rd,stage=render]{ergm:::.formatTocLatex(ergm:::.termToc("ergmTerm"))}}
#' \if{text}{\Sexpr[results=rd,stage=render]{ergm:::.formatTocText(ergm:::.termToc("ergmTerm"))}}
#' \if{html}{\Sexpr[results=rd,stage=render]{ergm:::.formatTocHtml(ergm:::.termToc("ergmTerm"))}}
#'
#' @seealso \code{\link[ergm:ergm-terms]{ergm-terms}} (from the
#' \code{\link[=ergm-package]{ergm}} package), \code{\link{ergm}},
#' \code{\link{network}}, \code{\link{%v%}}, \code{\link{%n%}}
#'
#' @references
#' - Krivitsky, P.N. (2012). Modeling of Dynamic Networks based on
#' Egocentric Data with Durational Information. *Pennsylvania State
#' University Department of Statistics Technical Report*, 2012(2012-01).
#' \url{https://web.archive.org/web/20170830053722/https://stat.psu.edu/research/technical-report-files/2012-technical-reports/TR1201A.pdf}
#'
#' - Krivitsky, P.N. (2012). Modeling Tie Duration in ERGM-Based Dynamic
#' Network Models. *Pennsylvania State University Department of Statistics
#' Technical Report*, 2012(2012-02).
#'
#' @keywords models
#'
NULL

#' Constraints and hints for Temporal Exponential Family Random Graph Models
#'
#' @name ergmConstraint
#' @aliases ergm-constraints constraints-ergm ergm.constraints constraints.ergm ergm-hints hints-ergm ergm.hints hints.ergm
#' @docType package
#' @description This page describes the hints and network sample
#'   space constraints that are included with the [`tergm`][tergm-package]
#'   package. For more information, and instructions for using hints and constraints, see
#'   [`ergmConstraint`][ergm::ergmConstraint] and [`ergm`].
#'
#' @section Constraints and hints to represent network statistics included in the [`tergm`][tergm-package] package:
#' ## Term index
#' \if{latex}{\Sexpr[results=rd,stage=render]{ergm:::.formatIndexLatex(ergm:::.buildTermsDataframe("ergmConstraint"))}}
#' \if{text}{\Sexpr[results=rd,stage=render]{ergm:::.formatIndexText(ergm:::.buildTermsDataframe("ergmConstraint"))}}
#' \if{html}{\Sexpr[results=rd,stage=render]{ergm:::.formatIndexHtml(ergm:::.buildTermsDataframe("ergmConstraint"))}}
#'
#' @references
#' - Krivitsky PN and Handcock MS (2014) A Separable Model for Dynamic Networks.
#'   *Journal of the Royal Statistical Society, Series B*, 76(1): 29-46. \doi{10.1111/rssb.12014}
#'
#' - Goodreau SM, Handcock MS, Hunter DR, Butts CT, Morris M (2008a).  A \pkg{statnet} Tutorial.
#'   *Journal of Statistical Software*, 24(8).  \url{https://www.jstatsoft.org/v24/i08/}.
#'
#' - Hunter, D. R. and Handcock, M. S. (2006) Inference in curved exponential family models for networks,
#'   *Journal of Computational and Graphical Statistics*.
#'
#' - Hunter DR, Handcock MS, Butts CT, Goodreau SM, Morris M (2008b).
#'   \pkg{ergm}: A Package to Fit, Simulate and Diagnose Exponential-Family Models for Networks.
#'   *Journal of Statistical Software*, 24(3). \url{https://www.jstatsoft.org/v24/i03/}.
#'
#' - Krivitsky PN (2012). Exponential-Family Random Graph Models for Valued Networks.
#'   *Electronic Journal of Statistics*, 2012, 6, 1100-1128. \doi{10.1214/12-EJS696}
#'
#' - Morris M, Handcock MS, Hunter DR (2008).  Specification of Exponential-Family Random Graph Models:
#'   Terms and Computational Aspects. *Journal of Statistical Software*, 24(4).
#'   \url{https://www.jstatsoft.org/v24/i04/}.
#'
#' @keywords models
#'
NULL
