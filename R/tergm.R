#  File R/tergm.R in package tergm, part of the Statnet suite of packages for
#  network analysis, https://statnet.org .
#
#  This software is distributed under the GPL-3 license.  It is free, open
#  source, and has the attribution requirements (GPL Section 7) at
#  https://statnet.org/attribution .
#
#  Copyright 2008-2025 Statnet Commons
################################################################################
################################################################################
# tergm --- fit Separable Temporal ERGMs.
################################################################################



#' Temporal Exponential-Family Random Graph Models
#' 
#' [tergm()] fits Temporal ERGMs' (TERGMs) and Separable Temporal ERGMs' (STERGMs)
#' Conditional MLE (CMLE) (Krivitsky and Handcock, 2010) and Equilibrium
#' Generalized Method of Moments Estimator (EGMME) (Krivitsky, 2009).
#' 
#' @param formula an ERGM formula.
#'
#' @template constraints
#'
#' @param estimate One of "EGMME" for Equilibrium Generalized Method of Moments
#' Estimation, based on a single network with some temporal information and
#' making an assumption that it is a product of a TERGM process running to its
#' stationary (equilibrium) distribution; "CMLE" for Conditional Maximum
#' Likelihood Estimation, modeling a transition between two networks, or
#' "CMPLE" for Conditional Maximum PseudoLikelihood Estimation, using MPLE
#' instead of MLE.  CMPLE is extremely inaccurate at this time.
#' 
#' @param times For CMLE and CMPLE estimation, times or indexes at
#'   which the networks whose transition is to be modeled are
#'   observed. This argument is mandatory if \code{nw} is a
#'   [`networkDynamic`] and defaults to \code{1:length(nw)} (all
#'   transitions) if \code{nw} is a [`network.list`] or a
#'   [`list`]. Ignored when estimating EGMME or if LHS is already a
#'   [`NetSeries`]. Note that at this time, the selected time points
#'   will be treated as temporally adjacent. Irregluarly spaced time
#'   series are not supported at this time.
#' 
#' @param offset.coef Numeric vector to specify offset parameters.
#' @param targets One-sided [ergm()]-style formula specifying
#' statistics whose moments are used for the EGMME. Unused for CMLE and CMPLE.
#' Targets is required for EGMME estimation. It may contain any valid ergm
#' terms.  Any offset terms are used only during the 
#' preliminary SAN run; they are removed automatically for the EGMME proper.
#' If \code{targets} is specified as a character
#' (one of \code{"formation"} and \code{"dissolution"}) then
#' the function [.extract.fd.formulae()] is used to determine the
#' corresponding formula; the user should be aware of its behavior and limitations.
#' @param SAN.offsets Offset coefficients (if any) to use during the SAN run.
#' @param target.stats A vector specifying the values of the \code{targets}
#' statistics that EGMME will try to match.  Defaults to the statistics of
#' \code{nw}. Unused for CMLE and CMPLE.
#' @param eval.loglik Whether or not to calculate the log-likelihood
#'   of a CMLE TERGM fit. See [ergm()] for details. Can be
#'   set globally via `option(tergm.eval.loglik=...)`, falling back to
#'   `getOption("ergm.eval.loglik")` if not set.
#' @param control A list of control parameters for algorithm tuning.
#' Constructed using [control.tergm()].
#' @template verbose
#' @param \dots Additional arguments, to be passed to lower-level functions.
#' @param basis optional network data overriding the left hand side of \code{formula}
#'
#' @return [tergm()] returns an object of class [`tergm`] that
#'   inherits from [`ergm`] and has the usual methods ([coef.ergm()],
#'   [summary.ergm()], [mcmc.diagnostics()], etc.) implemented for
#'   it. Note that [gof()] only works for the CMLE method.
#'
#' @seealso [`network`], [`networkDynamic`], and [NetSeries()] for the data structures,
#'   [ergm()] and [`ergmTerm`] for model specification,
#'   package vignette \code{browseVignettes(package='tergm')} for a
#'   short demonstration, the Statnet web site
#'   \url{https://statnet.org/workshop-tergm/} for a tutorial

#' @references
#'
#' Krackhardt, D and Handcock, MS (2006) Heider vs Simmel: Emergent
#' features in dynamic structures.  ICML Workshop on Statistical
#' Network Analysis. Springer, Berlin, Heidelberg, 2006.
#'
#' Hanneke S, Fu W, and Xing EP (2010). Discrete
#' Temporal Models of Social Networks. \emph{Electronic Journal of Statistics},
#' 2010, 4, 585-605.
#' \doi{10.1214/09-EJS548}
#'
#' Krivitsky P.N. and Handcock M.S. (2014) A Separable Model for Dynamic Networks. \emph{Journal of the Royal Statistical Society, Series B}, 76(1): 29-46. \doi{10.1111/rssb.12014}
#' 
#' Krivitsky, P.N. (2012). Modeling of Dynamic Networks based on
#' Egocentric Data with Durational Information. \emph{Pennsylvania State
#' University Department of Statistics Technical Report}, 2012(2012-01).
#' \url{https://arxiv.org/abs/2203.06866}
#' 
#'
#' @examples
#' \dontrun{
#' # EGMME Example
#' par(ask=FALSE)
#' n<-30
#' g0<-network.initialize(n,dir=FALSE)
#' 
#' #                     edges, degree(1), mean.age
#' target.stats<-c(      n*1/2,    n*0.6,        20)
#' 
#' dynfit<-tergm(g0 ~ Form(~edges + degree(1)) + Diss(~edges),
#'                targets = ~edges+degree(1)+mean.age,
#'                target.stats=target.stats, estimate="EGMME",
#'                control=control.tergm(SA.plot.progress=TRUE))
#' 
#' par(ask=TRUE)
#' mcmc.diagnostics(dynfit)
#' summary(dynfit)
#' }
#' \donttest{
#' # CMLE Example
#' data(samplk)
#'
#' # Fit a transition from Time 1 to Time 2
#' samplk12 <- tergm(list(samplk1, samplk2)~
#'                   Form(~edges+mutual+transitiveties+cyclicalties)+
#'                   Diss(~edges+mutual+transitiveties+cyclicalties),
#'                   estimate="CMLE")
#'
#' mcmc.diagnostics(samplk12)
#' summary(samplk12)
#'
#' samplk12.gof <- gof(samplk12)
#'
#' samplk12.gof
#'
#' plot(samplk12.gof)
#'
#' plot(samplk12.gof, plotlogodds=TRUE)
#'
#' # Fit a transition from Time 1 to Time 2 and from Time 2 to Time 3 jointly
#' samplk123 <- tergm(list(samplk1, samplk2, samplk3)~
#'                    Form(~edges+mutual+transitiveties+cyclicalties)+
#'                    Diss(~edges+mutual+transitiveties+cyclicalties),
#'                    estimate="CMLE")
#' 
#' mcmc.diagnostics(samplk123)
#' summary(samplk123)
#' }
#'
#' @import network
#' @import networkDynamic
#' @importFrom utils packageVersion
#' @export
tergm <- function(formula, constraints = ~., estimate, times=NULL, offset.coef=NULL,
                   targets=NULL, target.stats=NULL, SAN.offsets = NULL,
                   eval.loglik=NVL(getOption("tergm.eval.loglik"), getOption("ergm.eval.loglik")),
                   control=control.tergm(),
                   verbose=FALSE, ..., basis = eval_lhs.formula(formula)) {
  check.control.class("tergm", "tergm")

  tergm_call <- match.call(ergm)

  if(!is.null(control$seed))  set.seed(as.integer(control$seed))

  estimate <- match.arg(estimate,c("CMLE","CMPLE","EGMME"))
  
  if(!inherits(formula,"formula"))
    stop("Argument formula must be a formula.")
  
  out <- switch(estimate,
                CMLE=tergm.CMLE(formula=formula, times=times, constraints=constraints, estimate="MLE", offset.coef=offset.coef, target.stats=target.stats, eval.loglik=eval.loglik,control=control, verbose=verbose, ..., basis = basis),
                CMPLE=tergm.CMLE(formula=formula, times=times, constraints=constraints, estimate="MPLE", offset.coef=offset.coef, target.stats=target.stats, eval.loglik=eval.loglik,control=control, verbose=verbose, ..., basis = basis),
                EGMME=tergm.EGMME(formula, constraints,
                  offset.coef,
                  targets, target.stats, SAN.offsets, estimate, control, verbose, basis = basis)
                )

  out$call <- tergm_call
  out$formula <- formula
  out$estimate <- estimate
  out$estimate.desc <- switch(estimate,
                              CMPLE = ,
                              CMLE = sub("Maximum Pseudolikelihood", "Conditional Maximum Pseudolikelihood", sub("Maximum Likelihood", "Conditional Maximum Likelihood", out$estimate.desc)),
                              EGMME = switch(control$EGMME.main.method,
                                             `Gradient-Descent`="Gradient Descent Equilibrium Generalized Method of Moments Results"))
  out$control <- control
  out$constraints <- constraints
  out$tergm_version <- packageVersion("tergm")

  out
}

# Trivial functions taking care of the differences between EGMME and MLE-based methods.

#' @noRd
#' @export
gof.tergm_EGMME <- function (object, ...)
  stop("Goodness of fit for TERGM EGMME is not implemented at this time.")

#' @noRd
#' @importFrom ergm mcmc.diagnostics
#' @export
mcmc.diagnostics.tergm_EGMME <- function(object, ...) NextMethod("mcmc.diagnostics", object, esteq=FALSE, ...)
