#  File R/tergm.R in package tergm, part of the Statnet suite
#  of packages for network analysis, https://statnet.org .
#
#  This software is distributed under the GPL-3 license.  It is free,
#  open source, and has the attribution requirements (GPL Section 7) at
#  https://statnet.org/attribution
#
#  Copyright 2008-2019 Statnet Commons
#######################################################################
################################################################################
# tergm --- fit Separable Temporal ERGMs.
################################################################################



#' Temporal Exponential-Family Random Graph Models
#' 
#' \code{\link{tergm}} is used for finding Temporal ERGMs' (TERGMs) and Separable Temporal ERGMs' (TERGMs)
#' Conditional MLE (CMLE) (Krivitsky and Handcock, 2010) and Equilibrium
#' Generalized Method of Moments Estimator (EGMME) (Krivitsky, 2009).
#' 
#' \strong{Model Terms} See \code{\link{ergm}} and \code{\link{ergm-terms}} for
#' details. At this time, only linear ERGM terms are allowed.  \itemize{
#' \item For a brief demonstration, please see the tergm package vignette:
#' \code{browseVignettes(package='tergm')} \item A more detailed tutorial is
#' avalible on the statnet wiki:
#' \url{http://statnet.csde.washington.edu/workshops/SUNBELT/current/tergm/tergm_tutorial.pdf}
#' \item For more usage examples, see the wiki page at
#' \url{https://statnet.csde.washington.edu/trac/wiki/tergmUsage} }
#' 
#' @param formula an ERGM formula.
#' @param constraints A one-sided formula specifying one or more constraints on
#' the support of the distribution of the networks being modeled, using syntax
#' similar to the \code{formula} argument. Multiple constraints may be given,
#' separated by \dQuote{+} operators.  Together with the model terms in the
#' formula and the reference measure, the constraints define the distribution
#' of networks being modeled.
#' 
#' It is also possible to specify a proposal function directly by passing a
#' string with the function's name. In that case, arguments to the proposal
#' should be specified through the \code{prop.args} argument to
#' \code{\link{control.ergm}}.
#' 
#' The default is \code{~.}, for an unconstrained model.
#' 
#' See the [ERGM constraints][ergm-constraints] documentation for the
#' constraints implemented in the **[ergm][ergm-package]** package.
#' Other packages may add their own constraints.
#' 
#' Note that not all possible combinations of constraints are supported.
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
#'   observed. Default to \code{c(0,1)} if \code{nw} is a
#'   \code{\link[networkDynamic]{networkDynamic}} and to
#'   \code{1:length(nw)} (all transitions) if \code{nw} is a
#'   \code{\link{network.list}} or a \code{\link{list}}. Unused for
#'   EGMME. Note that at this time, the selected time points will be
#'   treated as temporally adjacent. Irregluarly spaced time series
#'   are not supported at this time.
#' 
#' @param offset.coef Numeric vector to specify offset parameters.
#' @param targets One-sided \code{\link{ergm}}-style formula specifying
#' statistics whose moments are used for the EGMME. Unused for CMLE and CMPLE.
#' Targets is required for EGMME estimation. It may contain any valid ergm
#' terms. Any offset terms are used only during the preliminary SAN run;
#' they are removed automatically for the EGMME proper.
#' @param SAN.offsets Offset coefficients (if any) to use during the SAN run.
#' @param target.stats A vector specifying the values of the \code{targets}
#' statistics that EGMME will try to match.  Defaults to the statistics of
#' \code{nw}. Unused for CMLE and CMPLE.
#' @param eval.loglik Whether or not to calculate the log-likelihood
#'   of a CMLE TERGM fit. See \code{\link{ergm}} for details. Can be
#'   set globally via `option(tergm.eval.loglik=...)`, falling back to
#'   `getOption("ergm.eval.loglik")` if not set.
#' @param control A list of control parameters for algorithm tuning.
#' Constructed using \code{\link{control.tergm}}.
#' @param verbose logical or integer; if TRUE or positive, the program will
#' print out progress information. Higher values result in more output.
#' @param \dots Additional arguments, to be passed to lower-level functions.
#' @return \code{\link{tergm}} returns an object of class \code{\link{tergm}}
#' that is a list consisting of the following elements:
#' @return \code{\link{tergm}} returns an object of class \code{\link{tergm}}
#' that is a list consisting of the following elements:
#' \item{formula}{The model formula.}
#' \item{coef}{The fitted model coefficients.}
#' \item{targets}{For EGMME, the targets formula.}
#' \item{target.stats}{For EGMME, the target statistics.}
#' \item{estimate}{The type of estimate.}
#' \item{opt.history}{For EGMME, a matrix containing the full trace of 
#' the optimization process: coefficients tried and target statistics simulated.}
#' \item{sample}{For EGMME, an \code{\link{mcmc}} object containing target
#' statistics sampled at the estimate.}
#' \item{covar}{For EGMME, the full estimated
#' variance-covariance matrix of the parameter estimates.}
#' \item{fit}{For CMLE and CMPLE, an \code{\link{ergm}} object from
#' fitting the model. For EGMME, stripped down \code{\link{ergm}}-like lists.}
#' \item{network}{For \code{estimate=="EGMME"}, the original network; for \code{estimate=="CMLE"}
#' or \code{estimate=="CMPLE"}, a \code{\link{network.list}} (a discrete series
#' of networks) to which the model was fit.}
#' \item{control}{The control parameters used to fit the model.}
#' 
#' See the method \code{\link{print.tergm}} for details on how an
#' \code{\link{tergm}} object is printed.  Note that the method
#' \code{\link{summary.tergm}} returns a summary of the relevant parts of the
#' \code{\link{tergm}} object in concise summary format.
#' @seealso ergm, network, \%v\%, \%n\%, \code{\link{ergm-terms}}
#' @references \itemize{
#'
#' \item Krivitsky P.N. and Handcock M.S. (2014) A Separable Model for Dynamic Networks. \emph{Journal of the Royal Statistical Society, Series B}, 76(1): 29-46. \doi{10.1111/rssb.12014}
#' 
#' \item Krivitsky, P.N. (2012). Modeling of Dynamic Networks based on
#' Egocentric Data with Durational Information. \emph{Pennsylvania State
#' University Department of Statistics Technical Report}, 2012(2012-01).
#' \url{http://stat.psu.edu/research/technical-report-files/2012-technical-reports/modeling-of-dynamic-networks-based-on-egocentric-data-with-durational-information}
#' 
#' }
#' @examples
#' 
#' # EGMME Example
#' par(ask=FALSE)
#' n<-30
#' g0<-network.initialize(n,dir=FALSE)
#' 
#' #                     edges, degree(1), mean.age
#' target.stats<-c(      n*1/2,    n*0.6,        20)
#' 
#' dynfit<-tergm(g0 ~ FormE(~edges + degree(1)) + DissE(~edges),
#'                targets = ~edges+degree(1)+mean.age,
#'                target.stats=target.stats, estimate="EGMME",
#'                control=control.tergm(SA.plot.progress=TRUE))
#' 
#' par(ask=TRUE)
#' mcmc.diagnostics(dynfit)
#' summary(dynfit)
#' 
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
#' # Fit a transition from Time 1 to Time 2 and from Time 2 to Time 3 jointly
#' samplk123 <- tergm(list(samplk1, samplk2, samplk3)~
#'                    Form(~edges+mutual+transitiveties+cyclicalties)+
#'                    Diss(~edges+mutual+transitiveties+cyclicalties),
#'                    estimate="CMLE")
#' 
#' mcmc.diagnostics(samplk123)
#' summary(samplk123)
#' 
#' @import network
#' @import networkDynamic
#' @export
tergm <- function(formula, constraints = ~., estimate, times=NULL, offset.coef=NULL,
                   targets=NULL, target.stats=NULL, SAN.offsets = NULL,
                   eval.loglik=NVL(getOption("tergm.eval.loglik"), getOption("ergm.eval.loglik")),
                   control=control.tergm(),
                   verbose=FALSE, ...) {
  check.control.class("tergm", "tergm")
  
  if(!is.null(control$seed))  set.seed(as.integer(control$seed))

  estimate <- match.arg(estimate,c("CMLE","CMPLE","EGMME"))
  
  if(!inherits(formula,"formula"))
    stop("Argument formula must be a formula.")
  
  out <- switch(estimate,
                CMLE=tergm.CMLE(formula=formula, times=times, constraints=constraints, estimate="MLE", offset.coef=offset.coef, target.stats=target.stats, eval.loglik=eval.loglik,control=control, verbose=verbose, ...),
                CMPLE=tergm.CMLE(formula=formula, times=times, constraints=constraints, estimate="MPLE", offset.coef=offset.coef, target.stats=target.stats, eval.loglik=eval.loglik,control=control, verbose=verbose, ...),
                EGMME=tergm.EGMME(formula, constraints,
                  offset.coef,
                  targets, target.stats, SAN.offsets, estimate, control, verbose)
                )
  
  out$estimate <- estimate
  out$control <- control
  out$constraints <- constraints
  
  class(out)<-c("tergm","ergm")
  out
}
