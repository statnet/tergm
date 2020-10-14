#  File R/stergm.R in package tergm, part of the Statnet suite
#  of packages for network analysis, https://statnet.org .
#
#  This software is distributed under the GPL-3 license.  It is free,
#  open source, and has the attribution requirements (GPL Section 7) at
#  https://statnet.org/attribution
#
#  Copyright 2008-2020 Statnet Commons
#######################################################################
################################################################################
# stergm --- fit Separable Temporal ERGMs.
################################################################################



#' Separable Temporal Exponential Family Random Graph Models
#' 
#' \code{\link{stergm}} is used for finding Separable Temporal ERGMs' (STERGMs)
#' Conditional MLE (CMLE) (Krivitsky and Handcock, 2010) and Equilibrium
#' Generalized Method of Moments Estimator (EGMME) (Krivitsky, 2009).
#' 
#' \strong{Model Terms} See \code{\link{ergm}} and \code{\link{ergm-terms}} for
#' details. At this time, only linear ERGM terms are allowed.  \itemize{
#' \item For a brief demonstration, please see the tergm package vignette:
#' \code{browseVignettes(package='tergm')} \item A more detailed tutorial is
#' avalible on the statnet wiki:
#' \url{https://statnet.org/Workshops/tergm_tutorial.html} }
#' 
#' @param nw A \code{\link[network]{network}} object (for EGMME); or
#' \code{\link[networkDynamic]{networkDynamic}} object, a
#' \code{\link{network.list}} object, or a \code{\link{list}} containing
#' networks (for CMLE and CMPLE).
#' 
#' \code{stergm} understands the \code{\link{lasttoggle}} "API".
#' @param formation,dissolution One-sided \code{\link{ergm}}-style formulas for
#' the formation and dissolution models, respectively.
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
#' For STERGMs in particular, the constraints apply to the post-formation and
#' the post-dissolution network, rather than the final network. This means, for
#' example, that if the degree of all vertices is constrained to be less than
#' or equal to three, and a vertex begins a time step with three edges, then,
#' even if one of its edges is dissolved during its time step, it won't be able
#' to form another edge until the next time step. This behavior may change in
#' the future.
#' 
#' Note that not all possible combinations of constraints are supported.
#' @param estimate One of "EGMME" for Equilibrium Generalized Method of Moments
#' Estimation, based on a single network with some temporal information and
#' making an assumption that it is a product of a STERGM process running to its
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
#' @param offset.coef.form Numeric vector to specify offset formation
#' parameters.
#' @param offset.coef.diss Numeric vector to specify offset dissolution
#' parameters.
#' @param targets One-sided \code{\link{ergm}}-style formula specifying
#' statistics whose moments are used for the EGMME. Unused for CMLE and CMPLE.
#' Targets is required for EGMME estimation. It may contain any valid ergm
#' terms. If specified as "formation" or "dissolution", it copies the formula
#' from the respective model. Any offset terms are removed automatically.
#' @param target.stats A vector specifying the values of the \code{targets}
#' statistics that EGMME will try to match.  Defaults to the statistics of
#' \code{nw}. Unused for CMLE and CMPLE.
#' @param eval.loglik Whether or not to calculate the log-likelihood
#'   of a CMLE STERGM fit. See \code{\link{ergm}} for details. Can be
#'   set globally via `option(tergm.eval.loglik=...)`, falling back to
#'   `getOption("ergm.eval.loglik")` if not set.
#' @param control A list of control parameters for algorithm tuning.
#' Constructed using \code{\link{control.stergm}}.
#' @param verbose logical or integer; if TRUE or positive, the program will
#' print out progress information. Higher values result in more output.
#' @param \dots Additional arguments, to be passed to lower-level functions.
#' @return \code{\link{stergm}} returns an object of class \code{\link{stergm}}
#' that is a list consisting of the following elements:
#' \item{formation, dissolution}{Formation and dissolution formulas,
#' respectively.}
#' \item{targets}{The targets formula.}
#' \item{target.stats}{The target statistics.}
#' \item{estimate}{The type of estimate.}
#' \item{opt.history}{A
#' matrix containing the full trace of the EGMME optimization process:
#' coefficients tried and target statistics simulated.}
#' \item{sample}{An \code{\link{mcmc}} object containing target
#' statistics sampled at the estimate.}
#' \item{covar}{The full estimated
#' variance-covariance matrix of the parameter estimates for EGMME. (Note that
#' although the CMLE formation parameter estimates are independent of the
#' dissolution parameter estimates due to the separability assumption, this is
#' not necessarily the case for EGMME.) }
#' \item{formation.fit, dissolution.fit}{For CMLE and CMPLE, \code{\link{ergm}} objects from
#' fitting formation and dissolution, respectively. For EGMME, stripped down
#' \code{\link{ergm}}-like lists.}
#' \item{network}{For
#' \code{estimate=="EGMME"}, the original network; for \code{estimate=="CMLE"}
#' or \code{estimate=="CMPLE"}, a \code{\link{network.list}} (a discrete series
#' of networks) to which the model was fit.}
#' \item{control}{The control
#' parameters used to fit the model.} 
#' 
#' See the method \code{\link{print.stergm}} for details on how an
#' \code{\link{stergm}} object is printed.  Note that the method
#' \code{\link{summary.stergm}} returns a summary of the relevant parts of the
#' \code{\link{stergm}} object in concise summary format.
#' @seealso ergm, network, \%v\%, \%n\%, \code{\link{ergm-terms}}
#' @references \itemize{
#'
#' \item Krivitsky P.N. and Handcock M.S. (2014) A Separable Model for Dynamic Networks. \emph{Journal of the Royal Statistical Society, Series B}, 76(1): 29-46. \doi{10.1111/rssb.12014}
#' 
#' \item Krivitsky, P.N. (2012). Modeling of Dynamic Networks based on
#' Egocentric Data with Durational Information. \emph{Pennsylvania State
#' University Department of Statistics Technical Report}, 2012(2012-01).
#' \url{https://web.archive.org/web/20170830053722/https://stat.psu.edu/research/technical-report-files/2012-technical-reports/TR1201A.pdf}
#' 
#' }
#' @examples
#' 
#' \donttest{
#' # EGMME Example
#' par(ask=FALSE)
#' n<-30
#' g0<-network.initialize(n,dir=FALSE)
#' 
#' #                     edges, degree(1), mean.age
#' target.stats<-c(      n*1/2,    n*0.6,        20)
#' 
#' dynfit<-stergm(g0,formation = ~edges+degree(1), dissolution = ~edges,
#'                targets = ~edges+degree(1)+mean.age,
#'                target.stats=target.stats, estimate="EGMME",
#'                control=control.stergm(SA.plot.progress=TRUE))
#' 
#' par(ask=TRUE)
#' mcmc.diagnostics(dynfit)
#' summary(dynfit)
#' 
#' # CMLE Example
#' data(samplk)
#' 
#' # Fit a transition from Time 1 to Time 2
#' samplk12 <- stergm(list(samplk1, samplk2),
#'                    formation=~edges+mutual+transitiveties+cyclicalties,
#'                    dissolution=~edges+mutual+transitiveties+cyclicalties,
#'                    estimate="CMLE")
#' 
#' mcmc.diagnostics(samplk12)
#' summary(samplk12)
#' 
#' # Fit a transition from Time 1 to Time 2 and from Time 2 to Time 3 jointly
#' samplk123 <- stergm(list(samplk1, samplk2, samplk3),
#'                     formation=~edges+mutual+transitiveties+cyclicalties,
#'                     dissolution=~edges+mutual+transitiveties+cyclicalties,
#'                     estimate="CMLE")
#' 
#' mcmc.diagnostics(samplk123)
#' summary(samplk123)
#' }
#' @import network
#' @import networkDynamic
#' @export
stergm <- function(nw, formation, dissolution, constraints = ~., estimate, times=NULL, offset.coef.form=NULL, offset.coef.diss=NULL,
                   targets=NULL, target.stats=NULL,
                   eval.loglik=NVL(getOption("tergm.eval.loglik"), getOption("ergm.eval.loglik")),
                   control=control.stergm(),
                   verbose=FALSE, ...) {
  check.control.class("stergm", "stergm")
  
  if(!is.null(control$seed))  set.seed(as.integer(control$seed))

  estimate <- match.arg(estimate,c("CMLE","CMPLE","EGMME"))
  
  if(!inherits(formation,"formula") || !inherits(dissolution,"formula"))
    stop("Arguments formation and dissolution must be formulas.")

  if(length(formation)==3){
    warning("Formation formula has an LHS, which will be ignored in favor of nw.")
    formation <- formation[c(1,3)] # in a formula f<-y~x, f[1]=~, f[2]=y, and f[3]=x
  }

  if(length(dissolution)==3){
    warning("Dissolution formula has an LHS, which will be ignored in favor of nw.")
    dissolution <- dissolution[c(1,3)] # in a formula f<-y~x, f[1]=~, f[2]=y, and f[3]=x
  }
  
  # lasttoggle
  if(estimate=="EGMME"){
  duration.dependent <- is.lasttoggle(nw,formation,dissolution,targets)
  
  if(duration.dependent)
    nw %n% "lasttoggle" <- NVL(nw %n% "lasttoggle",rep(round(-.Machine$integer.max/2), network.dyadcount(nw)))  else nw %n% "lasttoggle" <- NULL
  }
  
  out <- switch(estimate,
                CMLE=,
                CMPLE=stergm.CMLE(nw, formation, dissolution, constraints,
                  times, offset.coef.form, offset.coef.diss, eval.loglik,
                  estimate, control, verbose),
                EGMME=stergm.EGMME(nw, formation, dissolution, constraints,
                  offset.coef.form, offset.coef.diss,
                  targets, target.stats, estimate, control, verbose)
                  )
  
  
  out$formation <- formation
  out$dissolution <- dissolution
  out$control <- control
  
  class(out)<-"stergm"
  out
}
