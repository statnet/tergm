#  File R/stergm.R in package tergm, part of the
#  Statnet suite of packages for network analysis, https://statnet.org .
#
#  This software is distributed under the GPL-3 license.  It is free,
#  open source, and has the attribution requirements (GPL Section 7) at
#  https://statnet.org/attribution .
#
#  Copyright 2008-2023 Statnet Commons
################################################################################
################################################################################
# stergm --- fit Separable Temporal ERGMs.
################################################################################



#' Separable Temporal Exponential Family Random Graph Models (Deprecated)
#' 
#' \code{\link{stergm}} is used for finding Separable Temporal ERGMs'
#' (STERGMs) Conditional MLE (CMLE) (Krivitsky and Handcock, 2014) and
#' Equilibrium Generalized Method of Moments Estimator (EGMME)
#' (Krivitsky, 2009). This function is deprecated in favor of
#' [tergm()], whose special case it is, and may be removed in a future
#' version.
#' 
#' The \code{stergm} function uses a pair of formulas, \code{formation} and 
#' \code{dissolution} to model tie-dynamics.  The dissolution formula, however, is 
#' parameterized in terms of tie persistence: negative coefficients imply lower 
#' rates of persistence and postive coefficients imply higher rates.  
#' The dissolution effects are simply the negation of these coefficients, but
#' the discrepancy between the terminology and interpretation has always been
#' unfortunate, and we have fixed this in the new \code{tergm} function.
#' 
#' If you are making the transition from old \code{stergm} to new \code{tergm}, note that
#' the \code{dissolution} formula in \code{stergm} maps to the new \code{Persist()} 
#' operator in the \code{tergm} function, NOT the \code{Diss()} operator.
#'
#' @param nw A \code{\link[network]{network}} object (for EGMME); or
#' \code{\link[networkDynamic]{networkDynamic}} object, a
#' \code{\link{network.list}} object, or a \code{\link{list}} containing
#' networks (for CMLE and CMPLE).
#' 
#' \code{stergm} understands the \code{\link{lasttoggle}} "API".
#' @param formation,dissolution One-sided \code{\link{ergm}}-style formulas for
#' the formation and dissolution models, respectively.  In \code{stergm}, 
#' the dissolution formula is parameterized in
#' terms of tie persistence: negative coefficients imply lower rates of persistence
#' and postive coefficients imply higher rates.  The dissolution effects are simply the
#' negation of these coefficients.
#'
#' @template constraints
#'
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
#'   treated as temporally adjacent. Irregularly spaced time series
#'   are not supported at this time.
#' 
#' @param offset.coef.form Numeric vector to specify offset formation
#' parameters.
#' @param offset.coef.diss Numeric vector to specify offset dissolution
#' parameters.
#' @param targets One-sided \code{\link{ergm}}-style formula specifying
#' statistics whose moments are used for the EGMME. Unused for CMLE and CMPLE.
#' Targets is required for EGMME estimation. It may contain any valid ergm
#' terms.  Any offset terms are used only during the 
#' preliminary SAN run; they are removed automatically for the EGMME proper.
#' If \code{targets} is specified as a character
#' (one of \code{"formation"} and \code{"dissolution"}) then
#' the function \code{\link{.extract.fd.formulae}} is used to determine the
#' corresponding formula; the user should be aware of its behavior and limitations.
#' @param SAN.offsets Offset coefficients (if any) to use during the SAN run.
#' @param target.stats A vector specifying the values of the \code{targets}
#' statistics that EGMME will try to match.  Defaults to the statistics of
#' \code{nw}. Unused for CMLE and CMPLE.
#' @param eval.loglik Whether or not to calculate the log-likelihood
#'   of a CMLE STERGM fit. See \code{\link{ergm}} for details. Can be
#'   set globally via `option(tergm.eval.loglik=...)`, falling back to
#'   `getOption("ergm.eval.loglik")` if not set.
#' @param control A list of control parameters for algorithm tuning.
#' Constructed using \code{\link{control.stergm}}.  Remapped to 
#' \code{\link{control.tergm}}.
#' @template verbose
#' @param \dots Additional arguments, to be passed to lower-level functions.
#' @return \code{\link{stergm}} returns an object of class [`tergm`];
#'         see [tergm()] for details and methods.
#'
#' @seealso ergm, network, \%v\%, \%n\%, \code{\link{ergm-terms}}
#' @references
#'
#' Krivitsky P.N. and Handcock M.S. (2014) A Separable Model for Dynamic Networks. \emph{Journal of the Royal Statistical Society, Series B}, 76(1): 29-46. \doi{10.1111/rssb.12014}
#' 
#' Krivitsky, P.N. (2012). Modeling of Dynamic Networks based on
#' Egocentric Data with Durational Information. \emph{Pennsylvania State
#' University Department of Statistics Technical Report}, 2012(2012-01).
#' \url{https://web.archive.org/web/20170830053722/https://stat.psu.edu/research/technical-report-files/2012-technical-reports/TR1201A.pdf}
#' 
#' @import network
#' @import networkDynamic
#' @export
stergm <- function(nw, formation, dissolution, constraints = ~., estimate, times=NULL, offset.coef.form=NULL, offset.coef.diss=NULL,
                   targets=NULL, target.stats=NULL,
                   eval.loglik=NVL(getOption("tergm.eval.loglik"), getOption("ergm.eval.loglik")),
                   control=control.stergm(),
                   verbose=FALSE, ..., SAN.offsets = NULL) {
  .Deprecate_once("tergm")

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

  ## need to make sure offsets and inits are set up properly for the tergm call below
  if(!is.null(control$init.form) || 
     !is.null(control$init.diss) || 
     (estimate %in% c("CMLE", "CMPLE") && (!is.null(control$CMLE.form.ergm$init) || 
                                           !is.null(control$CMLE.diss.ergm$init)))) {
    if(estimate == "EGMME") {
      nw_stergm <- nw
      term.options <- control$term.options
      form_model <- ergm_model(formation, nw = nw_stergm, term.options = term.options, dynamic = TRUE, ...)
      diss_model <- ergm_model(dissolution, nw = nw_stergm, term.options = term.options, dynamic = TRUE, ...)
      init.form <- NVL(control$init.form, rep(NA, nparam(form_model, canonical = FALSE)))
      init.diss <- NVL(control$init.diss, rep(NA, nparam(diss_model, canonical = FALSE)))
    } else {
      if(!is(nw, "tergm_NetSeries")) {
        if(inherits(nw, "network.list") || (is.list(nw) && !is.network(nw) && is.network(nw[[1]]))) {
          nw_stergm <- NetSeries(nw, NA.impute=control$CMLE.NA.impute)
        } else if(inherits(nw,"networkDynamic")) {
          nw_stergm <- NetSeries(nw, times, NA.impute=control$CMLE.NA.impute)
        } else {
          stop("Unsupported specification for the network series. See help for ",sQuote("NetSeries")," for arguments.")
        }
      } else {
        nw_stergm <- nw
      }
      term.options <- control$CMLE.form.ergm$term.options
      form_model <- ergm_model(formation, nw = nw_stergm, term.options = term.options, ...)
      diss_model <- ergm_model(dissolution, nw = nw_stergm, term.options = term.options, ...)
      init.form <- NVL(control$CMLE.form.ergm$init, control$init.form, rep(NA, nparam(form_model, canonical = FALSE)))
      init.diss <- NVL(control$CMLE.diss.ergm$init, control$init.diss, rep(NA, nparam(diss_model, canonical = FALSE)))
    }

    if(length(init.form) != nparam(form_model, canonical = FALSE)) {
      stop("Incorrect length of init.form passed to stergm(); expected ", nparam(form_model, canonical = FALSE), ", got ", length(init.form), ".")
    }

    if(length(init.diss) != nparam(diss_model, canonical = FALSE)) {
      stop("Incorrect length of init.diss passed to stergm(); expected ", nparam(diss_model, canonical = FALSE), ", got ", length(init.diss), ".")
    }
    
    if(!is.null(offset.coef.form)) {
      if(length(offset.coef.form) != sum(form_model$etamap$offsettheta)) {
        stop("Incorrect length of offset.coef.form passed to stergm(); expected ", sum(form_model$etamap$offsettheta), ", got ", length(offset.coef.form), ".")
      }
      init.form[form_model$etamap$offsettheta] <- offset.coef.form
    }
    if(!is.null(offset.coef.diss)) {
      if(length(offset.coef.diss) != sum(diss_model$etamap$offsettheta)) {
        stop("Incorrect length of offset.coef.diss passed to stergm(); expected ", sum(diss_model$etamap$offsettheta), ", got ", length(offset.coef.diss), ".")
      }      
      init.diss[diss_model$etamap$offsettheta] <- offset.coef.diss
    }
    
    offset.coef.form <- init.form[form_model$etamap$offsettheta]
    offset.coef.diss <- init.diss[diss_model$etamap$offsettheta]
    
    if(any(is.na(offset.coef.form))) {
      stop("Formation model contains offsets whose coefficients have not been specified.")
    }

    if(any(is.na(offset.coef.diss))) {
      stop("Dissolution model contains offsets whose coefficients have not been specified.")
    }
    
    control$init.form <- init.form
    control$init.diss <- init.diss
  }
  
  control$init <- c(control$init.form, control$init.diss) 

  control$MCMC.prop <- control$MCMC.prop.form
  control$MCMC.prop.weights <- control$MCMC.prop.weights.form
  control$MCMC.prop.args <- control$MCMC.prop.args.form

  control$CMLE.ergm <- control$CMLE.form.ergm
  control$CMLE.ergm$init <- control$init
  
  control <- set.control.class("control.tergm")
  
  formula <- nw ~ Form(formation) + Persist(dissolution)
  
  tergm(formula=formula, constraints=constraints, estimate=estimate, times=times, offset.coef=c(offset.coef.form, offset.coef.diss), targets=targets, target.stats=target.stats, eval.loglik=eval.loglik, control=control, verbose=verbose, SAN.offsets = SAN.offsets, ...)
}
