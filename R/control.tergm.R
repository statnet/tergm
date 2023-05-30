#  File R/control.tergm.R in package tergm, part of the
#  Statnet suite of packages for network analysis, https://statnet.org .
#
#  This software is distributed under the GPL-3 license.  It is free,
#  open source, and has the attribution requirements (GPL Section 7) at
#  https://statnet.org/attribution .
#
#  Copyright 2008-2023 Statnet Commons
################################################################################
###########################################################################
# The <control.tergm> function allows the ergm fitting process to be tuned
# by returning a list of several control parameters
############################################################################



#' Auxiliary for Controlling Temporal ERGM Fitting
#' 
#' Auxiliary function as user interface for fine-tuning 'tergm' fitting.
#' 
#' This function is only used within a call to the \code{\link{tergm}}
#' function.  See the \code{usage} section in \code{\link{tergm}} for details.
#' 
#' @param init numeric or \code{NA} vector equal in
#'   length to the number of parameters in the
#'   model or \code{NULL} (the default); the initial values for the
#'   estimation and coefficient offset terms. If \code{NULL} is
#'   passed, all of the initial values are computed using the method
#'   specified by [`control$init.method`][control.ergm].
#'   If a numeric vector is given, the elements of the vector are
#'   interpreted as follows: \itemize{ \item Elements corresponding to
#'   terms enclosed in \code{offset()} are used as the fixed offset
#'   coefficients. These should match the offset values given in
#'   \code{offset.coef}.
#' 
#'   \item Elements that do not correspond to offset terms and are not
#'   \code{NA} are used as starting values in the estimation.
#' 
#'   \item Initial values for the elements that are \code{NA} are fit
#'   using the method specified by
#'   \code{\link[=control.ergm]{control$init.method}}.
#' 
#'   } Passing coefficients from a previous run can be used to
#'   "resume" an uncoverged \code{\link{tergm}} run.
#' @param init.method Estimation method used to acquire initial values
#'   for estimation. If \code{NULL} (the default), the initial values
#'   are computed using the edges dissolution approximation (Carnegie
#'   et al.) when appropriate; note that this relies on \code{\link{.extract.fd.formulae}}
#'   to identify the formation and dissolution parts of the formula; the user should
#'   be aware of its behavior and limitations.
#'   If \code{init.method} is set to "zeros", the initial values are set to zeros.
#' @param force.main Logical: If TRUE, then force MCMC-based
#'   estimation method, even if the exact MLE can be computed via
#'   maximum pseudolikelihood estimation.
#' @param MCMC.prop.weights Specifies the proposal weighting to use.
#' @param MCMC.prop.args A direct way of specifying arguments to the proposal.
#' @param MCMC.prop Hints and/or constraints for selecting and initializing the proposal.
#' @template control_MCMC_maxedges
#' @param MCMC.maxchanges Maximum number of changes permitted to occur during the simulation.
#' @template control_MCMC_packagenames
#' @param CMLE.MCMC.burnin Burnin used in CMLE fitting.
#' @param CMLE.MCMC.interval Number of Metropolis-Hastings steps
#'   between successive draws when running MCMC MLE.
#' @param CMLE.ergm Control parameters used
#'   to fit the CMLE.  See \code{\link{control.ergm}}.
#' @param CMLE.NA.impute In TERGM CMLE, missing dyads in
#'   transitioned-to networks are accommodated using methods of
#'   Handcock and Gile (2009), but a similar approach to
#'   transitioned-from networks requires much more complex methods
#'   that are not, currently, implemented.  \code{CMLE.NA.impute}
#'   controls how missing dyads in transitioned-from networks are be
#'   imputed. See argument \code{imputers} of
#'   \code{\link{impute.network.list}} for details.
#' 
#'   By default, no imputation is performed, and the fitting stops
#'   with an error if any transitioned-from networks have missing
#'   dyads.
#' @param CMLE.term.check.override The method
#'   \code{\link{tergm}{tergm}} uses at this time to fit a series of
#'   more than two networks requires certain assumptions to be made
#'   about the ERGM terms being used, which are tested before a fit is
#'   attempted. This test sometimes fails despite the model being
#'   amenable to fitting, so setting this option to \code{TRUE}
#'   overrides the tests.
#' @param EGMME.main.method Estimation method used to find the
#'   Equilibrium Generalized Method of Moments estimator.  Currently
#'   only "Gradient-Descent" is implemented.
#' @param EGMME.initialfit.control Control object for the ergm fit
#'   in tergm.EGMME.initialfit
#' @param EGMME.MCMC.burnin.min,EGMME.MCMC.burnin.max, Number of
#'   Metropolis-Hastings steps
#'   per time step used in EGMME fitting. By default, this is
#'   determined adaptively by keeping track of increments in the
#'   Hamming distance between the transitioned-from network and the
#'   network being sampled.
#'   Once \code{EGMME.MCMC.burnin.min} steps have elapsed,
#'   the increments are tested against 0, and when their average
#'   number becomes statistically indistinguishable from 0 (with the
#'   p-value being greater than \code{EGMME.MCMC.burnin.pval}), or
#'   \code{EGMME.MCMC.burnin.max} steps are proposed, whichever comes
#'   first, the simulation is stopped after an additional
#'   \code{EGMME.MCMC.burnin.add} times the number of elapsed steps
#'   had been taken. (Stopping immediately would bias the sampling.)
#' 
#'   To use a fixed number of steps, set 
#'   \code{EGMME.MCMC.burnin.min} and \code{EGMME.MCMC.burnin.max} to
#'   the same value.
#' @param EGMME.MCMC.burnin.pval,EGMME.MCMC.burnin.add Number of
#'   Metropolis-Hastings steps 
#'   per time step used in EGMME fitting. By default, this is
#'   determined adaptively by keeping track of increments in the
#'   Hamming distance between the transitioned-from network and the
#'   network being sampled.
#'   Once \code{EGMME.MCMC.burnin.min} steps have elapsed,
#'   the increments are tested against 0, and when their average
#'   number becomes statistically indistinguishable from 0 (with the
#'   p-value being greater than \code{EGMME.MCMC.burnin.pval}), or
#'   \code{EGMME.MCMC.burnin.max} steps are proposed, whichever comes
#'   first, the simulation is stopped after an additional
#'   \code{EGMME.MCMC.burnin.add} times the number of elapsed steps
#'   had been taken. (Stopping immediately would bias the sampling.)
#' 
#'   To use a fixed number of steps, set 
#'   \code{EGMME.MCMC.burnin.min} and \code{EGMME.MCMC.burnin.max} to
#'   the same value.
#' @param SAN.maxit When \code{target.stats} argument is passed to
#' [ergm()], the maximum number of attempts to use \code{\link{san}}
#' to obtain a network with statistics close to those specified.
#' @param SAN.nsteps.times Multiplier for \code{SAN.nsteps} relative to
#' \code{MCMC.burnin}. This lets one control the amount of SAN burn-in
#' (arguably, the most important of SAN parameters) without overriding the
#' other SAN defaults.
#' @param SAN SAN control parameters.  See
#'   \code{\link{control.san}}
#' @param SA.restarts Maximum number of times to restart a failed
#'   optimization process.
#' @param SA.burnin Number of time steps to advance the starting
#'   network before beginning the optimization.
#' @param SA.plot.progress,SA.plot.stats Logical: Plot information
#'   about the fit as it proceeds. If \code{SA.plot.progress==TRUE},
#'   plot the trajectories of the parameters and target statistics as
#'   the optimization progresses. If \code{SA.plot.stats==TRUE}, plot
#'   a heatmap representing correlations of target statistics and a
#'   heatmap representing the estimated gradient.
#' 
#'   Do NOT use these with non-interactive plotting devices like
#'   \code{\link{pdf}}. (In fact, it will refuse to do that with a
#'   warning.)
#' @param SA.max.plot.points If \code{SA.plot.progress==TRUE}, the
#'   maximum number of time points to be plotted. Defaults to 400. If
#'   more iterations elapse, they will be thinned to at most 400
#'   before plotting.
#' @param SA.init.gain Initial gain, the multiplier for the parameter
#'   update size.  If the process initially goes crazy beyond
#'   recovery, lower this value.
#' @param SA.gain.decay Gain decay factor.
#' @param SA.runlength Number of parameter trials and updates per C
#'   run.
#' @param SA.interval.mul The number of time steps between updates of
#'   the parameters is set to be this times the mean duration of
#'   extant ties.
#' @param SA.init.interval Initial number of time steps between
#'   updates of the parameters.
#' @param SA.min.interval,SA.max.interval Upper and lower bounds on
#'   the number of time steps between updates of the parameters.
#' @param SA.phase1.tries Number of runs trying to find a reasonable
#'   parameter and network configuration.
#' @param SA.phase1.jitter Initial jitter standard deviation of each
#'   parameter.
#' @param SA.phase1.max.q Q-value (false discovery rate) that a
#'   gradient estimate must obtain before it is accepted (since sign
#'   is what is important).
#' @param SA.phase1.backoff.rat,SA.phase2.backoff.rat If the run
#'   produces this relative increase in the approximate objective
#'   function, it will be backed off.
#' @param SA.phase1.minruns Number of runs during Phase 1 for
#'   estimating the gradient, before every gradient update.
#' @param SA.phase2.levels.min,SA.phase2.levels.max Range of gain
#'   levels (subphases) to go through.
#' @param SA.phase2.max.mc.se Approximate precision of the estimates
#'   that must be attained before stopping.
#' @param SA.phase2.repeats,SA.stepdown.maxn, A gain level may be
#'   repeated multiple times (up to \code{SA.phase2.repeats}) if the
#'   optimizer detects that the objective function is improving or the
#'   estimating equations are not centered around 0, so slowing down
#'   the parameters at that point is counterproductive. To detect this
#'   it looks at the the window controlled by \code{SA.keep.oh},
#'   thinning objective function values to get
#'   \code{SA.stepdown.maxn}, and 1) fitting a GLS model for a linear
#'   trend, with AR(2) autocorrelation and 2) conductiong an
#'   approximate Hotelling's T^2 test for equality of estimating
#'   equation values to 0. If there is no significance for either at
#'   \code{SA.stepdown.p} \code{SA.stepdown.ct} runs in a row, the
#'   gain level (subphase) is allowed to end. Otherwise, the process
#'   continues at the same gain level.
#' @param SA.stepdown.p,SA.stepdown.ct A gain level may be repeated
#'   multiple times (up to \code{SA.phase2.repeats}) if the optimizer
#'   detects that the objective function is improving or the
#'   estimating equations are not centered around 0, so slowing down
#'   the parameters at that point is counterproductive.  To detect
#'   this it looks at the the window controlled by \code{SA.keep.oh},
#'   thinning objective function values to get
#'   \code{SA.stepdown.maxn}, and 1) fitting a GLS model for a linear
#'   trend, with AR(2) autocorrelation and 2) conductiong an
#'   approximate Hotelling's T^2 test for equality of estimating
#'   equation values to 0. If there is no significance for either at
#'   \code{SA.stepdown.p} \code{SA.stepdown.ct} runs in a row, the
#'   gain level (subphase) is allowed to end. Otherwise, the process
#'   continues at the same gain level.
#' @param SA.stop.p At the end of each gain level after the minimum,
#'   if the precision is sufficiently high, the relationship between
#'   the parameters and the targets is tested for evidence of local
#'   nonlinearity. This is the p-value used.
#' 
#'   If that test fails to reject, a Phase 3 run is made with the new
#'   parameter values, and the estimating equations are tested for
#'   difference from 0. If this test fails to reject, the optimization
#'   is finished.
#' 
#'   If either of these tests rejects, at \code{SA.stop.p},
#'   optimization is continued for another gain level.
#' @param SA.keep.oh,SA.keep.min,SA.keep.min.runs Parameters
#'   controlling how much of optimization history to keep for gradient
#'   and covariance estimation.
#' 
#'   A history record will be kept if it's at least one of the
#'   following: \itemize{ \item Among the last \code{SA.keep.oh} (a
#'   fraction) of all runs.  \item Among the last \code{SA.keep.min} (a
#'   count) records.  \item From the last \code{SA.keep.min.runs} (a
#'   count) optimization runs.  }
#' @param SA.phase2.jitter.mul Jitter standard deviation of each
#'   parameter is this value times its standard deviation without
#'   jitter.
#' @param SA.phase2.maxreljump To keep the optimization from
#'   "running away" due to, say, a poor gradient estimate building on
#'   itself, if a magnitude of change (Mahalanobis distance) in
#'   parameters over the course of a run divided by average magnitude
#'   of change for recent runs exceeds this, the change is truncated
#'   to this amount times the average for recent runs.
#' @param SA.guard.mul The multiplier for the range of parameter and
#'   statistics values to compute the guard width.
#' @param SA.par.eff.pow Because some parameters have much, much
#'   greater effects than others, it improves numerical conditioning
#'   and makes estimation more stable to rescale the \eqn{k}th
#'   estimating function by \eqn{s_k = (\sum_{i=1}^{q}
#'   G_{i,k}^2/V_{i,i})^{-p/2}}, where \eqn{G_{i,k}} is the estimated
#'   gradient of the \eqn{i}th target statistics with respect to
#'   \eqn{k}th parameter. This parameter sets the value of \eqn{p}:
#'   \code{0} for no rescaling, \code{1} (default) for scaling by
#'   root-mean-square normalized gradient, and greater values for
#'   greater penalty.
#' @param SA.robust Whether to use robust linear regression (for
#'   gradients) and covariance estimation.
#' @param SA.oh.memory Absolute maximum number of data points per
#'   thread to store in the full optimization history.
#' @param SA.refine Method, if any, used to refine the point estimate
#'   at the end: "linear" for linear interpolation, "mean" for
#'   average, and "none" to use the last value.
#' @param SA.se Logical: If TRUE (the default), get an MCMC sample of
#'   statistics at the final estimate and compute the covariance
#'   matrix (and hence standard errors) of the parameters. This sample
#'   is stored and can also be used by
#'   [mcmc.diagnostics()] to assess convergence.
#' @param SA.phase3.samplesize.runs This many optimization runs will
#'   be used to determine whether the optimization has converged and
#'   to estimate the standard errors.
#' @param SA.restart.on.err Logical: if \code{TRUE} (the default) an
#'   error somewhere in the optimization process will cause it to
#'   restart with a smaller gain value. Otherwise, the process will
#'   stop. This is mainly used for debugging
#' @template term_options
#' @template seed
#' @template control_MCMC_parallel
#' @param MCMC.burnin,MCMC.burnin.mul No longer used. See
#'   \code{EGMME.MCMC.burnin.min}, \code{EGMME.MCMC.burnin.max},
#'   \code{EGMME.MCMC.burnin.pval}, \code{EGMME.MCMC.burnin.pval},
#'   \code{EGMME.MCMC.burnin.add} and \code{CMLE.MCMC.burnin} and
#'   \code{CMLE.MCMC.interval}.
#' @return A list with arguments as components.
#' @seealso \code{\link{tergm}}. The
#'   \code{\link{control.simulate.tergm}} function performs a similar
#'   function for \code{\link{simulate.tergm}}.
#' @references Boer, P., Huisman, M., Snijders,
#'   T.A.B., and Zeggelink, E.P.H. (2003), StOCNET User\'s
#'   Manual. Version 1.4.
#' 
#' Firth (1993), Bias Reduction in Maximum Likelihood Estimates.
#' Biometrika, 80: 27-38.
#' 
#' Hunter, D. R. and M. S. Handcock (2006), Inference in curved
#' exponential family models for networks. Journal of Computational and
#' Graphical Statistics, 15: 565-583.
#' 
#' Hummel, R. M., Hunter, D. R., and Handcock, M. S. (2010), A Steplength
#' Algorithm for Fitting ERGMs, Penn State Department of Statistics Technical
#' Report.
#' @export control.tergm
control.tergm<-function(init=NULL,
                         init.method=NULL,
                         force.main = FALSE,                         

                         MCMC.prop = ~discord + sparse,

                         MCMC.prop.weights="default",MCMC.prop.args=NULL,
                         MCMC.maxedges=Inf,
                         MCMC.maxchanges=1000000,
                         MCMC.packagenames=c(),

                         CMLE.MCMC.burnin = 1024*16,
                         CMLE.MCMC.interval = 1024,
                         CMLE.ergm=control.ergm(init=init, MCMC.burnin=CMLE.MCMC.burnin, MCMC.interval=CMLE.MCMC.interval, MCMC.prop=MCMC.prop, MCMC.prop.weights=MCMC.prop.weights, MCMC.prop.args=MCMC.prop.args, MCMC.maxedges=MCMC.maxedges, MCMC.packagenames=MCMC.packagenames, parallel=parallel, parallel.type=parallel.type, parallel.version.check=parallel.version.check, force.main=force.main, term.options=term.options),

                         CMLE.NA.impute=c(),
                         CMLE.term.check.override=FALSE,
                         
                         EGMME.main.method=c("Gradient-Descent"),

                         EGMME.initialfit.control = control.ergm(),

                         EGMME.MCMC.burnin.min=1000,
                         EGMME.MCMC.burnin.max=100000,
                         EGMME.MCMC.burnin.pval=0.5,
                         EGMME.MCMC.burnin.add=1,

                         MCMC.burnin=NULL, MCMC.burnin.mul=NULL,
                         
                         SAN.maxit=4,
                         SAN.nsteps.times=8,
                         SAN=control.san(
                           term.options=term.options,
                           SAN.maxit=SAN.maxit,
                           SAN.prop=MCMC.prop,
                           SAN.prop.weights=MCMC.prop.weights,
                           SAN.prop.args=MCMC.prop.args,

                           SAN.nsteps=round(sqrt(EGMME.MCMC.burnin.min*EGMME.MCMC.burnin.max))*SAN.nsteps.times,

                           SAN.packagenames=MCMC.packagenames,
                           
                           parallel=parallel,
                           parallel.type=parallel.type,
                           parallel.version.check=parallel.version.check,
                           parallel.inherit.MT=parallel.inherit.MT),

                         SA.restarts=10,
                         
                         SA.burnin=1000,

                         # Plot the progress of the optimization.
                         SA.plot.progress=FALSE,
                         SA.max.plot.points=400,
                         SA.plot.stats=FALSE,
                         
                         # Initial gain --- if the process initially goes
                         # crazy beyond recovery, lower this.
                         SA.init.gain=0.1,
                         SA.gain.decay=0.5, # Gain decay factor.
                         
                         SA.runlength=25, # Number of jumps per .C call.

                         # Interval --- number of steps between
                         # successive jumps --- is computed
                         # adaptively.
                         SA.interval.mul=2, # Set the mean duration of extant ties this to be the interval.
                         SA.init.interval=500, # Starting interval.
                         SA.min.interval=20, # The lowest it can go.
                         SA.max.interval=500, # The highest it can go.


                         SA.phase1.minruns=4, # Number of pure-jitter runs before gradient calculation and such begin.
                         SA.phase1.tries=20, # Number of iterations of trying to find a reasonable configuration. FIXME: nothing happens if it's exceeded.
                         SA.phase1.jitter=0.1, # Initial jitter sd of each parameter..
                         SA.phase1.max.q=0.1, # FDR q-value for considering a gradient to be significant.
                         SA.phase1.backoff.rat=1.05, # If a run produces this relative increase in the objective function, it will be backed off.
                         SA.phase2.levels.max=40, # Maximum number of gain levels to go through.
                         SA.phase2.levels.min=4, # Minimum number of gain levels to go through.
                         SA.phase2.max.mc.se=0.001, # Maximum Monte-Carlo variation error of the parameter estimates as a fraction of total variation.
                         SA.phase2.repeats=400, # Maximum number of times gain a subphase can be repeated if the optimization is "going somewhere".
                         SA.stepdown.maxn=200, # Thin the draws for trend detection to get this many.
                         SA.stepdown.p=0.05, # If the combined p-value for the trend in the parameters is less than this, reset the subphase counter.
                         SA.stop.p=0.1, # If the combined p-value of the stopping tests exceeds this, finish the optimization.
                         SA.stepdown.ct=5, # Baseline number of times in a row the p-value must be above SA.stepdown.p to reduce gain.
                         SA.phase2.backoff.rat=1.1, # If a run produces this relative increase in the objective function, it will be backed off.
                         SA.keep.oh=0.5, # Fraction of optimization history that is used for gradient and covariance calculation.
                         SA.keep.min.runs=8, # Minimum number of runs that are used for gradient and covariance calculation.
                         SA.keep.min=0, # Minimum number of observations that are used for gradient and covariance calculation.
                         SA.phase2.jitter.mul=0.2, # The jitter standard deviation of each parameter is this times its standard deviation sans jitter.
                         SA.phase2.maxreljump=4, # Maximum jump per run, relative to the magnitude of other jumps in the history.
                         SA.guard.mul = 4, # The multiplier for the range of parameter values to compute the guard width.
                         SA.par.eff.pow = 1, # How do we scale rows of the estimating equation as a function of G?
                         SA.robust = FALSE, # Should the (slower) robust linear regression and covariance estimation be used?
                         SA.oh.memory = 100000, # Maximum number of jumps to store.
                         
                         SA.refine=c("mean","linear","none"), # Method, if any, used to refine the point estimate: linear interpolation, average, and none for the last value.
                         
                         SA.se=TRUE, # Whether to run Phase 3 to compute the standard errors.
                         SA.phase3.samplesize.runs=10, # This times the interval is the number of steps to estimate the standard errors.
                         SA.restart.on.err=TRUE, # Whether to wrap certain routines in try() statements so that they that an error is handled gracefully. Set to FALSE to debug errors in those routines.

                         term.options=NULL,
                         
                         seed=NULL,
                         parallel=0,
                         parallel.type=NULL,
                         parallel.version.check=TRUE,
                         parallel.inherit.MT=FALSE){

  if(!is.null(MCMC.burnin) || !is.null(MCMC.burnin.mul)) stop("Control parameters MCMC.burnin and MCMC.burnin.mul are no longer used. See help for EGMME.MCMC.burnin.min, EGMME.MCMC.burnin.max, EGMME.MCMC.burnin.pval, EGMME.MCMC.burnin.pval, and CMLE.MCMC.burnin and CMLE.MCMC.interval for their replacements.")

  match.arg.pars=c("EGMME.main.method","SA.refine")
  
  control<-list()
  formal.args<-formals(sys.function())
  for(arg in names(formal.args))
    control[arg]<-list(get(arg))

  for(arg in match.arg.pars)
    control[arg]<-list(match.arg(control[[arg]][1],eval(formal.args[[arg]])))

  set.control.class("control.tergm")
}
