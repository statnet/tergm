\name{control.stergm}
\alias{control.stergm}
\title{Auxiliary for Controlling Separable Temporal ERGM Fitting
}
\description{  Auxiliary function as user interface for fine-tuning 'stergm' fitting.

}
\usage{
          control.stergm(init.form=NULL,
                         init.diss=NULL,
                         init.method=NULL,
                         force.main = FALSE,                         
                         MCMC.prop.weights.form="default",MCMC.prop.args.form=NULL,
                         MCMC.prop.weights.diss="default",MCMC.prop.args.diss=NULL,
                         MCMC.init.maxedges=20000,
                         MCMC.init.maxchanges=20000,
                         MCMC.packagenames=c(),
                         MCMC.burnin=10000,
                         MCMC.burnin.mul=7,
                         MCMC.interval=100,
                         CMLE.control=NULL,
                         CMLE.control.form=control.ergm(init=init.form, MCMC.burnin=MCMC.burnin, MCMC.interval=MCMC.interval, MCMC.prop.weights=MCMC.prop.weights.form, MCMC.prop.args=MCMC.prop.args.form, MCMC.init.maxedges=MCMC.init.maxedges, MCMC.packagenames=MCMC.packagenames, parallel=parallel, parallel.type=parallel.type, parallel.version.check=parallel.version.check, force.main=force.main),
                         CMLE.control.diss=control.ergm(init=init.diss, MCMC.burnin=MCMC.burnin, MCMC.interval=MCMC.interval, MCMC.prop.weights=MCMC.prop.weights.diss, MCMC.prop.args=MCMC.prop.args.diss, MCMC.init.maxedges=MCMC.init.maxedges, MCMC.packagenames=MCMC.packagenames, parallel=parallel, parallel.type=parallel.type, parallel.version.check=parallel.version.check, force.main=force.main),
                         CMLE.NA.impute=c("stop", "previous", "majority", "0", "1"),
                         EGMME.main.method=c("Gradient-Descent"),
                         SAN.maxit=10,
                         SAN.control=control.san(coef=init.form,
                           SAN.prop.weights=MCMC.prop.weights.form,
                           SAN.prop.args=MCMC.prop.args.form,
                           SAN.init.maxedges=MCMC.init.maxedges,
                           SAN.burnin=MCMC.burnin,
                           SAN.packagenames=MCMC.packagenames,
                           parallel=parallel,
                           parallel.type=parallel.type,
                           parallel.version.check=parallel.version.check),
                         SA.restarts=10,
                         SA.burnin=1000,
                         SA.plot.progress=FALSE,
                         SA.max.plot.points=400,
                         SA.plot.stats=FALSE,
                         SA.init.gain=0.1,
                         SA.gain.decay=0.5, 
                         SA.runlength=25, 
                         SA.interval.mul=2, 
                         SA.init.interval=500, 
                         SA.min.interval=20, 
                         SA.max.interval=500, 
                         SA.phase1.minruns=4, 
                         SA.phase1.tries=20, 
                         SA.phase1.jitter=0.1, 
                         SA.phase1.max.q=0.1, 
                         SA.phase1.backoff.rat=1.05, 
                         SA.phase2.levels.max=40, 
                         SA.phase2.levels.min=4, 
                         SA.phase2.max.mc.se=0.001, 
                         SA.phase2.repeats=400, 
                         SA.stepdown.maxn=200, 
                         SA.stepdown.p=0.05, 
                         SA.stop.p=0.1, 
                         SA.stepdown.ct=5, 
                         SA.phase2.backoff.rat=1.1, 
                         SA.keep.oh=0.5, 
                         SA.keep.min.runs=8, 
                         SA.keep.min=0, 
                         SA.phase2.jitter.mul=0.2, 
                         SA.phase2.maxreljump=4, 
                         SA.guard.mul = 4, 
                         SA.robust = FALSE, 
                         SA.refine=c("mean","linear","none"), 
                         SA.se=TRUE, 
                         SA.phase3.samplesize.runs=10, 
                         SA.restart.on.err=TRUE, 
                         seed=NULL,
                         parallel=0,
                         parallel.type=NULL,
                         parallel.version.check=TRUE)
}
\arguments{
  \item{init.form, init.diss}{numeric or \code{NA} vector equal in length to the number of
    parameters in the formation/dissolution model or \code{NULL} (the default); the initial
    values for the estimation and coefficient offset terms. If \code{NULL}
    is passed, all of the initial values are computed using the method specified by
    \code{\link[=control.ergm]{control$init.method}}.
    If a numeric vector is given, the elements of the vector are interpreted as follows:
    \itemize{
      \item{Elements corresponding to
        terms enclosed in \code{offset()} are used as the fixed offset
        coefficients.
      }
      \item{Elements that do not correspond to offset terms and are not
        \code{NA} are used as starting values in the estimation. 
      }
      \item{Initial values for the elements that are \code{NA} are fit
        using the method specified by
        \code{\link[=control.ergm]{control$init.method}}.
      }
    }
    Passing \code{control.ergm(init=coef(prev.fit))} can be used to ``resume''
    an uncoverged \code{\link{ergm}} run, but see
    \code{\link{enformulate.curved}}.
  }
  \item{init.method}{Estimation method
    used to acquire initial values for estimation. Unused at this time.
  }
  \item{force.main}{Logical:  If TRUE, then force MCMC-based estimation
    method, even if the exact MLE can be computed via maximum pseudolikelihood
    estimation.
  }
  \item{MCMC.prop.weights.form, MCMC.prop.weights.diss}{Specifies the method to allocate probabilities of
    being proposed to dyads in the formation/dissolution phase. Defaults to \code{"default"}, which picks a
    reasonable default for the specified constraint. Possible values include
    \code{"TNT"}, \code{"random"}, though not
    all values may be used
    with all possible constraints.
  }
  \item{MCMC.prop.args.form, MCMC.prop.args.diss}{
    An alternative, direct way of specifying additional arguments to
    the proposal in the formation/dissolution phase.
  }
  \item{MCMC.init.maxedges}{Maximum number 
    of edges for which to allocate space.
  }
  \item{MCMC.init.maxchanges}{Maximum number of changes in dynamic network
    simulation for which to allocate space.
  }
  \item{MCMC.packagenames}{Names of packages in which to look for change statistic
    functions in addition to those autodetected. This argument should not be
    needed outside of very strange setups.
  }
  \item{MCMC.burnin}{Maximum number of Metropolis-Hastings steps per phase
    (formation and dissolution) per time step.
  }
  \item{MCMC.burnin.mul}{Multiplier for the number of "changeables" to
    determine the Metropolis-Hastings steps per phase (formation and
    dissolution) per time steps. (I.e., for formation, the number of
    steps will be set to
    \code{(network.dyadcount(nw0)-network.edgecount(nw0))*control$MCMC.burnin.mul},
    while for dissolution, it will be set to
    \code{network.edgecount(nw0)*control$MCMC.burnin.mul}, or
    \code{control$MCMC.burnin}, whichever is smaller.)
  }
  \item{MCMC.interval}{Number of Metropolis-Hastings steps between
    successive draws when running MCMC MLE.}
  
  \item{CMLE.control}{A convenience argument for specifying both
    \code{CMLE.control.form} and \code{CMLE.control.diss} at once.  See
    \code{\link{control.ergm}}.
  }
  \item{CMLE.control.form, CMLE.control.diss}{Control parameters used to
    fit the CMLE for the formation/dissolution ERGM.  See
    \code{\link{control.ergm}}.
  }
  \item{CMLE.NA.impute}{In STERGM CMLE, missing dyads in transitioned-to
    network can be accommodated using methods of Handcock and Gile
    (2009), but transitioned-from network requires much more methods
    that are not, currently, implemented. \code{CMLE.NA.impute} controls
    the way in which they can be imputed:
    \describe{
      \item{\code{"stop"}}{If there are missing dyads in a
	transitioned-from network stop with an error.}
      \item{\code{"previous"}}{Impute the state of the same dyad in the
	previous time step's network (or earlier, if that one is also
	missing). If the first network has missing dyads, stop with an
	error. This imputation method is likely to lead to an
	underestimation of the formation and dissolution rates.}
      \item{\code{"majority"}}{Impute the missing dyad with the value of
	the majority among the non-missing dyads in that time step's
	network.}
      \item{\code{"0"}}{Assume missing dyads are all non-ties.}
      \item{\code{"1"}}{Assume missing dyads are all ties.}
    }
  }
    
  \item{EGMME.main.method}{Estimation method
    used to find the Equilibrium Generalized Method of Moments
    estimator.  Currently only "Robbins-Monro" is implemented.
  }
  \item{SAN.maxit}{ When \code{target.stats} argument is passed to
  \code{\link{ergm}}, the maximum number of attempts to use
  \code{\link{san}} to obtain a network with statistics close to those
  specified.
  }
  \item{SAN.control}{SAN control parameters.  See
    \code{\link{control.san}}
  }

  \item{SA.restarts}{
    Maximum number of times to restart a failed optimization process. 
  }
  
  \item{SA.burnin}{Number of time steps to advance the starting network
    before beginning the optimization.
  }
  \item{SA.plot.progress, SA.plot.stats}{Logical: Plot information about
    the fit as it proceeds. If \code{SA.plot.progress==TRUE}, plot the trajectories of the
    parameters and target statistics as the
    optimization progresses. If \code{SA.plot.stats==TRUE}, plot
    a heatmap reprsenting correlations of target statistics and a
    heatmap representing the estimated gradient.

    Do NOT use these with non-interactive plotting devices
    like \code{\link{pdf}}. (In fact, it will
    refuse to do that with a warning.)
  }

  \item{SA.max.plot.points}{
    If \code{SA.plot.progress==TRUE}, the maximum number of time points
    to be plotted. Defaults to 400. If more iterations elapse, they will
    be thinned to at most 400 before plotting.
    }
  \item{SA.init.gain}{Initial gain, the multiplier for the parameter update size.  If the process initially goes
    crazy beyond recovery, lower this value.
  }
  \item{SA.gain.decay}{Gain decay factor.
  }
  \item{SA.runlength}{Number of parameter trials and updates per C run.
  }
  \item{SA.interval.mul}{The number of time steps between updates of the
    parameters is set to be this times the mean duration of extant ties.}
  \item{SA.init.interval}{Initial number of time steps between updates
    of the parameters.
  }
  \item{SA.min.interval, SA.max.interval}{Upper and lower bounds on the number of time steps between updates
    of the parameters.
  }
  
  \item{SA.phase1.tries}{Number of runs trying to find a 
    reasonable parameter and network configuration.
  }
  \item{SA.phase1.jitter}{Initial jitter standard deviation
    of each parameter.
  }
  \item{SA.phase1.max.q}{Q-value (false discovery rate) that a gradient estimate must 
    obtain before it is accepted (since sign is what is important).
  }

  \item{SA.phase1.backoff.rat, SA.phase2.backoff.rat}{If the run produces this 
    relative increase in the approximate objective function, it will be backed off.
  }

  \item{SA.phase1.minruns}{Number of runs during Phase 1 for estimating
  the gradient, before every gradient update.
  }
  
  \item{SA.phase2.levels.min, SA.phase2.levels.max}{Range of gain levels (subphases) to go through.
  }

  \item{SA.phase2.max.mc.se}{Approximate precision of the estimates that
    must be attained before stopping. 
  }
  
  \item{SA.phase2.repeats, SA.stepdown.maxn,
    SA.stepdown.p, SA.stepdown.ct}{A gain level may be repeated multiple times (up to
    \code{SA.phase2.repeats}) if the optimizer detects that the
    objective function is improving or the estimating equations are not
    centered around 0, so slowing down the parameters at that point is counterproductive. To detect this
    it looks at the the window controlled by \code{SA.keep.oh}, thinning
    objective function values
    to get \code{SA.stepdown.maxn}, and 1) fitting a GLS model for a linear trend, with AR(2)
    autocorrelation and 2) conductiong an approximate Hotelling's T^2
    test for equality of estimating equation values to 0. If there is no
    significance for either at
    \code{SA.stepdown.p} \code{SA.stepdown.ct} runs in a row, the
    gain level (subphase) is allowed to end. Otherwise, the process
    continues at the same gain level.
  }

  \item{SA.stop.p}{At the end of each gain level after the minimum, if
    the precision is sufficiently high, the relationship between the
    parameters and the targets is tested for evidence of local
    nonlinearity. This is the p-value used.

    If that test fails to reject, a Phase 3 run is made with the new
    parameter values, and the estimating equations are tested for
    difference from 0. If this test fails to reject, the optimization is
    finished.

    If either of these tests rejects, at \code{SA.stop.p}, optimization
    is continued for another gain level.
  }
    
  
  \item{SA.keep.oh, SA.keep.min, SA.keep.min.runs}{Parameters
    controlling how much of
    optimization history to keep for gradient and covariance estimation.

    A history record will be kept if it's at least one of the following:
    \itemize{
      \item{Among the last \code{SA.keep.oh} (a fraction) of all runs.}
      \item{Among the last \code{SA.keep.min} (a count) records.}
      \item{From the last \code{SA.keep.min.runs} (a count) optimization
	runs.}
    }
  }
  \item{SA.phase2.jitter.mul}{Jitter standard deviation of each 
    parameter is this value times its standard deviation without jitter.
  }

  \item{SA.phase2.maxreljump}{To keep the optimization from "running
    away" due to, say, a poor gradient estimate building on itself, if a
    magnitude of change (Mahalanobis distance) in parameters over the course of a run divided by
    average magnitude of change for recent runs exceeds this, the change
    is truncated to this amount times the average for recent runs.
  }
  
  \item{SA.guard.mul}{
    The multiplier for the range of parameter and statistics values to compute the guard
    width.
  }

  \item{SA.robust}{
    Whether to use robust linear regression (for gradients) and
    covariance estimation.
  }
  
  \item{SA.refine}{Method, if any, used to 
    refine the point estimate at the end: "linear" for linear interpolation, 
    "mean" for average, and "none" to use the last value.
  }
  \item{SA.se}{Logical:  If TRUE (the default), get an MCMC sample of statistics at
    the final estimate and compute the 
    covariance matrix (and hence standard errors) of the
    parameters. This sample is stored and can also be used by
    \code{\link{mcmc.diagnostics.stergm}} to assess convergence.
  }
  \item{SA.phase3.samplesize.runs}{This many optimization runs will be
    used to determine whether the optimization has converged and to
    estimate the standard errors.
  }

  \item{SA.restart.on.err}{Logical: if \code{TRUE} (the default) an error
  somewhere in the optimization process will cause it to restart with a
  smaller gain value. Otherwise, the process will stop. This is mainly
  used for debugging
  }
  
  \item{seed}{Seed value (integer) for the random number generator.  
    See \code{\link[base]{set.seed}}
  }
  \item{parallel}{Number of threads in which to run the
    sampling. Defaults to 0 (no parallelism). See the entry on
    \link[=ergm-parallel]{parallel processing} for
    details and troubleshooting.
  }
  \item{parallel.type}{API to use for parallel
    processing. Supported values are \code{"MPI"} and
    \code{"SOCK"}. Defaults to using the \code{snow} package
    default.
  }
  \item{parallel.version.check}{Logical:  If TRUE, check that the version of
    \code{\link[=ergm-package]{ergm}} running on the slave nodes is the
    same as that running on the master node.
  }
}
\details{
      This function is only used within a call to the \code{\link{stergm}} function.
    See the \code{usage} section in \code{\link{stergm}} for details.
}
\value{
  A list with arguments as components.
}
\references{
    \itemize{
    \item
    Boer, P., Huisman, M., Snijders, T.A.B., and Zeggelink, E.P.H. (2003),
    StOCNET User\'s Manual. Version 1.4.
    
    \item Firth (1993),
    Bias Reduction in Maximum Likelihood Estimates.
    Biometrika,
    80: 27-38.

    \item Hunter, D. R. and M. S. Handcock (2006), Inference in curved
    exponential family models for networks. Journal of Computational
    and Graphical Statistics, 15: 565-583.
    
    \item Hummel, R. M., Hunter, D. R., and Handcock, M. S. (2010),
    A Steplength Algorithm for Fitting ERGMs, Penn State Department
    of Statistics Technical Report.
  }
}
\seealso{
  \code{\link{stergm}}. The \code{\link{control.simulate.stergm}} 
function performs a 
similar function for
\code{\link{simulate.stergm}}.  
}
