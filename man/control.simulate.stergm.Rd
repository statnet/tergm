%  File man/control.simulate.stergm.Rd in package tergm, part of the Statnet suite
%  of packages for network analysis, http://statnet.org .
%
%  This software is distributed under the GPL-3 license.  It is free,
%  open source, and has the attribution requirements (GPL Section 7) at
%  http://statnet.org/attribution
%
%  Copyright 2003-2013 Statnet Commons
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
\name{control.simulate.stergm}
\alias{control.simulate.stergm}
\alias{control.simulate.network}
\title{ Auxiliary for Controlling Separable Temporal ERGM Simulation }
\description{
  Auxiliary function as user interface for fine-tuning STERGM simulation.
}
\usage{
control.simulate.stergm(MCMC.burnin.min=NULL, MCMC.burnin.max=NULL,
                        MCMC.burnin.pval=NULL, MCMC.burnin.add=NULL,
                        MCMC.burnin=NULL, MCMC.burnin.mul=NULL,
                        MCMC.prop.weights.form=NULL,MCMC.prop.args.form=NULL,
                        MCMC.prop.weights.diss=NULL,MCMC.prop.args.diss=NULL,
                        MCMC.init.maxedges=NULL,
                        MCMC.packagenames=NULL,
                        MCMC.init.maxchanges=NULL)

control.simulate.network(MCMC.burnin.min=1000, MCMC.burnin.max=100000,
                         MCMC.burnin.pval=0.5, MCMC.burnin.add=1,
                         MCMC.burnin=NULL, MCMC.burnin.mul=NULL,
                         MCMC.prop.weights.form="default",MCMC.prop.args.form=NULL,
                         MCMC.prop.weights.diss="default",MCMC.prop.args.diss=NULL,
                         MCMC.init.maxedges=20000,
                         MCMC.packagenames=c(),
                         MCMC.init.maxchanges=1000000)

}
\details{
    This function is only used within a call to the \code{\link{simulate}} function.
    See the \code{usage} section in \code{\link{simulate.stergm}} for details.
}
\arguments{
  \item{MCMC.burnin.min, MCMC.burnin.max, MCMC.burnin.pval, MCMC.burnin.add}{Number of
    Metropolis-Hastings steps per phase (formation and dissolution) per
    time step used in simulation. By default, this is
    determined adaptively by keeping track of increments in the Hamming
    distance between the transitioned-from network and the network being
    sampled (formation network or dissolution network). Once
    \code{MCMC.burnin.min} steps have elapsed, the increments are tested
    against 0, and when their average number becomes statistically
    indistinguishable from 0 (with the p-value being greater than
    \code{MCMC.burnin.pval}), or \code{MCMC.burnin.max} steps are
    proposed, whichever comes first, the simulation is stopped after an
    additional \code{MCMC.burnin.add} times the number of elapsed
    steps had been taken. (Stopping immediately would bias the sampling.)

    To use a fixed number of steps, set both \code{MCMC.burnin.min} and
    \code{MCMC.burnin.max} to the desired number of steps.
  }

  \item{MCMC.prop.weights.form, MCMC.prop.weights.diss}{Specifies the proposal distribution used in the MCMC
    Metropolis-Hastings algorithm for formation and dissolution, respectively. Possible choices are \code{"TNT"} or 
    \code{"random"}; the \code{"default"}.  
    The \code{TNT} (tie / no tie) option puts roughly equal weight on selecting a
    dyad with or without a tie as a candidate for toggling, whereas the
    \code{random} option puts equal weight on all possible dyads, though the
    interpretation of \code{random} may change according to the constraints in
    place.  When no constraints are in place, the default is TNT, which 
    appears to improve Markov chain mixing particularly for networks with a
    low edge density, as is typical of many realistic social networks.
  }
  \item{MCMC.prop.args.form, MCMC.prop.args.diss}{An alternative, direct way of specifying additional arguments to proposals. 
  }
  \item{MCMC.init.maxchanges}{Maximum number 
    of toggles changes for which to allocate space.
  }
  \item{MCMC.packagenames}{Names of packages in which to look for change statistic
    functions in addition to those autodetected. This argument should not be
    needed outside of very strange setups.
  }
  \item{MCMC.init.maxedges}{Maximum number of edges expected in network.
  }

  \item{MCMC.burnin, MCMC.burnin.mul}{No longer used. See
    \code{MCMC.burnin.min}, \code{MCMC.burnin.max},
    \code{MCMC.burnin.pval}, \code{MCMC.burnin.pval}, and
    \code{MCMC.burnin.add}.
  }
}
\value{
  A list with arguments as components.
}
\seealso{ \code{\link{simulate.stergm}}, \code{\link{simulate.formula}}. 
\code{\link{control.stergm}} performs a 
similar function for
\code{\link{stergm}}.   }
%\keyword{simulation}
\keyword{models}
