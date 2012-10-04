\name{simulate.stergm}
\alias{simulate.stergm}
\alias{simulate.network}
\alias{simulate.networkDynamic}
\title{
Draw from the distribution of an Separable Temporal Exponential Family Random Graph Model
}
\description{
\code{\link[stats]{simulate}} is used to draw from separable temporal exponential family 
  random network models in their natural parameterizations.
  See \code{\link{stergm}} for more information on these models.
}
\usage{
\method{simulate}{stergm}(object, nsim=1, seed=NULL,
                          coef.form = object$formation.fit$coef,
                          coef.diss = object$dissolution.fit$coef,
                          monitor = object$targets,
                          time.slices, time.burnin=0, time.interval=1,
                          control=control.simulate.stergm(),
                          statsonly=NULL,
                          output=c("networkDynamic", "stats", "changes"),
                          stats.form = FALSE,
                          stats.diss = FALSE,
                          verbose=FALSE, 
                          \dots)
\method{simulate}{network}(object, nsim=1, seed=NULL,
                           formation, dissolution,
                           coef.form, coef.diss,
                           monitor = NULL,
                           time.slices, time.burnin=0, time.interval=1,
                           control=control.simulate.network(),
                           statsonly=NULL,
                           output=c("networkDynamic", "stats", "changes"),
                           stats.form = FALSE,
                           stats.diss = FALSE,
                           verbose=FALSE,
                           \dots)
\method{simulate}{networkDynamic}(object, nsim=1, seed=NULL,
                                  formation = attr(object, "formation"),
                                  dissolution = attr(object, "dissolution"),
                                  coef.form = attr(object, "coef.form"),
                                  coef.diss = attr(object, "coef.diss"),
                                  constraints = attr(object, "constraints"),
                                  monitor = attr(object, "monitor"),
                                  time.slices, time.burnin=0, time.interval=1,
                                  control=control.simulate.network(),
                                  statsonly=NULL,
                                  output=c("networkDynamic", "stats", "changes"),
                                  stats.form = FALSE,
                                  stats.diss = FALSE,
                                  verbose=FALSE,
                                  \dots)
}
\arguments{
  \item{object}{an \R object of type \code{\link{stergm}} giving a model
    fit or of type \code{\link{network}} giving the initial network.
  }
  \item{formation, dissolution}{
    One-sided \code{\link{ergm}}-style formulas for the formation and
    dissolution models, respectively.
  }
  
  \item{nsim}{Number of replications (separate chains of networks) of
    the process to run and return. The \code{\link{networkDynamic}}
    method only supports \code{nsim=1}.
  }
  \item{seed}{Random number integer seed.  
    See \code{\link[base]{set.seed}}.
  }
  
  \item{coef.form}{Parameters for the model from which the 
    post-formation network is drawn.
  }
  \item{coef.diss}{As \code{coef.form}, but for the post-dissolution network.
  }
  \item{monitor}{Either a one-sided formula specifying one or more terms whose
    value is to be monitored, or a string containing \code{"formation"} or
    \code{"dissolution"}, to monitor their respective terms, or
    \code{"all"} to monitor distinct terms from both.
  }
  \item{time.slices}{
    Number of networks or network statistics to return from each replication of the dynamic process. See below for
    return types.
    }
  \item{time.burnin}{
    Number of time steps to discard before starting to collect network
    statistics. Actual network will only be returned if \code{time.burnin==0}.
  }
  \item{time.interval}{
    Number of time steps between successive recordings of network
    statistics. Actual network will only be returned if \code{time.interval==1}.
  }
  \item{control}{
    A list of control parameters for algorithm
    tuning. Constructed using \code{\link{control.simulate.stergm}} or \code{\link{control.simulate.network}}.
  }
  \item{statsonly}{Deprecated in favor of \code{output}.
  }
  \item{output}{A character vector specifying output type: one of
    "networkDynamic" (the default), "stats", and "changes", with partial
    matching allowed. See Value section for details.}
  \item{stats.form, stats.diss}{Logical: Whether to return
    formation/dissolution model statistics. This is not the recommended
    method: use \code{monitor} argument instead.}

  \item{verbose}{Logical:  If TRUE, extra information is printed
    as the Markov chain progresses. 
  }
  \item{\dots}{Further arguments passed to or used by methods.
  }
}
\details{
  The dynamic process is run forward and the results are returned. For
  the method for \code{\link{networkDynamic}}, the simulation is resumed
  from the last generated time point of \code{object}, by default with
  the same model and parameters.
}
\value{
  Depends on the \code{output} argument:
  \item{"stats"}{If \code{stats.form==FALSE} and
    \code{stats.diss==FALSE}, an \code{\link{mcmc}} matrix with
    monitored statistics, and if either of them is \code{TRUE}, a
    list containing elements \code{stats} for statistics specified in the
    \code{monitor} argument, and \code{stats.form} and \code{stats.diss}
    for the respective formation and dissolution statistics.

    If \code{stats.form==FALSE} and
    \code{stats.diss==FALSE} and no monitored statistics are specified,
    an empty list is returned, with a warning.

    When \code{nsim>1}, an \code{\link{mcmc.list}} (or list of them) of
    the statistics is returned instead.}

  \item{"networkDynamic"}{A \code{\link[networkDynamic]{networkDynamic}}
    object representing the simulated process, with ties present in the
    initial network having onset \code{-Inf} and ties present at the end
    of the simulation having terminus \code{+Inf}. The method for
    \code{\link[networkDynamic]{networkDynamic}} returns the initial
    \code{\link[networkDynamic]{networkDynamic}} with simulated changes
    applied to it.

    Additionally, attributes (\code{\link{attr}}, not network
    attributes) are attached as follows:
    \describe{
      \item{formation, dissolution, monitor:}{Formation, dissolution,
	and monitoring formulas	used in the simulation,	respectively.}
      \item{stats, stats.form, stats.diss:}{Network statistics as above.}
      \item{coef.form, coef.diss:}{Coefficients used in the simulation.}
      \item{start, end:}{Start and end times of the simulation.}
    }

    When \code{nsim>1}, a \code{\link{network.list}} of these
    \code{\link{networkDynamic}}s is returned.
  }
  
  \item{"changes"}{An integer matrix with four columns (\code{time},
      \code{tail}, \code{head}, and \code{dir}), giving the timestampled
      changes relative to the initial network. \code{dir} is \code{+1}
      if a tie was formed and \code{-1} if a tie was dissolved. The
      convention for \code{time} is that it gives the time point during
      which the change is effective. For example, a row
      \code{c(1,2,3,+1)} indicates that between time \eqn{0} and
      \eqn{1}, a tie from node \eqn{2} to node \eqn{3} was formed, so
      that it is present at time point \eqn{1}.

      Additionally, same attributes (\code{\link{attr}}, not network
      attributes) as with \code{output=="networkDynamic"} are attached.

      When \code{nsim>1}, a list of these change matrices is returned.}
}

\examples{
logit<-function(p)log(p/(1-p))
coef.form.f<-function(coef.diss,density) -log(((1+exp(coef.diss))/(density/(1-density)))-1)

# Construct a network with 20 nodes and 20 edges
n<-20
target.stats<-edges<-20
g0<-network.initialize(n,dir=TRUE)
g1<-san(g0~edges,target.stats=target.stats,verbose=TRUE)

S<-10

# To get an average duration of 10...
duration<-10
coef.diss<-logit(1-1/duration)

# To get an average of 20 edges...
dyads<-network.dyadcount(g1)
density<-edges/dyads
coef.form<-coef.form.f(coef.diss,density)

# ... coefficients.
print(coef.form)
print(coef.diss)

# Simulate a networkDynamic
dynsim<-simulate(g1,formation=~edges,dissolution=~edges,coef.form=coef.form,coef.diss=coef.diss,time.slices=S,verbose=TRUE)

# "Resume" the simulation.
dynsim2<-simulate(dynsim,time.slices=S,verbose=TRUE)
}