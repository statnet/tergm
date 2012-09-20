\name{tergm-terms}
\alias{tergm-terms}
\alias{tergm.terms}
\alias{terms-tergm}
\alias{terms.tergm}
\alias{edge.ages}
\alias{mean.age}
\alias{edges.ageinterval}
\alias{edgecov.ages}
\alias{edgecov.mean.age}
\alias{InitErgmTerm.edge.ages}
\alias{InitErgmTerm.edgecov}
\alias{InitErgmTerm.edgecov.ages}
\alias{InitErgmTerm.edges.ageinterval}
\docType{methods}
\title{Temporally-Sensitive Terms used in Exponential Family Random Graph Models}
\description{
  Unlike ordinary \code{\link{ergm-terms}}, which take only a single
  network as an argument, the terms documented here also take into
  account the "ages" of extant ties in the network: the time elapsed since
  their formation.

  As implemented, many of these terms cannot be used to "drive" the process
  of network evolution, but they can be used as target statistics to
  infer the terms that do. More concretely, they may appear in
  \code{targets=} or \code{monitor=} formulas of the
  \code{\link{stergm}} or \code{\link{simulate.stergm}}, but they may
  not appear in their \code{\link{formation=}} and
  \code{\link{dissolution=}} formulas. These terms are marked as "\strong{Target-only}".
}

\section{Possible terms to represent network statistics}{
  The terms currently available are:
  \describe{

    \item{\code{edges.ageinterval(from, to=+Inf)}}{\emph{Number of edges
      with age falling into a specified range:}
      This term adds one statistic equaling the number of edges for
      which the time elapsed since formation is greater than or equal to
      \code{from} but strictly less than \code{to}. In other words, it
      is in the semiopen interval \code{[from, to)}.}

    \item{\code{edge.ages}}{\emph{Sum of ages of extant
	ties:} \strong{Target-only.}
      This term adds one statistic equaling sum, over all ties
      present in the network, of the amount of time elapsed since
      formation.

      Unlike \code{\link{mean.age}}, this statistic is well-defined on
      an empty network. However, if used as a target, it appears to
      produce highly biased dissolution parameter estimates if the goal
      is to get an intended average duration.
    }

    \item{\code{edgecov.ages(x, attrname=NULL)}}{\emph{Weighted sum of ages of extant
	ties:} \strong{Target-only.}
      This term adds one statistic equaling sum, over all ties
      present in the network, of the amount of time elapsed since
      formation, multiplied by a dyadic covariate. See the help for the
      \code{\link[=ergm-terms]{edgecov}} term for details for
      specifying the covariate.
      
      "Weights" can be negative.
      
      Unlike \code{\link{edgecov.mean.age}}, this statistic is well-defined on
      an empty network. However, if used as a target, it appears to
      produce highly biased dissolution parameter estimates if the goal
      is to get an intended average duration.
    }
    
    \item{\code{edgecov.mean.age(x, attrname=NULL,
	emptyval=0)}}{\emph{Weighted average age of an extant
	tie:} \strong{Target-only.}
      This term adds one statistic equaling the average, over all ties
      present in the network, of the amount of time elapsed since
      formation, weighted by a (nonnegative) dyadic covariate. See the help for the
      \code{\link[=ergm-terms]{edgecov}} term for details for
      specifying the covariate.
      
      The behavior when there are negative weights is undefined.
      
      Because this average is undefined for an empty network (or a
      network all of whose extant edges have been weighted 0), the
      argument \code{emptyval} can be used to specify the value returned
      if this is the case. This is, technically, an arbitrary value, but it should
      not have a substantial effect unless a non-negligible fraction of
      networks at the parameter configuration of interest is empty
      and/or if only a few dyads have nonzero weights.  }
    
    
    \item{\code{mean.age(emptyval=0)}}{\emph{Average age of an extant
	tie:} \strong{Target-only.}
      This term adds one statistic equaling the average, over all ties
      present in the network, of the amount of time elapsed since
      formation.
      
      Because this average is undefined for an empty network, the
      argument \code{emptyval} can be used to specify the value returned
      if it is. This is, technically, an arbitrary value, but it should
      not have a substantial effect unless a non-negligible fraction of
      networks at the parameter configuration of interest is empty.  }
    
    
  }
}
\references{
  
  Krivitsky, P.N. (2012). Modeling of Dynamic Networks based on
  Egocentric Data with Durational Information. \emph{Pennsylvania State
  University Department of Statistics Technical Report},
  2012(2012-01). \url{http://stat.psu.edu/Research/2012-technical-reports}

  Krivitsky, P.N. (2012). Modeling Tie Duration in ERGM-Based Dynamic
  Network Models. \emph{Pennsylvania State University Department of
  Statistics Technical Report},
  2012(2012-02). \url{http://stat.psu.edu/Research/2012-technical-reports}
  
}
\seealso{ergm-terms, stergm, simulate.network, network, \%v\%, \%n\%}
\keyword{models}
