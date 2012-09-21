\name{ergm-terms}
\alias{ergm-terms}
\alias{ergm.terms}
\alias{terms-ergm}
\alias{terms.ergm}
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
  Unlike ordinary \code{\link[ergm]{ergm-terms}}, which take only a single
  network as an argument, the terms documented here also take into
  account the "ages" of extant ties in the network: the time elapsed since
  their formation.

  As implemented, many of these terms cannot be used to "drive" the process
  of network evolution, but they can be used as target statistics to
  infer the terms that do. More concretely, they may appear in
  \code{targets=} or \code{monitor=} formulas of
  \code{\link{stergm}}, \code{\link{simulate.stergm}}, or
  \code{\link[ergm]{summary}} (with an ERGM formula), but they may
  not appear in their \code{\link{formation=}} and
  \code{\link{dissolution=}} formulas. These terms are marked with
  "(target-only)".

  All terms listed here are binary.
}

\section{Terms to represent network statistics included in the
  \code{\link[=tergm-package]{tergm}} pacakge}{
  \describe{

    \item{\code{edges.ageinterval(from, to=+Inf)}}{\emph{Number of edges
      with age falling into a specified range:}
      This term adds one statistic equaling the number of edges for
      which the time elapsed since formation is greater than or equal to
      \code{from} but strictly less than \code{to}. In other words, it
      is in the semiopen interval \code{[from, to)}.}

    \item{\code{edge.ages} (target-only)}{\emph{Sum of ages of extant
	ties:} 
      This term adds one statistic equaling sum, over all ties
      present in the network, of the amount of time elapsed since
      formation.

      Unlike \code{\link{mean.age}}, this statistic is well-defined on
      an empty network. However, if used as a target, it appears to
      produce highly biased dissolution parameter estimates if the goal
      is to get an intended average duration.
    }

    \item{\code{edgecov.ages(x, attrname=NULL)} (target-only)}{\emph{Weighted sum of ages of extant
	ties:} 
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
	emptyval=0)} (target-only)}{\emph{Weighted average age of an extant
	tie:} 
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
    
    
    \item{\code{mean.age(emptyval=0)} (target-only)}{\emph{Average age of an extant
	tie:} 
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
  \itemize{
    \item Handcock M. S., Hunter D. R., Butts C. T., Goodreau S. G.,
      Krivitsky P. N. and Morris M.
(2012). _Fit, Simulate and Diagnose Exponential-Family Models for
Networks_. Version 3.1.  Project home
page at <URL: http://www.statnet.org>, <URL:
CRAN.R-project.org/package=ergm>.
  
\item  Krivitsky, P.N. (2012). Modeling of Dynamic Networks based on
  Egocentric Data with Durational Information. \emph{Pennsylvania State
  University Department of Statistics Technical Report},
  2012(2012-01). \url{http://stat.psu.edu/Research/2012-technical-reports}

\item  Krivitsky, P.N. (2012). Modeling Tie Duration in ERGM-Based Dynamic
  Network Models. \emph{Pennsylvania State University Department of
  Statistics Technical Report},
  2012(2012-02). \url{http://stat.psu.edu/Research/2012-technical-reports}
  }
}
\seealso{\code{\link[ergm:ergm-terms]{ergm-terms} (from the \code{\link[=ergm-package]{ergm}} package), \code{\link{ergm}}, \code{\link{network}}, \code{\link{\%v\%}}, \code{\link{\%n\%}}}
\keyword{models}
