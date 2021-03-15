\name{ergm-constraints}
\alias{ergm-constraints}
\alias{constraints-ergm}
\alias{ergm.constraints}
\alias{constraints.ergm}

\alias{blockdiag}

\docType{methods}

\title{Sample Space Constraints for Exponential-Family Random Graph Models}
\description{
  \code{\link{ergm}} is used to fit exponential-family random graph
  models (ERGMs), in which
  the probability of a given network, \eqn{y}, on a set of nodes is
  \eqn{h(y) \exp\{\eta(\theta) \cdot
    g(y)\}/c(\theta)}, where
  \eqn{h(y)} is the reference measure (usually \eqn{h(y)=1}),
  \eqn{g(y)} is a vector of network statistics for \eqn{y},
  \eqn{\eta(\theta)} is a natural parameter vector of the same
  length (with \eqn{\eta(\theta)=\theta} for most terms), and \eqn{c(\theta)} is the
  normalizing constant for the distribution.

  This page describes the constraints (the networks \eqn{y} for which \eqn{h(y)>0})
  that are included with the
  \code{\link[=ergm-package]{ergm}} package. Other packages may add new
  constraints.

}

\section{Constraints implemented in  the \code{\link[=tergm-package]{tergm}} package}{
  \describe{
    \item{\code{blockdiag(attrname)} (dyad-independent)}{
      Force a block-diagonal structure (and its bipartite analogue) on
      the network. Only dyads \eqn{(i,j)} for which
      \code{attrname(i)==attrname(j)} can have edges.

      Note that the current implementation requires that blocks be
      contiguous for \dQuote{unipartite} graphs, and for bipartite
      graphs, they must be contiguous within a partition and must have
      the same ordering in both partitions. (They do not, however,
      require that all blocks be represented in both partitions, but
      those that overlap must have the same order.)

      If multiple block-diagonal constraints are given, or if
      \code{attrname} is a vector with multiple attribute names, blocks
      will be constructed on \emph{all} attributes matching.
    }
  }
}

\references{

Goodreau SM, Handcock MS, Hunter DR, Butts CT, Morris M (2008a).
 A \pkg{statnet} Tutorial.
 \emph{Journal of Statistical Software}, 24(8).
 \url{http://www.jstatsoft.org/v24/i08/}.

Hunter, D. R. and Handcock, M. S. (2006)
    \emph{Inference in curved exponential family models for networks},
   Journal of Computational and Graphical Statistics.

Hunter DR, Handcock MS, Butts CT, Goodreau SM, Morris M (2008b).
 \pkg{ergm}: A Package to Fit, Simulate and Diagnose
  Exponential-Family Models for Networks.
 \emph{Journal of Statistical Software}, 24(3).
 \url{http://www.jstatsoft.org/v24/i03/}.

Karwa V, Krivitsky PN, and Slavkovi\'{c} AB (2016). Sharing Social Network
 Data: Differentially Private Estimation of Exponential-Family Random
 Graph Models. \emph{Journal of the Royal Statistical Society, Series
 C}, to appear. \url{http://arxiv.org/abs/1511.02930}

Krivitsky PN (2012). Exponential-Family Random Graph Models for Valued
Networks. \emph{Electronic Journal of Statistics}, 2012, 6,
1100-1128. \doi{10.1214/12-EJS696}

Morris M, Handcock MS, Hunter DR (2008).
 Specification of Exponential-Family Random Graph Models:
  Terms and Computational Aspects.
 \emph{Journal of Statistical Software}, 24(4).
 \url{http://www.jstatsoft.org/v24/i04/}.
  
}

\keyword{models}
