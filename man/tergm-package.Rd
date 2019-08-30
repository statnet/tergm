%  File man/tergm-package.Rd in package tergm, part of the Statnet suite
%  of packages for network analysis, https://statnet.org .
%
%  This software is distributed under the GPL-3 license.  It is free,
%  open source, and has the attribution requirements (GPL Section 7) at
%  https://statnet.org/attribution
%
%  Copyright 2008-2019 Statnet Commons
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
\name{tergm-package}
\alias{tergm-package}
\alias{tergm}
\docType{package}
\title{
Fit, Simulate and Diagnose Dynamic Network Models derived from
Exponential-Family Random Graph Models
}
\description{
  
  \code{\link[=tergm-package]{tergm}} is a collection of extensions to
  the \code{\link[=ergm-package]{ergm}} package to fit, diagnose, and
  simulate models for dynamic networks --- networks that evolve over
  time --- based on exponential-family random graph models (ERGMs). For
  a list of functions type \code{help(package='tergm')}

  When publishing results obtained using this package, please cite the
  original authors as described in \code{citation(package="tergm")}.

  All programs derived from this package must cite it.
}
\details{

  An exponential-family random graph model (ERGM) postulates an
  exponential family over the sample space of networks of interest, and
  \code{\link[=ergm-package]{ergm}} package implements a suite of tools
  for modeling single networks using ERGMs.

  More recently, there has been a number of extensions of ERGMs to model
  evolution of networks, including the temporal ERGM (TERGM) of Hanneke
  et al. (2010) and the separable termporal ERGM (STERGM) of Krivitsky
  and Handcock (2013). The latter model allows familiar ERGM terms and
  statistics to be reused in a dynamic context, interpreted in terms of
  formation and dissolution of ties. Krivitsky (2012) suggested a method
  for fitting dyanmic models when only a cross-sectional network is
  available, provided some temporal information for it is available as
  well.

  This package aims to implement these and other ERGM-based models for
  network evoluation. At this time, it implements, via the
  \code{\link{stergm}} function, the STERGMs, both a conditional MLE
  (CMLE) fitting to a series of networks and an Equilibrium Generalized
  Method of Moments Estimation (EGMME) for fitting to a single network
  with temporal information. For further development, see the referenced
  papers.
  
  For detailed information on how to download and install the software,
  go to the Statnet project website:
  \url{http://statnet.org}.
  A tutorial, support newsgroup, references and links to further resources are provided there.
}
\author{
  Pavel N. Krivitsky \email{krivitsky@stat.psu.edu} and \cr
  Mark S. Handcock \email{handcock@stat.ucla.edu},\cr
  with contributions from \cr
  David R. Hunter \email{dhunter@stat.psu.edu},\cr
  Steven M. Goodreau \email{goodreau@u.washington.edu},\cr
  Martina Morris \email{morrism@u.washington.edu},\cr
  Nicole Bohme Carnegie \email{nicole.carnegie@nyu.edu}, and\cr
  Ayn Leslie-Cook \email{aynlc3@uw.edu}

  Maintainer: Pavel N. Krivitsky \email{krivitsky@stat.psu.edu}
}
\references{
  \itemize{
  \item Hanneke S, Fu W, and Xing EP (2010). Discrete Temporal Models
  of Social Networks. \emph{Electronic Journal of Statistics}, 2010, 4,
  585-605. \href{http://dx.doi.org/10.1214/09-EJS548}{\code{doi:10.1214/09-EJS548}}
  
  \item Krivitsky PN & Handcock MS (2014) A Separable Model for Dynamic Networks. \emph{Journal of the Royal Statistical Society, Series B}, 76(1): 29-46. \doi{10.1111/rssb.12014}

\item  Krivitsky, P.N. (2012). Modeling of Dynamic Networks based on
  Egocentric Data with Durational Information. \emph{Pennsylvania State
  University Department of Statistics Technical Report},
  2012(2012-01). \url{http://stat.psu.edu/research/technical-report-files/2012-technical-reports/modeling-of-dynamic-networks-based-on-egocentric-data-with-durational-information}


\item Butts CT (2008).
 \pkg{network}: A Package for Managing Relational Data in \R.
 \emph{Journal of Statistical Software}, 24(2).
 \url{http://www.jstatsoft.org/v24/i02/}.

 \item Goodreau SM, Handcock MS, Hunter DR, Butts CT, Morris M (2008a).
 A \pkg{statnet} Tutorial.
 \emph{Journal of Statistical Software}, 24(8).
 \url{http://www.jstatsoft.org/v24/i08/}.

 \item Handcock MS, Hunter DR, Butts CT, Goodreau SM, Krivitsky P, and Morris M (2012).
 \pkg{ergm}: A Package to Fit, Simulate and Diagnose
 Exponential-Family Models for Networks.
 Statnet Project, Seattle, WA.
 Version 3, \url{http://statnet.org}.

 \item Handcock MS, Hunter DR, Butts CT, Goodreau SM, Krivitsky P, Morris M (2012).
 \pkg{statnet}: Software Tools for the Statistical Modeling of
 Network Data.
 Statnet Project, Seattle, WA.
 Version 3, \url{http://statnet.org}.

 \item Hunter, D. R. and Handcock, M. S. (2006)
 Inference in curved exponential family models for networks,
 \emph{Journal of Computational and Graphical Statistics},
 15: 565-583
 
 \item Hunter DR, Handcock MS, Butts CT, Goodreau SM, Morris M (2008b).
 \pkg{ergm}: A Package to Fit, Simulate and Diagnose
 Exponential-Family Models for Networks.
 \emph{Journal of Statistical Software}, 24(3).
 \url{http://www.jstatsoft.org/v24/i03/}.
 
 \item Morris M, Handcock MS, Hunter DR (2008).
 Specification of Exponential-Family Random Graph Models:
 Terms and Computational Aspects.
 \emph{Journal of Statistical Software}, 24(4).
 \url{http://www.jstatsoft.org/v24/i04/}.
}  
}
\keyword{ package }
\keyword{ models }
