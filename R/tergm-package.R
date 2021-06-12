#' Fit, Simulate and Diagnose Dynamic Network Models derived from
#' Exponential-Family Random Graph Models
#'
#'
#' \code{\link[=tergm-package]{tergm}} is a collection of extensions to the
#' \code{\link[=ergm-package]{ergm}} package to fit, diagnose, and simulate
#' models for dynamic networks --- networks that evolve over time --- based on
#' exponential-family random graph models (ERGMs). For a list of functions type
#' \code{help(package='tergm')}
#'
#' When publishing results obtained using this package, please cite the
#' original authors as described in \code{citation(package="tergm")}.
#'
#' All programs derived from this package must cite it.
#'
#'
#' An exponential-family random graph model (ERGM) postulates an exponential
#' family over the sample space of networks of interest, and
#' \code{\link[=ergm-package]{ergm}} package implements a suite of tools for
#' modeling single networks using ERGMs.
#'
#' There have been a number of extensions of ERGMs for modeling the evolution
#' of networks, including the temporal ERGM (TERGM) of Hanneke et al. (2010)
#' and the separable termporal ERGM (STERGM) of Krivitsky and Handcock (2014).
#' The latter model allows familiar ERGM terms and statistics to be reused in a
#' dynamic context, interpreted in terms of formation and dissolution
#' (persistence) of ties. Krivitsky (2012) suggested a method for fitting
#' dynamic models when only a cross-sectional network is available, provided
#' some temporal information for it is available as well.
#'
#' This package aims to implement these and other ERGM-based models for network
#' evoluation. At this time, it implements, via the \code{\link{tergm}}
#' function, a general framework for modeling tie dynamics in temporal networks
#' with flexible model specification (including (S)TERGMs).  Estimation options
#' include a conditional MLE (CMLE) approach for fitting to a series of
#' networks and an Equilibrium Generalized Method of Moments Estimation (EGMME)
#' for fitting to a single network with temporal information. For further
#' development, see the referenced papers.
#'
#' If you previously used the [stergm()] function in this package, please
#' note that [stergm()] has been superceded by the new \code{\link{tergm}}
#' function, and is likely to be deprecated in future releases.  The
#' \code{dissolution} formula in [stergm()] maps to the new \code{Persist()}
#' operator in the [tergm()] function, **not** the \code{Diss()} operator.
#'
#' For detailed information on how to download and install the software, go to
#' the Statnet project website: \url{https://statnet.org}.  A tutorial, support
#' newsgroup, references and links to further resources are provided there.
#'
#' @name tergm-package
#' @docType package
#' @references \itemize{ \item Hanneke S, Fu W, and Xing EP (2010). Discrete
#' Temporal Models of Social Networks. \emph{Electronic Journal of Statistics},
#' 2010, 4, 585-605.
#' \doi{10.1214/09-EJS548}
#'
#' \item Krivitsky PN & Handcock MS (2014) A Separable Model for Dynamic
#' Networks. \emph{Journal of the Royal Statistical Society, Series B}, 76(1):
#' 29-46. \doi{10.1111/rssb.12014}
#'
#' \item Krivitsky, P.N. (2012). Modeling of Dynamic Networks based on
#' Egocentric Data with Durational Information. \emph{Pennsylvania State
#' University Department of Statistics Technical Report}, 2012(2012-01).
#' \url{https://web.archive.org/web/20170830053722/https://stat.psu.edu/research/technical-report-files/2012-technical-reports/TR1201A.pdf}
#'
#' \item Butts CT (2008).  \pkg{network}: A Package for Managing Relational
#' Data in .  \emph{Journal of Statistical Software}, 24(2).
#' \url{https://www.jstatsoft.org/v24/i02/}.
#'
#' \item Goodreau SM, Handcock MS, Hunter DR, Butts CT, Morris M (2008a).  A
#' \pkg{statnet} Tutorial.  \emph{Journal of Statistical Software}, 24(8).
#' \url{https://www.jstatsoft.org/v24/i08/}.
#'
#' \item Handcock MS, Hunter DR, Butts CT, Goodreau SM, Krivitsky P, and Morris
#' M (2012).  \pkg{ergm}: A Package to Fit, Simulate and Diagnose
#' Exponential-Family Models for Networks.  Statnet Project, Seattle, WA.
#' Version 3, \url{https://statnet.org}.
#'
#' \item Handcock MS, Hunter DR, Butts CT, Goodreau SM, Krivitsky P, Morris M
#' (2012).  \pkg{statnet}: Software Tools for the Statistical Modeling of
#' Network Data.  Statnet Project, Seattle, WA.  Version 3,
#' \url{https://statnet.org}.
#'
#' \item Hunter, D. R. and Handcock, M. S. (2006) Inference in curved
#' exponential family models for networks, \emph{Journal of Computational and
#' Graphical Statistics}, 15: 565-583
#'
#' \item Hunter DR, Handcock MS, Butts CT, Goodreau SM, Morris M (2008b).
#' \pkg{ergm}: A Package to Fit, Simulate and Diagnose Exponential-Family
#' Models for Networks.  \emph{Journal of Statistical Software}, 24(3).
#' \url{https://www.jstatsoft.org/v24/i03/}.
#'
#' \item Morris M, Handcock MS, Hunter DR (2008).  Specification of
#' Exponential-Family Random Graph Models: Terms and Computational Aspects.
#' \emph{Journal of Statistical Software}, 24(4).
#' \url{https://www.jstatsoft.org/v24/i04/}. }
#' @keywords package models
NULL