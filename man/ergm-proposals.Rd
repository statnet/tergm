%  File man/ergm-proposals.Rd in package tergm, part of the
%  Statnet suite of packages for network analysis, https://statnet.org .
%
%  This software is distributed under the GPL-3 license.  It is free,
%  open source, and has the attribution requirements (GPL Section 7) at
%  https://statnet.org/attribution .
%
%  Copyright 2008-2021 Statnet Commons
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
\name{ergm-proposals}
\alias{ergm-proposals}
\alias{InitErgmProposal}
\alias{InitErgmProposal.staticDiscordTNT}
\alias{InitErgmProposal.discordTNT}
\alias{InitErgmProposal.discordBDStratTNT}

\docType{methods}

\title{Metropolis-Hastings Proposal Methods for TERGM MCMC}

\description{
  \code{\link{tergm}} uses a Metropolis-Hastings (MH) algorithm provided by 
  \code{\link{ergm}} to control the behavior of the Markov Chain Monte Carlo 
  (MCMC) for sampling networks.  The MCMC chain is intended to step around the
  sample space of possible networks, selecting a network at regular intervals 
  to evaluate the statistics in the model.  For each MCMC step, \eqn{n} 
  (\eqn{n=1} in the simple case) toggles are proposed to change the dyad(s) 
  to the opposite value. The probability of accepting the proposed change is
  determined by the MH acceptance ratio.  The role of the different MH methods
  implemented in \code{\link{tergm}} is to vary how the sets of dyads are 
  selected for toggle proposals.  This is used in some cases to improve the 
  performance (speed and mixing) of the algorithm, and in other cases to 
  constrain the sample space.

  The user may specify proposals directly if need be, but is encouraged to 
  utilize the hints and constraints UI instead.
}
  
\details{
  \subsection{MH proposal methods implemented in the \code{\link[=tergm-package]{tergm}} package}{
    \itemize{
      \item{\code{InitErgmProposal.staticDiscordTNT}} {
      
        A version of \code{\link[ergm:ergm-proposals]{TNT}} appropriate for CMLE fitting, 
        with proposals stratified both by discordance status and edge status.  The argument
        \code{ref} specifies the data relative to which discordance status is defined.
      }

      \item{\code{InitErgmProposal.discordTNT}} {
      
        A temporal version of \code{\link[ergm:ergm-proposals]{TNT}}, with approximately
        \code{discordance_fraction} of proposed toggles being made on the set of discordant dyads,
        and approximately \code{1 - discordance_fraction} of proposed toggles being TNT proposals from
        the network.  The value of \code{discordance_fraction} can be set by the user as a proposal argument,
        and defaults to \code{0.5}.
      }
                  
      \item{\code{InitErgmProposal.discordBDStratTNT}} {

        A temporal version of \code{\link[ergm:ergm-proposals]{BDStratTNT}}.  Within each 
        strat mixing type, approximately 50\% of proposed toggles are made on
        discordant dyads, and approximately 50\% of proposed toggles are
        TNT proposals from the network, all subject to the bounded degree
        and mixing type constraints.  The degree bound constraint is imposed 
        on the instantaneous network state 
        (rather than the temporal operator networks).
        
        Proposal arguments are the same as for \code{\link[ergm:ergm-proposals]{BDStratTNT}},
        and should be passed in via the \code{\link[ergm:ergm-constraints]{bd}} and
        \code{\link[ergm:ergm-constraints]{blocks}} constraints and 
        the \code{\link[ergm:ergm-hints]{Strat}} hint.
      }    
    }  
  }
}

\seealso{\code{\link[=tergm-package]{tergm}} package, \code{\link{ergm}}, 
         \code{\link{ergm-constraints}}, and ergm's 
         \code{\link[ergm]{ergm_proposal}} help page}


\references{
  
  Goodreau SM, Handcock MS, Hunter DR, Butts CT, Morris M (2008a).
  A \pkg{statnet} Tutorial.
  \emph{Journal of Statistical Software}, 24(8).
  \url{https://www.jstatsoft.org/v24/i08/}.
  
  Hunter, D. R. and Handcock, M. S. (2006)
  \emph{Inference in curved exponential family models for networks},
  Journal of Computational and Graphical Statistics.
  
  Hunter DR, Handcock MS, Butts CT, Goodreau SM, Morris M (2008b).
  \pkg{ergm}: A Package to Fit, Simulate and Diagnose
  Exponential-Family Models for Networks.
  \emph{Journal of Statistical Software}, 24(3).
  \url{https://www.jstatsoft.org/v24/i03/}.
  
  Krivitsky PN (2012). Exponential-Family Random Graph Models for Valued
  Networks. \emph{Electronic Journal of Statistics}, 2012, 6,
  1100-1128. \doi{10.1214/12-EJS696}
  
  Morris M, Handcock MS, Hunter DR (2008).
  Specification of Exponential-Family Random Graph Models:
    Terms and Computational Aspects.
  \emph{Journal of Statistical Software}, 24(4).
  \url{https://www.jstatsoft.org/v24/i04/}.
  
}

\keyword{models}
\keyword{internal}

