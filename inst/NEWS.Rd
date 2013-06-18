%  File inst/NEWS.Rd in package tergm, part of the Statnet suite
%  of packages for network analysis, http://statnet.org .
%
%  This software is distributed under the GPL-3 license.  It is free,
%  open source, and has the attribution requirements (GPL Section 7) at
%  http://statnet.org/attribution
%
%  Copyright 2003-2013 Statnet Commons
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
\name{NEWS}
\title{NEWS file for the \code{\link[=tergm-package]{tergm}} package}
\encoding{UTF-8}

% When a release is forked, copy the sections to a new release and
% delete unused sections.
\section{Changes in version 3.1.1}{

  

  \subsection{DEPRECATED & DEFUNCT}{
    \itemize{
      \item The previously deprecated start and end attr valeus attached 
           to networks returned by \code{\link{simulate.network}} and 
           \code{\link{simulate.networkDynamic}} have been removed and 
           replaced with observation spells recorded as a 
           \code{\link{net.obs.period}} network attribute. 
    }
  }

  \subsection{BUG FIXES}{
    \itemize{
      
      \item No longer generates deprecation warnings about start and end
       attrs when various internal functions use as.data.frame. Resolved
      by replacing start and end with a \code{\link{net.obs.period}} object.  
    }
  }
}
