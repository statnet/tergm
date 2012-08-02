\name{NEWS}
\title{NEWS file for the \code{\link[=ergm-package]{ergm}} package}
\encoding{UTF-8}

% When a release is forked, copy the sections to a new release and
% delete unused sections.
\section{Changes in version 3.1}{

  \subsection{THIS IS AN INITIAL RELEASE}{
    This package consists of the dynamic network modeling code that has been
    split out of the \code{\link[=ergm-package]{ergm}} package.

    Changes listed in the following sections are relative to
    \code{\link[=ergm-package]{ergm}} version 3.0. 
  }
  
  \subsection{NEW FEATURES}{
    \itemize{
      
      \item Although fitting the EGMME for dissolution was possible
      before, it was impractical due to nonidentifiability: for example,
      one cannot fit both \code{\link{edges}} formation and
      \code{\link{edges}} dissolution with only one \code{\link{edges}}
      target statistic. Four new statistics have been added and
      documented that focus on targeting observed duration:
      \code{\link{mean.age}}, \code{\link{edgecov.mean.age}},
      \code{\link{edge.ages}}, and \code{\link{edgecov.ages}}. This
      allows jointly fitting formation and dissolution.

      \item In addition to the progress plot, EGMME fitting routines can
      now plot the estimated gradient matrix and the matrix of
      correlations among the target statistics.

      \item EGMME fitting can now take advantage of multiple CPUs,
      cores, or cluster nodes for faster and more robust fitting.
      
    } }
  
  \subsection{OTHER USER-VISIBLE CHANGES}{
    \itemize{
      \item The EGMME fitting algorithm has been vastly improved and is
      now a lot more adaptive and able to recover from problems. Many
      bugs have also been fixed.

    }
  }

  % \subsection{C-LEVEL FACILITIES}{
  %   \itemize{
  %     \item It is not clear that we need a subsection called "C-LEVEL FACILITIES".
  %   }
  % }

  % \subsection{UTILITIES}{
  %   \itemize{
  %     \item Likewise, "UTILITIES" is probably unnecessary.
  %   }
  % }

  % \subsection{INSTALLATION}{
  %   \itemize{
  %     \item ... as is "INSTALLATION".
  %   }
  % }

  % \subsection{PACKAGE INSTALLATION}{
  %   \itemize{
  %     \item ditto.
  %   }
  % }

  \subsection{DEPRECATED & DEFUNCT}{
    \itemize{
      \item A number of \code{\link{control.stergm}} parameters had been
      renamed and otherwise chaged.
    }
  }

  \subsection{BUG FIXES}{
    \itemize{
      
      \item The conditional MPLE (CMPLE) for formation and dissolution
      is now fit correctly. This also means that the starting values for
      the CMLE are much better.

      \item Bugs in EGMME code related to handling of bipartite networks
      have been fixed. 
    }
  }
}
