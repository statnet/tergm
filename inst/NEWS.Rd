%  File inst/NEWS.Rd in package tergm, part of the Statnet suite
%  of packages for network analysis, http://statnet.org .
%
%  This software is distributed under the GPL-3 license.  It is free,
%  open source, and has the attribution requirements (GPL Section 7) at
%  http://statnet.org/attribution
%
%  Copyright 2008-2017 Statnet Commons
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
\name{NEWS}
\title{NEWS file for the \code{\link[=tergm-package]{tergm}} package}
\encoding{UTF-8}

% When a release is forked, copy the following template to a new release
% and delete unused sections.

%% \section{Changes in version X.Y.Z}{
   %%  \subsection{NEW FEATURES}{
   %%   \itemize{
   %%   }
     
   %% }

  %% \subsection{BUG FIXES}{
  %%   \itemize{
  %%   }
  %% }

  % \subsection{OTHER USER-VISIBLE CHANGES}{
  %   \itemize{
  %   }
  % }

  % \subsection{DEPRECATED & DEFUNCT}{
  %   \itemize{
  %   }
  % }
  
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
%% }

\section{Changes in version 3.4.1}{
    \subsection{NEW FEATURES}{
      \itemize{
	\item Development has been migrated to GitHub.
      }
   }

  \subsection{BUG FIXES}{
    \itemize{
      \item Constraint violations in dynamic MCMC are now handled
      correctly.
      \item Native C routine registration has been added.
    }
  }

  \subsection{OTHER USER-VISIBLE CHANGES}{
    \itemize{
      By request from Sam Jennes,
      \code{\link{stergm.getMCMCsample.slave}} has been split out of
      \code{\link{stergm.getMCMCsample}} and exported.
    }
  }
}


\section{Changes in version 3.4.0}{
  \subsection{BUG FIXES}{
    \itemize{
      \item Miscellaneous bugfixes.

      \item The list of functions exported is now explicitly defined.

      \item \code{\link{robust.inverse}} has been replaced by
      \code{\link{ginv}}.
    }
  }
}


\section{Changes in version 3.3.1}{
    \subsection{NEW FEATURES}{
     \itemize{
       \item STERGM conditional MLE can now be fit to series of more
       than two bipartite networks. (Previously, for bipartite networks,
       only a series of two (i.e., one transition) could be fit.)
     }
     
   }

  \subsection{BUG FIXES}{
    \itemize{
      \item Miscellaneous CMLE bug fixes.

      \item Some low-priority warnings are no longer printed by
	default.

      \item A bug in formation TNT proposals has been fixed. The only
      analyses that may be affected are those on tiny (about 8 nodes) or
      extremely dense (i.e., high probability of a complete graph)
      networks.
    }
  }
}

\section{Changes in version 3.3}{
    \subsection{NEW FEATURES}{
     \itemize{
       \item A new control parameter, \code{SA.par.eff.pow}, now controls the degree to which the magnitude of the gradient afects the rescaling of the estimating functions (see help): higher means that "stronger" parameters effectively get lower gain, making the optimization more stable.
       \item Running with verbose=TRUE now prints \code{\link{gc}} output periodically. 
       \item Most durational functions can now take a binary option log=, which will result in them using log-ages instead of raw ages. (In a discrete-time model, there are no 0-age ties, so log transformation is safe.) 
       \item tergm.godfather example now tests for a few more terms, degrange ages in particular. 
       \item Implemented reference to \code{MH_listTNT} in applicable constraints in tergm as well.
       \item Durational statistics can now take \code{NULL} as their "role", disabling role testing.
     }
     
   }

  \subsection{BUG FIXES}{
    \itemize{
      \item  Running STERGM EGMME with parallel > 0 no longer causes it to eat up all the memory. 
      \item Namespace prefixing has been added to some functions called
      from external packages to deal with CRAN check warnings, replaced
      \code{\link{require}} and \code{\link{library}} calls with
      requireNamespace() 
      \item Updated \code{CITATION} file to deal with CRAN check warnings (\code{\link{packageDescription}} is no longer used).
      \item  Added S3 method registration in \code{NAMESPACE} for all S3
      functions.
      \item Handling of \code{\link{networkDynamic}} inputs has been improved.
    }
  }

  \subsection{OTHER USER-VISIBLE CHANGES}{
    \itemize{
      \item Changed CMLE MCMC default parameters to match those of ergm.
      \item Additions and corrections to documentation
      \item Most formerly undocumented 'internal' functions have been made internal by the \code{NAMESPACE}. As far as we know, this avoids functions in use by reverse-Depending packages.  These currently appear as commented items in the \code{NAMESPACE} file.
    }
  }
}

\section{Changes in version 3.2.4}{

  \subsection{BUG FIXES}{
    \itemize{
      \item{Some potential memory issues revealed by Valgrind have been fixed.}
    }
  }
}

\section{Changes in version 3.2}{
    \subsection{NEW FEATURES}{
     \itemize{
       \item Speed improvements. Internal data handling has been modified to avoid allocating (and passing around) a large matrix of tie ages unless the formula/model explicitly requires tie duration information. This increases speed of model initialization, which is especially helpful when simulating a network where the model must be re-initialized at each timestep (i.e. models with vital dynamics).
      \item{Simplified specification of offset values in formation model, when using stergm.EGMME with target statistics. Edges dissolution approximation is used to initialize the coefficients for estimation when there is an offset. See \link{http://statnet.csde.washington.edu/trac/wiki/tergmUsage} for examples.}
      \item{Offset terms in the target formula are automatically removed. This includes the case where the target is set equal to formation: \code{targets = "formation"}}
       \item The simulate.networkDynamic and related function now make use of networkDynamics persistant id (PID) functionality. If the network does not have any \code{\link[networkDynamic]{persistent.ids}}  defined for vertices, a vertex.pid will be attached in a vertex attribute  named \code{'tergm_pid'} to facilitate 'bookkeeping' between the networkDynamic argument and the simulated network time step. 
       \item{Parallel functionality: tergm can take user-created clusters as the control.stergm(parallel) control parameter. This is the recommended method for using tergm on a high-performance computing cluster. See \link{ergm-parallel}. Functionality is now implemented via the \code{parallel} package, direct dependence on \code{snow} package is removed}
       \item{Addition of bipartite formation and dissolution Metropolis-Hastings proposals}
       \item{Addition of a TNT dissolution proposal constraint}
       \item{Inclusion of a Stergm vignette}
       \item{-Inf offsets are recoded to a large negative number to be compatible with C code. The fitted coefficients will show the large negative number instead of the -Inf offset.}
     }
     
   }

  \subsection{BUG FIXES}{
    \itemize{
      \item{Namespace prefixing has been added to some functions called from external packages to deal with CRAN check warnings}
      \item{Bug fix where stergm constraints were not being passed from all the way down into san}
    }
  }

  \subsection{OTHER USER-VISIBLE CHANGES}{
    \itemize{
      \item Changes to \code{control.stergm} parameters:
      \itemize{
        \item \code{init.method} If \code{NULL} (the default), the initial values
    are computed using the edges dissolution approximation (Carnegie et al.) when appropriate. 
    If set to "zeros", the initial values are set to zeros.
        \item \code{SA.oh.memory = 100000} Absolute maximum number of data points per thread to store in the full optimization history.
      }
      \item  Changes to \code{simulate.stergm} parameters: addition of \code{duration.dependent} parameter: {Logical: Whether the model terms in formula or model are duration dependent. E.g., if a duration-dependent term  is used in estimation/simulation model, the probability of forming or dissolving a tie may dependent on the age the dyad status. If TRUE, the matrix of tie ages will be allocated.}
    }
  }

  % \subsection{NEW FEATURES}{
  %   \itemize{
  %   }
  % }
  
  % \subsection{DEPRECATED & DEFUNCT}{
  %   \itemize{
  %   }
  % }
  
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

}
\section{Changes in version 3.1.3}{

  \subsection{BUG FIXES}{
    \itemize{
      \item Some Metropolis-Hastings proposal functions would sometimes
      return incorrect acceptance probabilities when combined with the
      \code{\link[=ergm-constraints]{bd}} constraint. This has been
      fixed.

      \item In \code{\link{simulate.networkDynamic}}, an error in which
      vertices were queried and updated in the simulation has been
      fixed.
    }
  }

  \subsection{OTHER USER-VISIBLE CHANGES}{
    \itemize{
      \item \code{\link{simulate.networkDynamic}} and
      \code{\link{simulate.network}} functions now take an additional
      argument, \code{time.offset}. See help for those functions for
      details.
    }
  }

  % \subsection{NEW FEATURES}{
  %   \itemize{
  %   }
  % }
  
  % \subsection{DEPRECATED & DEFUNCT}{
  %   \itemize{
  %   }
  % }
  
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

}


\section{Changes in version 3.1.2}{

  \subsection{NEW FEATURES}{
    \itemize{
      \item The Tie-NonTie (TNT) proposal has been implemented for
      dissolution phase models. This should improve mixing and inference
      for these models.

      \item \code{stergm} and \code{simulate.stergm} now determine the
      number of Metropolis-Hastings steps per time step adaptively,
      stopping when the formation/dissolution process appears to have
      converged.
    }
  }
  
  \subsection{DEPRECATED & DEFUNCT}{
    \itemize{
      \item The previously deprecated start and end attr values attached 
      to networks returned by \code{\link{simulate.network}} and 
      \code{\link{simulate.networkDynamic}} have been removed and 
      replaced with observation spells recorded as a 
      \code{\link{net.obs.period}} network attribute.

      \item \code{MCMC.burnin} and \code{MCMC.interval} arguments to
      \code{control.stergm}, \code{control.simulate.network}, and
      \code{control.simulate.stergm} have been replaced by a different
      mechanism. See the help for the respective control functions for
      more information.
    }
  }

  \subsection{BUG FIXES}{
    \itemize{
      
      \item No longer generates deprecation warnings about start and end
       attrs when various internal functions use as.data.frame. Resolved
      by replacing start and end with a \code{\link{net.obs.period}}
      object.  

      \item A bug in the check for whether terms were amenable to being
      fit using STERGM CMLE for multiple transitions has been fixed.

      \item Some minor documentation typos have been fixed.
    }
  }
}

\section{Changes in version 3.1.1}{
  \itemize{
    \item Version 3.1.1 has been skipped to ease upgrading for those using a
    preview release.
  }
}

\section{Changes in version 3.1}{

  \subsection{THIS IS AN INITIAL RELEASE}{
    \itemize{
      \item This package consists of the dynamic network modeling code that has been
      split out of the \code{\link[=ergm-package]{ergm}} package.
      
      \item Changes listed in the following sections are relative to
      \code{\link[=ergm-package]{ergm}} version 3.0.
      }
  }
  
  \subsection{NEW FEATURES}{
    \itemize{
      
      \item Although fitting the EGMME for dissolution was possible
      before, it was impractical due to nonidentifiability: for example,
      one cannot fit both \code{\link{edges}} formation and
      \code{\link{edges}} dissolution with only one \code{\link{edges}}
      target statistic. Four new statistics have been added and
      documented that focus on targeting observed and hazard:
      \code{\link{mean.age}}, \code{\link{edgecov.mean.age}},
      \code{\link{degree.mean.age}}, \code{\link{degrange.mean.age}},
      \code{\link{edge.ages}}, \code{\link{edgecov.ages}}, and
      \code{\link{edges.ageinterval}}. This allow jointly fitting
      formation and dissolution.

      \item In addition to the progress plot, EGMME fitting routines can
      now plot the estimated gradient matrix and the matrix of
      correlations among the target statistics.

      \item EGMME fitting is now more adaptive when determining when to
      stop the optimization and return the result.

      \item EGMME fitting can now take advantage of multiple CPUs,
      cores, or cluster nodes for faster and more robust fitting.

      \item EGMME and CMLE fitting can now accommodate a
      \code{constraints} argument. However, note that the constraints
      apply to post-formation (\eqn{y^+=y^0\cup\y^1}) and
      post-dissolution (\eqn{y^-=y^0\cap\y^1}) networks, not to the
      final network (\eqn{y^1}). This may change in the future.

      \item CMLE can now be fit to more than two networks. (Not all ERGM
      terms and constraints can be used in this mode, however.)

      \item A new (sort of) function, \code{\link{tergm.godfather}} has
      been fixed and documented; it can be used to apply a specific set
      of changes to a network, returning statistics of interest as it
      evolves. In particular, it can be used to \dQuote{retrace} the evolution
      of a \code{\link{networkDynamic}}, calculating statistics at
      discrete time points along the way.

      \item \code{\link[=tergm-package]{tergm}} now implements a
      formula-based \code{\link[=summary.networkDynamic]{summary}}
      method for \code{\link{networkDynamic}} LHS, to compute dynamic
      network statistics at specified time points.

      \item Fitting CMLE to network series with transitioned-from
      networks having missing dyads is now possible, using automatic
      imputation. See \code{\link{impute.network.list}} and
      \code{\link{control.tergm}}.

      STERGM \code{\link[=simulate.stergm]{simulate}} can now be used
      for CMLE fits, takes a number of new arguments, and can be used to
      \dQuote{resume} a simulation from a \code{\link{networkDynamic}}
      object.

      \code{\link{gof}} methods have been implemented for CMLE fits.      
      
    }
  }
  
  \subsection{OTHER USER-VISIBLE CHANGES}{
    \itemize{
      \item \code{\link{stergm}} EGMME initial fitting code has been
      vastly improved.      
      
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
      renamed and otherwise changed.

      \item Argument \code{statsonly} to the
      \code{\link{simulate.stergm}} and related functions has been
      deprecated in favor of \code{output}.
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
