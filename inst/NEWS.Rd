%  File inst/NEWS.Rd in package tergm, part of the Statnet suite of packages
%  for network analysis, https://statnet.org .
%
%  This software is distributed under the GPL-3 license.  It is free, open
%  source, and has the attribution requirements (GPL Section 7) at
%  https://statnet.org/attribution .
%
%  Copyright 2008-2025 Statnet Commons
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
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


\section{Changes in version 4.2.2}{
  \subsection{C-LEVEL FACILITIES}{
    \itemize{
      \item{
        Update to work with \CRANpkg{ergm} 4.9 ABI.
      }
    }
  }
}


\section{Changes in version 4.2.1}{
  \subsection{NEW FEATURES}{
    \itemize{
      \item{
        In TERGM simulation, space for toggle history is now dynamically allocated; this is mainly visible in more efficient memory usage.
      }
     }
   }

  \subsection{BUG FIXES}{
    \itemize{
      \item{
        \code{Calloc()}, \code{Realloc()}, and \code{Free()} C calls now use their \code{R_}-prefixed counterparts.
      }
      \item{
        \code{\link{tergm}()} with CMLE for \code{\link[networkDynamic]{networkDynamic}} now stops with an informative error message if \code{times=} is not passed.
      }
      \item{
        Some corner cases in the discordant bounded-degree stratified TNT proposal have been fixed.
      }
      \item{
        General documentation fixes and improvements.
      }
    }
  }

  \subsection{OTHER USER-VISIBLE CHANGES}{
    \itemize{
      \item{
        Error messages in TERGM MCMC sampling are now more informative.
      }
    }
  }
}


\section{Changes in version 4.2.0}{
    \subsection{NEW FEATURES}{
      \itemize{
        \item Operators \code{Form()}, \code{Persist()}, \code{Diss()}, \code{Cross()}, and \code{Change()} can now specify time-varying network-level covariates using the syntax of \code{\link[ergm.multi:N-ergmTerm]{N}}.
        \item \code{\link{tergm}()} now supports a \code{basis=} argument like \code{\link{ergm}()} does.
        \item \pkg{tergm}'s functions now support \CRANpkg{networkLite} objects in place of \CRANpkg{network} objects.
     }
   }

  \subsection{BUG FIXES}{
    \itemize{
      \item \code{\link{tergm}()}'s term options are now propagated correctly.
      \item \pkg{tergm}'s \code{\link{simulate_formula.network}()} no longer propagates \code{dynamic} argument onto the \CRANpkg{ergm}'s back-end when \code{dynamic=FALSE}.
    }
  }

  \subsection{OTHER USER-VISIBLE CHANGES}{
    \itemize{
      \item \CRANpkg{tergm} now relies on \CRANpkg{ergm.multi} for much of its back-end.
    }
  }
}


\section{Changes in version 4.1.1}{
  \subsection{BUG FIXES}{
    \itemize{
      \item \code{networkDynamic} now handles dimnames more consistently.
    }
  }

  \subsection{OTHER USER-VISIBLE CHANGES}{
    \itemize{
      \item Documentation improvements.

      \item Terms \code{\link[=edgecov.mean.age-ergmTerm]{edgecov.mean.age}()} and \code{\link[=edgecov.ages-ergmTerm]{edgecov.ages}()} now handle their arguments more consistently with \code{\link[ergm:edgecov-ergmTerm]{edgecov}()} and stop with error if the network does not have the specified attribute.
    }
  }

  \subsection{DEPRECATED & DEFUNCT}{
    \itemize{
      \item \code{stergm()} and related facilities have been deprecated.

      \item Some vestigial files have been removed.

      \item Workaround storing \code{ergm::simulate_formula()} is no longer needed.
    }
  }
}



\section{Changes in version 4.1.0}{
  \subsection{NEW FEATURES}{
    \itemize{
      \item Term documentation has been ported to the \CRANpkg{ergm} 4.2 documentation API.

      \item New term operator \code{\link[=EdgeAges-ergmTerm]{EdgeAges}}.

      \item Term and proposal keyword/concept \code{durational} is registered by \CRANpkg{tergm}.
    }
  }

  \subsection{BUG FIXES}{
    \itemize{
      \item Fixes to handling of \code{\link{networkDynamic}} objects.

      \item Fixes to documentation URLs.
    }
  }

  \subsection{OTHER USER-VISIBLE CHANGES}{
    \itemize{
      \item Miscellaneous speed-ups.

      \item Documentation for some arguments and control parameters is now more standardized and detailed.

      \item Dissolution parameter estimation in the simple (EDApprox) case is a bit more accurate.
    }
  }

  \subsection{DEPRECATED & DEFUNCT}{
    \itemize{
      \item \code{\link{stergm}()} has been deprecated in favor of [tergm()], whose special case it is. It, and its associated functions may be removed in a future release.
    }
  }
}



\section{Changes in version 4.0.2}{
  \subsection{BUG FIXES}{
    \itemize{
      \item{
        Fixes to acommodate \code{\link[ergm:ergm-package]{ergm}}'s \code{BDStratTNT} API changes (which were not as backwards-compatible as hoped).
      }
      \item{
        Better error checking in \code{\link[tergm]{simulate}()} calls.
      }
      \item{
        Some memory leaks fixed.
      }
    }
  }

  \subsection{OTHER USER-VISIBLE CHANGES}{
    \itemize{
      \item{
        \code{\link{Cross}()} operator term now labels its coefficient names as such.
      }
    }
  }
}


\section{Changes in version 4.0.1}{
    \subsection{NEW FEATURES}{
      \itemize{
        \item The package has undergone a major rewrite, making use of operator terms introduced in \pkg{\link[ergm:ergm-packge]{ergm}} 4.0 to facilitate model specification. See the \dQuote{Conversions for \pkg{tergm} version 4.0} vignette for some examples.

        \item \code{\link{stergm}} has been superseded by \code{\link{tergm}}, which can fit non-separable discrete temporal ERGMs as well as separable.

        \item In temporal estimation and simulation, operator terms \code{\link{Form}()}, \code{\link{Persist}()}, and \code{\link{Diss}()} can be used to specify STERGMs, with \code{Form()} evaluating statistics on a post-formation network, \code{Persist()} on post-dissolution (i.e., persistence) network, and \code{Diss()} on the post-dissolution network, but with statistics negated for convenience. \code{\link{Change}()} and \code{\link{Cross}()} can be used to specify models for changes (i.e., a \dQuote{virtual} network that has edges where the previous network differs from the current one) and models for cross-section.

        \item \code{\link{NetSeries}()} constructor can be used to specify network series for conditional estimation and simulation (though \code{\link{tergm}} does it automatically for appropriate input).

        \item Constraints applied via \code{constraints=} argument for \code{\link{tergm}()} and \code{\link{simulate.tergm}()} are now applied to the sampled or estimated networks rather than to formation and dissolution phases separately.

        \item Control parameters in \pkg{\link[=tergm-package]{tergm}}'s functions can now be specified using \code{\link{snctrl}()}.
      }
   }

  \subsection{DEPRECATED & DEFUNCT}{
    \itemize{
      \item \code{\link{stergm}} and its family of functions has been deprecated.
    }
  }  
}

\section{Changes in version 3.7.0}{
  \subsection{BUG FIXES}{
    \itemize{
      \item Durational terms \code{edges.ageinterval} (\code{s_} function) and \code{edgecov.mean.age} (\code{s_} and \code{d_} functions) had bugs fixed.

      \item A bug involving assumptions based on the first entry of the \code{lasttoggle} vector has been corrected.

      \item A number of broken links has been fixed.
    }
  }

  \subsection{OTHER USER-VISIBLE CHANGES}{
    \itemize{
      \item The package no longer (incorrectly) requires \R 3.6.0 or the very-very latest version of \code{nlme}.

      \item Progress reports are now printed to \code{stderr} rather than \code{stdout}.
    }
  }
}


\section{Changes in version 3.6.0}{
    \subsection{NEW FEATURES}{
     \itemize{
	
       \item Compatibility with \code{\link[=ergm:ergm-package]{ergm}} 3.10 API changes.
       
       \item Optimizations to make use of \code{\link[=ergm:ergm-package]{ergm}} 3.10 APIs and parallel facilities.

       \item Improved utilization of \code{\link[ergm]{san}}.

       \item \code{\link{stergm}}'s \code{eval.loglik=} argument's default can now be controlled by a \code{\link{options}}-style option.

     }
     
   }

  \subsection{BUG FIXES}{
    \itemize{

      \item Miscellaneous compiler warnings fixed.
      
    }
  }
}


\section{Changes in version 3.5.2}{
    \subsection{NEW FEATURES}{
      \itemize{
	
	\item Compatibility with \code{\link[=ergm:ergm-package]{ergm}}
	3.9 API changes.
	
     }
     
   }

  \subsection{BUG FIXES}{
    \itemize{
      
      \item Documentation fixes.

      \item A number of methods are no longer called directly but
      through generics.
      
    }
  }

  \subsection{DEPRECATED & DEFUNCT}{
    \itemize{
      
      \item \code{\link{stergm.getMCMCsample}} and
      \code{\link{stergm.mcmcslave}} have been deprecatedin favor of
      \code{\link{stergm_MCMC_sample}} and
      \code{\link{stergm_MCMC_slave}}, respectively.

    }
  }
}


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
      \item By request from Sam Jennes,
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
      \item{Simplified specification of offset values in formation model, when using stergm.EGMME with target statistics. Edges dissolution approximation is used to initialize the coefficients for estimation when there is an offset. See \link{https://statnet.org/Workshops/tergm_tutorial.html} for examples.}
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
      \item{
        Changes to \code{control.stergm} parameters:
        \itemize{
          \item \code{init.method} If \code{NULL} (the default), the initial values
          are computed using the edges dissolution approximation (Carnegie et al.) when appropriate.
          If set to "zeros", the initial values are set to zeros.
          \item \code{SA.oh.memory = 100000} Absolute maximum number of data points per thread to store in the full optimization history.
        }
      }
      \item  Changes to \code{simulate.stergm} parameters: addition of \code{duration.dependent} parameter to control whether the model terms in formula or model are duration dependent. E.g., if a duration-dependent term  is used in estimation/simulation model, the probability of forming or dissolving a tie may dependent on the age the dyad status.
    }
  }
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
