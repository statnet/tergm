/*  File src/DynSA.h in package tergm, part of the Statnet suite
 *  of packages for network analysis, https://statnet.org .
 *
 *  This software is distributed under the GPL-3 license.  It is free,
 *  open source, and has the attribution requirements (GPL Section 7) at
 *  https://statnet.org/attribution
 *
 *  Copyright 2008-2019 Statnet Commons
 */
#ifndef DYNSA_H
#define DYNSA_H

#include "MCMCDyn.h"
#include "ergm_MHproposal.h"
                 
MCMCDynStatus MCMCDynSArun(ErgmState *s,
                  StoreDyadMapInt *discord,
                  int nstatsmonitor,
                  // Model fitting.
                  double *eta, 
                  double *inputdev, // DEViation of the current network's targeted statistics from the target statistics.
                  int runlength,
                  double *WinvGradient, double *jitter, double *dejitter, 
                  double *dev_guard, double *par_guard,
                  
                  // Space for output.
                  int maxedges, int maxchanges,
                  double *opt_history,
                  // MCMC settings.
                  unsigned int SA_burnin, unsigned int SA_interval, unsigned int min_MH_interval, unsigned int max_MH_interval, double MH_pval, double MH_interval_add,
                  // Verbosity.
                  int verbose);
#endif
