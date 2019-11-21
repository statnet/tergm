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

void MCMCDynSArun_wrapper(// Observed network.
                 int *tails, int *heads, int *time, int *lasttoggle, int *n_edges,
                 int *n_nodes, int *dflag, int *bipartite, 
                 // Formation terms and proposals.
                 int *nterms, char **funnames, char **sonames,
                 char **MHProposaltype, char **MHProposalpackage,
                 double *inputs, int *nstatmonitor,
                 // Parameter fitting.
                 double *eta0,
                 double *init_dev,
                 int *runlength,
                 double *WinvGradient,
                 double *jitter, double *dejitter,
                 double *dev_guard,
                 double *par_guard,
                 // Degree bounds.
                 int *attribs, int *maxout, int *maxin, int *minout,
                 int *minin, int *condAllDegExact, int *attriblength,
                 // MCMC settings.
                 int *SA_burnin, int *SA_interval, int *min_MH_interval, int *max_MH_interval, double *MH_pval, double *MH_interval_add,
                 // Space for output.
                 int *maxedges, int *maxchanges,
                 int *newnetworktail, int *newnetworkhead, 
                 double *opt_history,
                 // Verbosity.
                 int *fVerbose,
                 int *status);
MCMCDynStatus MCMCDynSArun(// Observed and discordant network.
                  Network *nwp, StoreDyadMapInt *discord,
                  // Formation terms and proposals.
                  Model *m, MHProposal *MHp,
                  int nstatsmonitor,
                  // Model fitting.
                  double *eta, 
                  double *inputdev, // DEViation of the current network's targeted statistics from the target statistics.
                  int runlength,
                  double *WinvGradient, double *jitter, double *dejitter, 
                  double *dev_guard, double *par_guard,
                  
                  // Space for output.
                  Edge maxedges, Edge maxchanges,
                  Vertex *difftime, Vertex *difftail, Vertex *diffhead,
                  double *opt_history,
                  // MCMC settings.
                  unsigned int SA_burnin, unsigned int SA_interval, unsigned int min_MH_interval, unsigned int max_MH_interval, double MH_pval, double MH_interval_add,
                  // Verbosity.
                  int fVerbose);
#endif
