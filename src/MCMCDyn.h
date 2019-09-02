/*  File src/MCMCDyn.h in package tergm, part of the Statnet suite
 *  of packages for network analysis, https://statnet.org .
 *
 *  This software is distributed under the GPL-3 license.  It is free,
 *  open source, and has the attribution requirements (GPL Section 7) at
 *  https://statnet.org/attribution
 *
 *  Copyright 2008-2019 Statnet Commons
 */
#ifndef MCMCDYN_H
#define MCMCDYN_H
#include<string.h>
#include "ergm_edgetree.h"
#include "ergm_MHproposal.h"
#include "ergm_changestat.h"
#include "ergm_model.h"

// TODO: This might be worth moving into a common "constants.h".
typedef enum MCMCDynStatus_enum {
  MCMCDyn_OK = 0,
  MCMCDyn_TOO_MANY_EDGES = 1,
  MCMCDyn_MH_FAILED = 2,
  MCMCDyn_TOO_MANY_CHANGES = 3
} MCMCDynStatus;


void MCMCDyn_init_common(int *tails, int *heads, int time, int *lasttoggle, int n_edges,
			 int n_nodes, int dflag, int bipartite, Network **nwp,
			 
			 int nterms, char *funnames, char *sonames, double *inputs, Model **m,
			 
			 int *attribs, int *maxout, int *maxin, int *minout,
			 int *minin, int condAllDegExact, int attriblength,
			 
			 char *MHProposaltype, char *MHProposalpackage, MHProposal **MH,
                         StoreDyadSet **discord,
			 int fVerbose);


void MCMCDyn_finish_common(Network *nwp,
			   Model *m,
			   MHProposal *MH,
                           StoreDyadSet *discord);

MCMCDynStatus MCMCSampleDyn(// Observed and discordant network.
			    Network *nwp, StoreDyadSet *discord,
			    // terms and proposals.
			    Model *m, MHProposal *MH,
			    double *eta,
			    // Space for output.
			    double *stats,
			    Edge maxedges,
			    Edge maxchanges,
			    int log_changes,
			    Vertex *difftime, Vertex *difftail, Vertex *diffhead, int *diffto,		    
			    // MCMC settings.
			    unsigned int nsteps, unsigned int min_MH_interval, unsigned int max_MH_interval, double MH_pval, double MH_interval_add,
			    unsigned int burnin, unsigned int interval, 
			    // Verbosity.
			    int fVerbose);

MCMCDynStatus MCMCDyn1Step(// Observed and discordant network.
                           Network *nwp, StoreDyadSet *discord,
		  // terms and proposals.
		  Model *m, MHProposal *MH, double *eta,
		  // Space for output.
		  unsigned log_changes,
		  double *stats,
		  unsigned int maxchanges, Edge *nextdiffedge,
		  Vertex *difftime, Vertex *difftail, Vertex *diffhead, int *diffto,
		  // MCMC settings.
		  unsigned int min_MH_interval, unsigned int max_MH_interval, double MH_pval, double MH_interval_add,
		  // Verbosity.
		  int fVerbose);

#endif
