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
#include "tergm_model.h"
#include "ergm_state.h"

// TODO: This might be worth moving into a common "constants.h".
typedef enum MCMCDynStatus_enum {
  MCMCDyn_OK = 0,
  MCMCDyn_TOO_MANY_EDGES = 1,
  MCMCDyn_MH_FAILED = 2,
  MCMCDyn_TOO_MANY_CHANGES = 3
} MCMCDynStatus;

SEXP MCMCDyn_wrapper(ARGS_STATE, // ergm_state
                SEXP eta,      // double
                SEXP nsteps,   // integer
                SEXP min_MH_interval, // integer
                SEXP max_MH_interval, // integer
                SEXP MH_pval,  // double
                SEXP MH_interval_add, // double
                SEXP burnin, // integer
                SEXP interval, // integer
                SEXP collect, // integer (logical)
                SEXP maxedges, // integer
                SEXP maxchanges, // integer
                SEXP log_changes, // integer (logical)
                SEXP fVerbose);

MCMCDynStatus MCMCSampleDyn(ErgmState *s,
                double *eta,
                // Space for output.
                double *stats,
                int maxedges,
                int maxchanges,
                int log_changes,
                Vertex *difftime, Vertex *difftail, Vertex *diffhead, int *diffto,            
                // MCMC settings.
                unsigned int nsteps, unsigned int min_MH_interval, unsigned int max_MH_interval, double MH_pval, double MH_interval_add,
                unsigned int burnin, unsigned int interval, 
                // Verbosity.
                int fVerbose);

MCMCDynStatus MCMCDyn1Step(ErgmState *s,
                           double *eta,
                           // Space for output.
                           double *stats,
                           unsigned int maxchanges, Edge *nextdiffedge,
                           Vertex *difftime, Vertex *difftail, Vertex *diffhead, int *diffto,
                           // MCMC settings.
                           unsigned int min_MH_interval, unsigned int max_MH_interval, double MH_pval, double MH_interval_add,
                           // Verbosity.
                           int fVerbose);

MCMCDynStatus MCMCDyn1Step_advance(ErgmState *s,
                                   // Space for output.
                                   double *stats,
                                   unsigned int maxchanges, Edge *nextdiffedge,
                                   Vertex *difftime, Vertex *difftail, Vertex *diffhead, int *diffto,
                                   // Verbosity.
                                   int fVerbose);
#endif
