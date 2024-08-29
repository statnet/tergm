/*  File src/MCMCDyn.h in package tergm, part of the
 *  Statnet suite of packages for network analysis, https://statnet.org .
 *
 *  This software is distributed under the GPL-3 license.  It is free,
 *  open source, and has the attribution requirements (GPL Section 7) at
 *  https://statnet.org/attribution .
 *
 *  Copyright 2008-2023 Statnet Commons
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
#include "changestats_lasttoggle.h"
/*#include "diff_vect.h"*/

#include "ergm_kvec.h"
typedef kvec_t(int) kvint;

// TODO: This might be worth moving into a common "constants.h".
typedef enum MCMCDynStatus_enum {
  MCMCDyn_OK = 0,
  MCMCDyn_TOO_MANY_EDGES = 1,
  MCMCDyn_MH_FAILED = 2,
  MCMCDyn_TOO_MANY_CHANGES = 3
} MCMCDynStatus;

MCMCDynStatus MCMCSampleDyn(ErgmState *s,
                StoreTimeAndLasttoggle *dur_inf,
                double *eta,
                // Space for output.
                double *stats,
                int maxedges,
                int maxchanges,
                int log_changes,
                kvint *difftime, kvint *difftail, kvint *diffhead, kvint *diffto,
                // MCMC settings.
                unsigned int nsteps, unsigned int min_MH_interval, unsigned int max_MH_interval, double MH_pval, double MH_interval_add,
                unsigned int burnin, unsigned int interval,
                // Verbosity.
                int verbose);

MCMCDynStatus MCMCDyn1Step(ErgmState *s,
                           StoreTimeAndLasttoggle *dur_inf,
                           double *eta,
                           // Space for output.
                           double *stats,
                           unsigned int maxchanges, Edge *nextdiffedge,
                           kvint *difftime, kvint *difftail, kvint *diffhead, kvint *diffto,
                           // MCMC settings.
                           unsigned int min_MH_interval, unsigned int max_MH_interval, double MH_pval, double MH_interval_add,
                           // Verbosity.
                           int verbose);

MCMCDynStatus MCMCDyn1Step_advance(ErgmState *s,
                                   StoreTimeAndLasttoggle *dur_inf,
                                   // Space for output.
                                   double *stats,
                                   unsigned int maxchanges, Edge *nextdiffedge,
                                   kvint *difftime, kvint *difftail, kvint *diffhead, kvint *diffto,
                                   // Verbosity.
                                   int verbose);
#endif
