/*  File src/godfather.h in package tergm, part of the Statnet suite
 *  of packages for network analysis, https://statnet.org .
 *
 *  This software is distributed under the GPL-3 license.  It is free,
 *  open source, and has the attribution requirements (GPL Section 7) at
 *  https://statnet.org/attribution
 *
 *  Copyright 2008-2019 Statnet Commons
 */
#ifndef godfather_H 
#define godfather_H

#include <R.h>
#include "ergm_model.h"
#include "MCMCDyn.h"

/* Function prototypes */
SEXP godfather_wrapper(SEXP stateR,
               SEXP total_toggles_arg,
               SEXP toggletimes_arg, 
               SEXP toggletails_arg,
               SEXP toggleheads_arg,
               SEXP start_time_arg,
               SEXP end_time_arg,
               SEXP fVerbose_arg);
      
#endif

