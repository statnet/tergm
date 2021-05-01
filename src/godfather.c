/*  File src/godfather.c in package tergm, part of the Statnet suite
 *  of packages for network analysis, https://statnet.org .
 *
 *  This software is distributed under the GPL-3 license.  It is free,
 *  open source, and has the attribution requirements (GPL Section 7) at
 *  https://statnet.org/attribution
 *
 *  Copyright 2008-2020 Statnet Commons
 */
#include "godfather.h"
#include "ergm_util.h"
#include "changestats_lasttoggle.h"

/*****************
 void godfather_wrapper

 ...we'll make them an offer (of toggles) they can't refuse.
 This function takes a list of toggles, each with a time stamp,
 then produces a matrix of changestats (with one row for each unique
 time stamp value) that result from performing all the toggles at
 each time step.  For instance, one might use this function to 
 find the changestats that result from starting from an empty network
 and then adding all of the edges to make up an observed network of interest.
*****************/
SEXP godfather_wrapper(SEXP stateR,
               SEXP total_toggles_arg,
               SEXP toggletimes_arg, 
               SEXP toggletails_arg,
               SEXP toggleheads_arg,
               SEXP start_time_arg,
               SEXP end_time_arg,
               SEXP verbose_arg){
  GetRNGstate();  /* R function enabling uniform RNG */
  ErgmState *s = ErgmStateInit(stateR, ERGM_STATE_NO_INIT_PROP);
  
  Network *nwp = s->nwp;
  Model *m = s->m;  
  
  /* Each ModelTerm in termarray has an aux_storage pointer,
     regardless of whether it asked for one; and the index of the
     lasttoggle auxiliary is in model$slots.extra.aux$system . Once we
     grab that, cast it to the lasttoggle data structure and extract
     the discord hashtable. */
  StoreTimeAndLasttoggle *dur_inf = (StoreTimeAndLasttoggle *)m->termarray->aux_storage[asInteger(getListElement(getListElement(m->R, "slots.extra.aux"), "system"))];

  int total_toggles = asInteger(total_toggles_arg);
  int *toggletimes = INTEGER(toggletimes_arg); 
  int *toggletails = INTEGER(toggletails_arg);
  int *toggleheads = INTEGER(toggleheads_arg);
  int start_time = asInteger(start_time_arg);
  int end_time = asInteger(end_time_arg);
  int verbose = asInteger(verbose_arg);
  
  /*********************
  changestats are modified in groups of m->n_stats, and they
  reflect the CHANGE in the values of the statistics from the
  original network.  Thus, when we begin, the initial 
  values of the first group of m->n_stats changestats should 
  all be zero
  *********************/
  
  SEXP changestatsRV = PROTECT(allocVector(REALSXP, (end_time - start_time + 1)*m->n_stats));
  double *changestats = REAL(changestatsRV);
  memset(changestats, 0, m->n_stats*sizeof(double));
  
  /* Now start obtaining change statistics */

  unsigned int pos = 0;
  // The reason it's = start_time but < end_time is that by the time t_stat arrives at end_time, the toggles for end_time will have already been applied.
  for(int t_stat = start_time; t_stat < end_time; t_stat++){
    changestats += m->n_stats;
    memcpy(changestats, changestats-m->n_stats, m->n_stats*sizeof(double));
    
    /* If the term has an extension, send it a "TICK" signal. */
    memset(m->workspace, 0, m->n_stats*sizeof(double)); /* Zero all change stats. */
    SEND_X_SIGNAL_INTO(nwp, m, NULL, m->workspace, TICK, NULL);
    /* Record network statistics for posterity. */
    addonto(changestats, m->workspace, m->n_stats);

    while(pos < total_toggles && toggletimes[pos]==t_stat+1){
      Vertex tail = toggletails[pos], head = toggleheads[pos];
      Rboolean edgestate = IS_OUTEDGE(tail, head);
      ChangeStats1(tail, head, nwp, m, edgestate);
      addonto(changestats, m->workspace, m->n_stats);
    
      ToggleKnownEdge(tail, head, nwp, edgestate);

      pos++;
    }
    
    MCMCDyn1Step_advance(s, dur_inf, changestats,
                         0, NULL, NULL, NULL, NULL, NULL,
                         verbose);
  }

  SEXP status = PROTECT(ScalarInteger(MCMCDyn_OK));
  const char *outn[] = {"status", "s", "state", ""};
  SEXP outl = PROTECT(mkNamed(VECSXP, outn));
  SET_VECTOR_ELT(outl, 0, status);
  SET_VECTOR_ELT(outl, 1, changestatsRV);
  
  /* record new generated network to pass back to R */
  if(asInteger(status) == MCMCDyn_OK){
    s->stats = REAL(changestatsRV) + (end_time - start_time)*m->n_stats;
    SET_VECTOR_ELT(outl, 2, ErgmStateRSave(s));
  }

  // save state for output as in MCMCDyn

  ErgmStateDestroy(s);  
  PutRNGstate();  /* Disable RNG before returning */
  UNPROTECT(3);  
  return outl;
}
