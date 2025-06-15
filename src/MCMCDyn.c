/*  File src/MCMCDyn.c in package tergm, part of the Statnet suite of packages
 *  for network analysis, https://statnet.org .
 *
 *  This software is distributed under the GPL-3 license.  It is free, open
 *  source, and has the attribution requirements (GPL Section 7) at
 *  https://statnet.org/attribution .
 *
 *  Copyright 2008-2025 Statnet Commons
 */
#include "MCMCDyn.h"
#include "ergm_util.h"
#include "changestats_lasttoggle.h"
/*****************
 Note on undirected networks:  For j<k, edge {j,k} should be stored
 as (j,k) rather than (k,j).  In other words, only directed networks
 should have (k,j) with k>j.
*****************/

/*****************
 void MCMCDyn_wrapper

 Wrapper for a call from R.
*****************/
SEXP MCMCDyn_wrapper(SEXP stateR, // ergm_state
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
                SEXP verbose){  // integer
  GetRNGstate();  /* R function enabling uniform RNG */
  ErgmState *s = ErgmStateInit(stateR, 0);

  Model *m = s->m;
  MHProposal *MHp = s->MHp;

  /* Each ModelTerm in termarray has an aux_storage pointer,
     regardless of whether it asked for one; and the index of the
     lasttoggle auxiliary is in model$slots.extra.aux$system . Once we
     grab that, cast it to the lasttoggle data structure and extract
     the discord hashtable. */
  StoreTimeAndLasttoggle *dur_inf = (StoreTimeAndLasttoggle *)m->termarray->aux_storage[asInteger(getListElement(getListElement(m->R, "slots.extra.aux"), "system"))];

  SEXP sample = PROTECT(allocVector(REALSXP, (asInteger(nsteps) + 1)*m->n_stats));
  memset(REAL(sample), 0, (asInteger(nsteps) + 1)*m->n_stats*sizeof(double));
  memcpy(REAL(sample), s->stats, m->n_stats*sizeof(double));
  
  kvint difftime, difftail, diffhead, diffto;
  kv_init(difftime);
  kv_init(difftail);
  kv_init(diffhead);
  kv_init(diffto);

  // pre-allocate the 0th element for size
  kv_push(int, difftime, 0);
  kv_push(int, difftail, 0);
  kv_push(int, diffhead, 0);
  kv_push(int, diffto, 0);

  SEXP status;
  if(MHp) status = PROTECT(ScalarInteger(MCMCSampleDyn(s,
              dur_inf,
              REAL(eta),
              asInteger(collect)?REAL(sample):NULL, asInteger(maxedges), asInteger(maxchanges), asInteger(log_changes), &difftime, &difftail, &diffhead, &diffto,
              asInteger(nsteps), asInteger(min_MH_interval), asInteger(max_MH_interval), asReal(MH_pval), asReal(MH_interval_add), asInteger(burnin), asInteger(interval),
              asInteger(verbose))));
  else status = PROTECT(ScalarInteger(MCMCDyn_MH_FAILED));
   
  const char *outn[] = {"status", "s", "state", "diffnwtime", "diffnwtails", "diffnwheads", "diffnwdirs", ""};
  SEXP outl = PROTECT(mkNamed(VECSXP, outn));
  SET_VECTOR_ELT(outl, 0, status);
  SET_VECTOR_ELT(outl, 1, sample);
  
  /* record new generated network to pass back to R */
  if(asInteger(status) == MCMCDyn_OK){
    s->stats = REAL(sample) + asInteger(nsteps)*m->n_stats;
    SET_VECTOR_ELT(outl, 2, ErgmStateRSave(s));
  }
  
  SET_VECTOR_ELT(outl, 3, PROTECT(kvint_to_SEXP(difftime)));
  SET_VECTOR_ELT(outl, 4, PROTECT(kvint_to_SEXP(difftail)));
  SET_VECTOR_ELT(outl, 5, PROTECT(kvint_to_SEXP(diffhead)));
  SET_VECTOR_ELT(outl, 6, PROTECT(kvint_to_SEXP(diffto)));

  kv_destroy(difftime);
  kv_destroy(difftail);
  kv_destroy(diffhead);
  kv_destroy(diffto);

  ErgmStateDestroy(s);  
  PutRNGstate();  /* Disable RNG before returning */
  UNPROTECT(7);
  return outl;
}

/*********************
 MCMCDynStatus MCMCSampleDyn

 Using the parameters contained in the array eta, obtain the
 network statistics for a sample of size nsteps.  burnin is the
 initial number of Markov chain steps before sampling anything
 and interval is the number of MC steps between successive 
 networks in the sample.  Put all the sampled statistics into
 the statistics array. 
*********************/
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
                int verbose){
  Network *nwp = s->nwp;
  Model *m = s->m;

  int i, j;
  Edge nextdiffedge=1;


  /*if (verbose)
    Rprintf("Total m->n_stats is %i; total nsteps is %d\n",
    m->n_stats,nsteps);*/
  
  
  /* Burn in step. */

  for(i=0;i<burnin;i++){
    R_CheckUserInterrupt();
    MCMCDynStatus status = MCMCDyn1Step(s,
                    dur_inf,
                    eta,
                    stats,
                    maxchanges, log_changes ? &nextdiffedge : NULL, difftime, difftail, diffhead, diffto,
                    min_MH_interval, max_MH_interval, MH_pval, MH_interval_add, verbose);
    // Check that we didn't run out of log space.
    if(status==MCMCDyn_TOO_MANY_CHANGES)
      return MCMCDyn_TOO_MANY_CHANGES;
  
    // If we need to return a network, then stop right there, since the network is too big to return, so stop early.
    if(maxedges!=0 && EDGECOUNT(nwp) >= maxedges-1)
      return MCMCDyn_TOO_MANY_EDGES;
  }
  
  //Rprintf("MCMCSampleDyn post burnin numdissolve %d\n", *numdissolve);
  
  if (verbose){
    Rprintf("Returned from STERGM burnin\n");
  }
  
  /* Now sample networks */
  for (i=0; i < nsteps; i++){
    /* Set current vector of stats equal to previous vector */
    if(stats){
      for (j=0; j<m->n_stats; j++){
    stats[j+m->n_stats] = stats[j];
      }
      stats += m->n_stats;
    }

    /* This then adds the change statistics to these values */
    for(j=0;j<interval;j++){
      R_CheckUserInterrupt();
      MCMCDynStatus status = MCMCDyn1Step(s,
                      dur_inf,
                      eta,
                      stats,
                      maxchanges, log_changes ? &nextdiffedge : NULL, difftime, difftail, diffhead, diffto,
                      min_MH_interval, max_MH_interval, MH_pval, MH_interval_add, verbose);
      
      // Check that we didn't run out of log space.
      if(status==MCMCDyn_TOO_MANY_CHANGES)
        return MCMCDyn_TOO_MANY_CHANGES;
      
      // If we need to return a network, then stop right there, since the network is too big to return, so stop early.
      if(maxedges!=0 && EDGECOUNT(nwp) >= maxedges-1)
    return MCMCDyn_TOO_MANY_EDGES;
    }
    
    //Rprintf("MCMCSampleDyn loop numdissolve %d\n", *numdissolve);
    if (verbose){
      if( ((3*i) % nsteps)<3 && nsteps > 500){
        Rprintf("Advanced %d time steps.\n", i);
      }
    }
  }

  if(log_changes) {
      kv_A(*difftime, 0) = nextdiffedge - 1;
      kv_A(*difftail, 0) = nextdiffedge - 1;
      kv_A(*diffhead, 0) = nextdiffedge - 1;
      kv_A(*diffto, 0) = nextdiffedge - 1;
  }
  return MCMCDyn_OK;
}


/* The following outlines what happens during each time step:

   1. All statistics with x_ functions are sent a TICK signal. This triggers _lasttoggle auxiliary to increment the timer.
   2. Any updates to statistics are recorded.
   3. A toggle proposal is made.
   4. Change statistics are calculated.
   5. If proposal is rejected, continue to 9.
   7. Relevant dyad is toggled in nwp. This triggers _lasttoggle to do the following:
      * If it hasn't been toggled during this time step, its lasttoggle information, if any, is saved in discord, and its lasttoggle information is updated.
      * If it has been toggled during this time step, its lasttoggle information, if any, is restored from discord. Discord's backup is deleted.
   8. Continue from 4 until sufficiently long run.
   9. Toggles in discord are logged.
   10. All statistics with x_ functions are sent a TOCK signal. This triggers _lasttoggle auxiliary to clear discord.
   11. Any updates to statistics are recorded.

 */

/*********************
 void MCMCDyn1Step

 Simulate evolution of a dynamic network for 1 step.
*********************/
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
                           int verbose){                              
  Network *nwp = s->nwp;
  Model *m = s->m;
  MHProposal *MHp = s->MHp;
  StoreDyadMapInt *discord = dur_inf->discord;

  /* If the term has an extension, send it a "TICK" signal. */
  memset(m->workspace, 0, m->n_stats*sizeof(double)); /* Zero all change stats. */
  SEND_X_SIGNAL_INTO(nwp, m, MHp, m->workspace, TICK, NULL);
  /* Record network statistics for posterity. */
  if(stats) addonto(stats, m->workspace, m->n_stats);

  /* Run the process. */
  
  double cutoff;
  double 
    si = 0, // sum of increments
    si2 = 0, // sum of squared increments
    sw = 0, // sum of weights 
    sw2 = 0 // sum of squared weights
    ;
  double sdecay = 1 - 1.0/min_MH_interval;
  
  unsigned int step=0; // So that we could print out the number of steps later.
  for(unsigned int finished = 0, extrasteps = 0; step < max_MH_interval && finished < extrasteps+1; step++) {
    unsigned int prev_discord = kh_size(discord);
    
    MHp->logratio = 0;
    (*(MHp->p_func))(MHp, nwp); /* Call MHp function to propose toggles */
    //      Rprintf("Back from proposal; step=%d\n",step);

    // Proposal failed.
    if(MHp->toggletail[0]==MH_FAILED){
      switch(MHp->togglehead[0]){
      case MH_UNRECOVERABLE:
    error("Something very bad happened during proposal. Memory has not been deallocated, so restart R soon.");
    
      case MH_IMPOSSIBLE:
    Rprintf("MH MHProposal function encountered a configuration from which no toggle(s) can be proposed.\n");
    return MCMCDyn_MH_FAILED;
    
      case MH_UNSUCCESSFUL:
      case MH_CONSTRAINT:
    continue;
      }
    }

    ChangeStats(MHp->ntoggles, MHp->toggletail, MHp->togglehead, nwp, m);
    
    //  Rprintf("change stats:"); 
    /* Calculate inner product */
    double ip = dotprod(eta, m->workspace, m->n_stats);
      //  Rprintf("%f ", m->workspace[i]); 
    //}
    //  Rprintf("\n ip %f dedges %f\n", ip, m->workspace[0]); 
    /* The logic is to set exp(cutoff) = exp(ip) * qratio ,
       then let the MHp probability equal min{exp(cutoff), 1.0}.
       But we'll do it in log space instead.  */
    cutoff = ip + MHp->logratio;
    
    /* if we accept the proposed network */
    if (cutoff >= 0.0 || log(unif_rand()) < cutoff) { 
      /* Hold off updating timesteamps until the changes are committed,
         which doesn't happen until later. */
      for (unsigned int i=0; i < MHp->ntoggles; i++){
        ToggleEdge(MHp->toggletail[i], MHp->togglehead[i], nwp);
      }
      /* Record network statistics for posterity. */
      if(stats) {
        addonto(stats, m->workspace, m->n_stats);
      }
    }

    int i = kh_size(discord) - prev_discord;
    sw *= sdecay; si *= sdecay;
    sw++; si += i;
    sw2 *= sdecay*sdecay; si2 *= sdecay;
    sw2++; si2 += i*i;
         
    if(step >= min_MH_interval && !finished) { 
      // Now, perform the test:
      double mi = (double)si / sw, mi2 = (double)si2 / sw;
      
      double vi = mi2 - mi*mi;
      double zi = mi / sqrt(vi * sw2/(sw*sw)); // denom = sqrt(sum(w^2 * v)/sum(w)^2)
      double pi = pnorm(zi, 0, 1, FALSE, FALSE); // Pr(Z > zi)

      if(verbose>=5) Rprintf("%u: sw=%2.2f sw2=%2.2f d=%d i=%d si=%2.2f si2=%2.2f mi=%2.2f vi=%2.2f ni=%2.2f zi=%2.2f pi=%2.2f\n", step, sw, sw2, kh_size(discord), i, si, si2, mi, vi, (sw*sw)/sw2, zi, pi);
  
      if(pi > MH_pval){
    extrasteps = step*MH_interval_add+round(runif(0,1));
    finished++;
      }
    }

    if(finished) finished++;
  }

  /* Step finished: record changes. */
  
  if(verbose>=4){
    if(step>=max_MH_interval ) Rprintf("Convergence not achieved after %u M-H steps.\n",step);
    else Rprintf("Convergence achieved after %u M-H steps.\n",step);
  }

  return MCMCDyn1Step_advance(s, dur_inf, stats,
                              maxchanges, nextdiffedge, difftime, difftail, diffhead, diffto,
                              verbose);
}


MCMCDynStatus MCMCDyn1Step_advance(ErgmState *s,
                                   StoreTimeAndLasttoggle *dur_inf,
                                   // Space for output.
                                   double *stats,
                                   unsigned int maxchanges, Edge *nextdiffedge,
                                   kvint *difftime, kvint *difftail, kvint *diffhead, kvint *diffto,
                                   // Verbosity.
                                   int verbose){
  StoreDyadMapInt *discord = dur_inf->discord;
  int t = dur_inf->time;
  
  Network *nwp = s->nwp;
  Model *m = s->m;
  MHProposal *MHp = s->MHp;

  if(nextdiffedge) {
    TailHead dyad;
    kh_foreach_key(discord, dyad,{    
        if(*nextdiffedge<maxchanges){
          // and record the toggle.
          if(difftime) kv_push(int, *difftime, t);
          if(difftail) kv_push(int, *difftail, dyad.tail);
          if(diffhead) kv_push(int, *diffhead, dyad.head);
          if(diffto) kv_push(int, *diffto, GetEdge(dyad.tail,dyad.head,nwp));
          (*nextdiffedge)++;
        }else{
          return(MCMCDyn_TOO_MANY_CHANGES);
        }
      });
  }

  /* If the term has an extension, send it a "TOCK" signal and the set
     of dyads that changed. */
  memset(m->workspace, 0, m->n_stats*sizeof(double)); /* Zero all change stats. */
  SEND_X_SIGNAL_INTO(nwp, m, MHp, m->workspace, TOCK, NULL);
  /* Record network statistics for posterity. */
  if(stats) {
    addonto(stats, m->workspace, m->n_stats);
  }

  return MCMCDyn_OK;
}
