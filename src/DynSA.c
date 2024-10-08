/*  File src/DynSA.c in package tergm, part of the
 *  Statnet suite of packages for network analysis, https://statnet.org .
 *
 *  This software is distributed under the GPL-3 license.  It is free,
 *  open source, and has the attribution requirements (GPL Section 7) at
 *  https://statnet.org/attribution .
 *
 *  Copyright 2008-2024 Statnet Commons
 */
#include "DynSA.h"
#include "changestats_lasttoggle.h"

SEXP MCMCDynSArun_wrapper(SEXP stateR,
                 SEXP nstatsmonitor,
                 SEXP eta0,
                 SEXP init_dev,
                 SEXP runlength,
                 SEXP WinvGradient,
                 SEXP jitter,
                 SEXP dejitter,
                 SEXP dev_guard,
                 SEXP par_guard,
                 // MCMC settings.
                 SEXP SA_burnin, 
                 SEXP SA_interval,
                 SEXP min_MH_interval,
                 SEXP max_MH_interval,
                 SEXP MH_pval, 
                 SEXP MH_interval_add,
                 SEXP maxedges,
                 SEXP maxchanges,
                 SEXP verbose){    
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

  double *inputdev = R_calloc(m->n_stats, double);
  memcpy(inputdev + m->n_stats - asInteger(nstatsmonitor), REAL(init_dev), asInteger(nstatsmonitor)*sizeof(double));
  
  SEXP opt_history = PROTECT(allocVector(REALSXP, (2*m->n_stats - asInteger(nstatsmonitor))*asInteger(runlength)*asInteger(SA_interval)));
  memset(REAL(opt_history), 0, (2*m->n_stats - asInteger(nstatsmonitor))*asInteger(runlength)*asInteger(SA_interval));
  
  SEXP eta = PROTECT(allocVector(REALSXP, m->n_stats));
  memcpy(REAL(eta), REAL(eta0), m->n_stats*sizeof(double));
  
  SEXP status;
  if(MHp) status = PROTECT(ScalarInteger(MCMCDynSArun(s,
             dur_inf,
  
             asInteger(nstatsmonitor),
    
             REAL(eta),
             inputdev, 
             asInteger(runlength),
             REAL(WinvGradient), REAL(jitter), REAL(dejitter), REAL(dev_guard), REAL(par_guard),
             
             asInteger(maxedges), asInteger(maxchanges),
             REAL(opt_history),
             
             asInteger(SA_burnin), asInteger(SA_interval), asInteger(min_MH_interval), asInteger(max_MH_interval), asReal(MH_pval), asReal(MH_interval_add),
             asInteger(verbose))));
  else status = PROTECT(ScalarInteger(MCMCDyn_MH_FAILED));
  
  SEXP nw_diff = PROTECT(allocVector(REALSXP, asInteger(nstatsmonitor)));
  memcpy(REAL(nw_diff), inputdev + m->n_stats - asInteger(nstatsmonitor), asInteger(nstatsmonitor)*sizeof(double));
  
  const char *outn[] = {"status", "opt.history", "state", "eta", "nw.diff", ""};
  SEXP outl = PROTECT(mkNamed(VECSXP, outn));
  SET_VECTOR_ELT(outl, 0, status);
  SET_VECTOR_ELT(outl, 1, opt_history);
  
  /* record new generated network to pass back to R */
  if(asInteger(status) == MCMCDyn_OK){
    SET_VECTOR_ELT(outl, 2, ErgmStateRSave(s));
  }
  
  SET_VECTOR_ELT(outl, 3, eta);
  SET_VECTOR_ELT(outl, 4, nw_diff);

  ErgmStateDestroy(s);  
  PutRNGstate();  /* Disable RNG before returning */
  UNPROTECT(5);
  return outl;
}


/*********************
 void MCMCSampleDynPhase12
*********************/
MCMCDynStatus MCMCDynSArun(ErgmState *s,
                  StoreTimeAndLasttoggle *dur_inf,
                  int nstatsmonitor,
                  // Model fitting.
                  double *eta, 
                  double *inputdev, // DEViation of the current network's targeted statistics from the target statistics.
                  int runlength,
                  double *WinvGradient, double *jitter, double *dejitter, double *dev_guard, double *par_guard,
                  
                  // Space for output.
                  int maxedges, int maxchanges,
                  double *opt_history,
                  // MCMC settings.
                  unsigned int SA_burnin, unsigned int SA_interval, unsigned int min_MH_interval, unsigned int max_MH_interval, double MH_pval, double MH_interval_add,
                  // Verbosity.
                  int verbose){
  Network *nwp = s->nwp;
  Model *m = s->m;

  unsigned int hist_pos=0, p=m->n_stats - nstatsmonitor, n, rowsize = p*2 + nstatsmonitor;
  double *meandev=(double*)R_alloc(nstatsmonitor,sizeof(double)), *last_jitter=(double*)R_alloc(p,sizeof(double)), *init_eta=(double*)R_alloc(p,sizeof(double));
  memcpy(init_eta, eta, (m->n_stats - nstatsmonitor)*sizeof(double));
  
  double *dev = inputdev + p;
  
  for(unsigned int i=0; i < runlength; i++){
    n = 0;
    for(unsigned int j=0; j<nstatsmonitor; j++){
      meandev[j]=0;
    }

    // Jitter parameters
    for(unsigned int j=0; j<p; j++){
      if(jitter[j]!=0){
        last_jitter[j] = rnorm(0,jitter[j]);
        eta[j] += last_jitter[j];
      }else last_jitter[j] = 0;
    }

    // Burn in
    for(unsigned int j=0;j < SA_burnin;j++){
      R_CheckUserInterrupt();
      MCMCDynStatus status = MCMCDyn1Step(s,
                      dur_inf,
                      eta,
                      inputdev,
                      maxchanges, NULL,
                      NULL, NULL, NULL, NULL,
                      min_MH_interval, max_MH_interval, MH_pval, MH_interval_add,
                      verbose);

      if(status==MCMCDyn_TOO_MANY_CHANGES)
        return MCMCDyn_TOO_MANY_CHANGES;
      
      if(EDGECOUNT(nwp) >= maxedges-1)
        return MCMCDyn_TOO_MANY_EDGES;
    }

    // Sampling run
    for(unsigned int j=0;j < SA_interval;j++){
      R_CheckUserInterrupt();
      MCMCDynStatus status = MCMCDyn1Step(s,
                      dur_inf,
                      eta,
                      inputdev,
                      maxchanges, NULL,
                      NULL, NULL, NULL, NULL,
                      min_MH_interval, max_MH_interval, MH_pval, MH_interval_add,
                      verbose);

      if(status==MCMCDyn_TOO_MANY_CHANGES)
        return MCMCDyn_TOO_MANY_CHANGES;
      
      if(EDGECOUNT(nwp) >= maxedges-1)
        return MCMCDyn_TOO_MANY_EDGES;
    
      for(unsigned int k=0;k<nstatsmonitor; k++){
        meandev[k]+=dev[k]*1;
        n+=1;
      }
      if (verbose>2){
        for(unsigned int k=0; k<p; k++){
          Rprintf("eta[%d] = %f\n", k, eta[k]);
        }
        for(unsigned int k=0; k<nstatsmonitor; k++){
          Rprintf("M_dev[%d] = %f\n", k, dev[k]);
        }

        Rprintf("\n");
      }

      // Record configurations and estimating equation values.
      for(unsigned int j=0; j<p; j++){
        opt_history[hist_pos*rowsize+j] = eta[j];
      }
      for(unsigned int j=0; j<p; j++){
        opt_history[hist_pos*rowsize+p+j] = last_jitter[j];
      }
      for(unsigned int j=0; j<nstatsmonitor; j++){
        opt_history[hist_pos*rowsize+p+p+j] = dev[j];
      }
      hist_pos++;
    }
    
    if(verbose>1){
      for(unsigned int k=0; k<p; k++){
        Rprintf("eta[%d] = %f\n", k, eta[k]);
      }
      for(unsigned int k=0; k<nstatsmonitor; k++){
        Rprintf("meandev[%d] = %f\n", k, meandev[k]/n);
      }
      
      Rprintf("\n");
    }

    // Evaluate mean deviations.
    for(unsigned int j=0; j<nstatsmonitor; j++){
      meandev[j]/=n;
    }
    
    // If the statistics are getting worse by too much, stop updating eta and collect data for the rest of the run.
    for(unsigned int j=0; j<nstatsmonitor; j++){
      if(fabs(meandev[j]) > dev_guard[j]){
        memset(WinvGradient, 0, nstatsmonitor*p*sizeof(double));
        memset(dejitter, 0, p*p*sizeof(double));
      }
    }

    
    // Update formation and dissolution parameters, and cancel the effect of jitter.
    // eta[t+1] = eta[t] - a*(G^-1)*W*(d[t] - G*jit[t])
    //          = eta[t] - a*(G^-1)*W*d[t] + a*(G^-1)*W*G*jit[t]
    for(unsigned int k=0; k<nstatsmonitor; k++){
      for(unsigned int j=0; j<p; j++){
        eta[j] -= WinvGradient[k*p+j] * meandev[k];
      }
    }
    for(unsigned int k=0; k<p; k++){
      for(unsigned int j=0; j<p; j++){
        eta[j] += dejitter[k*p+j] * last_jitter[k];
      }
    }

    // Undo jitter
    for(unsigned int j=0; j<p; j++){
      if(jitter[j]!=0){
        eta[j] -= last_jitter[j];
      }
    }

    // If a parameter has moved suspiciously far, keep it from moving any farther.
    for(unsigned int j=0; j<p; j++){
      double change = eta[j] - init_eta[j];
      if(fabs(change) > par_guard[j]){
        eta[j] = init_eta[j] + (change>0?+1:-1)*par_guard[j]; 
      }
    }
  }

  return MCMCDyn_OK;
}
