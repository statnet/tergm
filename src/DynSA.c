/*  File src/DynSA.c in package tergm, part of the Statnet suite
 *  of packages for network analysis, https://statnet.org .
 *
 *  This software is distributed under the GPL-3 license.  It is free,
 *  open source, and has the attribution requirements (GPL Section 7) at
 *  https://statnet.org/attribution
 *
 *  Copyright 2008-2019 Statnet Commons
 */
#include "DynSA.h"

void MCMCDynSArun_wrapper(// Observed network.
                 int *tails, int *heads, int *time, int *lasttoggle, int *n_edges,
                 int *n_nodes, int *dflag, int *bipartite, 
                 // Formation terms and proposals.
                 int *nterms, char **funnames, char **sonames,
                 char **MHProposaltype, char **MHProposalpackage,
                 double *inputs, int *nstatsmonitor,
                 // Parameter fittig.
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
                 int *status){
  ErgmState *s;
  StoreDyadMapInt *discord;
    
  Vertex *difftime, *difftail, *diffhead;
  difftime = (Vertex *) Calloc(*maxchanges,Vertex);
  difftail = (Vertex *) Calloc(*maxchanges,Vertex);
  diffhead = (Vertex *) Calloc(*maxchanges,Vertex);

  memset(newnetworktail,0,*maxedges*sizeof(int));
  memset(newnetworkhead,0,*maxedges*sizeof(int));

  MCMCDyn_init_common(tails, heads, *time, lasttoggle, *n_edges,
                      *n_nodes, *dflag, *bipartite,
                      *nterms, *funnames, *sonames, inputs,
                      attribs, maxout, maxin, minout,
                      minin, *condAllDegExact, *attriblength,
                      *MHProposaltype, *MHProposalpackage,
                      &s, &discord);
  Network *nwp = s->nwp;
  Model *m = s->m;

  double *inputdev = Calloc(m->n_stats, double);
  memcpy(inputdev + m->n_stats - *nstatsmonitor, init_dev, (*nstatsmonitor)*sizeof(double));

  *status = MCMCDynSArun(s, discord,
  
             *nstatsmonitor,
    
             eta0,
             inputdev, 
             *runlength,
             WinvGradient, jitter, dejitter, dev_guard, par_guard,
             
             *maxedges, *maxchanges,
             difftime, difftail, diffhead,
             opt_history,
             
             *SA_burnin, *SA_interval, *min_MH_interval, *max_MH_interval, *MH_pval, *MH_interval_add,
             *fVerbose);
  
  memcpy(init_dev, inputdev + m->n_stats - *nstatsmonitor, (*nstatsmonitor)*sizeof(double));
  
  /* record the final network to pass back to R */

  if(*status==MCMCDyn_OK){
    newnetworktail[0]=newnetworkhead[0]=EdgeTree2EdgeList((Vertex*)newnetworktail+1,(Vertex*)newnetworkhead+1,nwp,*maxedges);
    *time = nwp->duration_info->time;

    if(nwp->duration_info){
      lasttoggle[0] = kh_size(nwp->duration_info->lasttoggle);
      TailHead dyad;
      int ts;
      unsigned int i=1;
      kh_foreach(nwp->duration_info->lasttoggle, dyad, ts, {
        lasttoggle[i] = dyad.tail;
        lasttoggle[i+lasttoggle[0]] = dyad.head;
        lasttoggle[i+lasttoggle[0]+lasttoggle[0]] = ts;
        i++;
      });
    }
  }

  Free(difftime);
  Free(difftail);
  Free(diffhead);
  Free(inputdev);
  MCMCDyn_finish_common(s, discord);
}


/*********************
 void MCMCSampleDynPhase12
*********************/
MCMCDynStatus MCMCDynSArun(ErgmState *s, StoreDyadMapInt *discord,
                  int nstatsmonitor,
                  // Model fitting.
                  double *eta, 
                  double *inputdev, // DEViation of the current network's targeted statistics from the target statistics.
                  int runlength,
                  double *WinvGradient, double *jitter, double *dejitter, double *dev_guard, double *par_guard,
                  
                  // Space for output.
                  Edge maxedges, Edge maxchanges,
                  Vertex *difftime, Vertex *difftail, Vertex *diffhead,
                  double *opt_history,
                  // MCMC settings.
                  unsigned int SA_burnin, unsigned int SA_interval, unsigned int min_MH_interval, unsigned int max_MH_interval, double MH_pval, double MH_interval_add,
                  // Verbosity.
                  int fVerbose){
  Network *nwp = s->nwp;
  Model *m = s->m;

  Edge nextdiffedge=1;
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
      MCMCDynStatus status = MCMCDyn1Step(s, discord,
                      eta,
                      inputdev,
                      maxchanges, &nextdiffedge,
                      difftime, difftail, diffhead, NULL,
                      min_MH_interval, max_MH_interval, MH_pval, MH_interval_add,
                      fVerbose);

      if(status==MCMCDyn_TOO_MANY_CHANGES)
        return MCMCDyn_TOO_MANY_CHANGES;
      
      if(nwp->nedges >= maxedges-1)
        return MCMCDyn_TOO_MANY_EDGES;
    }

    // Sampling run
    for(unsigned int j=0;j < SA_interval;j++){
      // Periodically expire older timestamps.
      if(nwp->duration_info->time%TIMESTAMP_HORIZON_FREQ==0) ExpireTimestamps(TIMESTAMP_HORIZON_EDGE, TIMESTAMP_HORIZON_NONEDGE, nwp);

      MCMCDynStatus status = MCMCDyn1Step(s, discord,
                      eta,
                      inputdev,
                      maxchanges, &nextdiffedge,
                      difftime, difftail, diffhead, NULL,
                      min_MH_interval, max_MH_interval, MH_pval, MH_interval_add,
                      fVerbose);

      if(status==MCMCDyn_TOO_MANY_CHANGES)
        return MCMCDyn_TOO_MANY_CHANGES;
      
      if(nwp->nedges >= maxedges-1)
        return MCMCDyn_TOO_MANY_EDGES;
    
      for(unsigned int k=0;k<nstatsmonitor; k++){
        meandev[k]+=dev[k]*1;
        n+=1;
      }
      if (fVerbose>2){
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
    
    if(fVerbose>1){
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
