/*  File src/DynSA.c in package tergm, part of the Statnet suite
 *  of packages for network analysis, https://statnet.org .
 *
 *  This software is distributed under the GPL-3 license.  It is free,
 *  open source, and has the attribution requirements (GPL Section 7) at
 *  https://statnet.org/attribution
 *
 *  Copyright 2008-2020 Statnet Commons
 */
#include "DynSA.h"

void MCMCDynSArun_wrapper(// Observed network.
			     int *tails, int *heads, int *time, int *lasttoggle_flag, int *lasttoggle, int *n_edges,
			     int *n_nodes, int *dflag, int *bipartite, 
			     // Formation terms and proposals.
			     int *F_nterms, char **F_funnames, char **F_sonames,
			     char **F_MHProposaltype, char **F_MHProposalpackage,
			     double *F_inputs, 
			     // Dissolution terms and proposals.
			     int *D_nterms, char **D_funnames, char **D_sonames, 
			     char **D_MHProposaltype, char **D_MHProposalpackage,
			     double *D_inputs,
			     // Parameter fittig.
			     double *eta0, 
			     int *M_nterms, char **M_funnames, char **M_sonames, double *M_inputs,
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

  Network *nwp;
  Model *F_m, *D_m, *M_m;
  MHProposal *F_MH, *D_MH;
  
  if(!*lasttoggle_flag) lasttoggle = NULL;

  Vertex *difftime, *difftail, *diffhead;
  difftime = (Vertex *) Calloc(*maxchanges,Vertex);
  difftail = (Vertex *) Calloc(*maxchanges,Vertex);
  diffhead = (Vertex *) Calloc(*maxchanges,Vertex);

  memset(newnetworktail,0,*maxedges*sizeof(int));
  memset(newnetworkhead,0,*maxedges*sizeof(int));

  MCMCDyn_init_common(tails, heads, *time, lasttoggle, *n_edges,
		      *n_nodes, *dflag, *bipartite, &nwp,
		      *F_nterms, *F_funnames, *F_sonames, F_inputs, &F_m,
		      *D_nterms, *D_funnames, *D_sonames, D_inputs, &D_m,
		      *M_nterms, *M_funnames, *M_sonames, M_inputs, &M_m,
		      attribs, maxout, maxin, minout,
		      minin, *condAllDegExact, *attriblength,
		      *F_MHProposaltype, *F_MHProposalpackage, &F_MH,
		      *D_MHProposaltype, *D_MHProposalpackage, &D_MH,
		      *fVerbose);

  *status = MCMCDynSArun(nwp,
			    
			 F_m, F_MH,
			 D_m, D_MH,
			 
			 eta0, M_m,
			 init_dev, 
			 *runlength,
			 WinvGradient, jitter, dejitter, dev_guard, par_guard,
			 
			 *maxedges, *maxchanges,
			 difftime, difftail, diffhead,
			 opt_history,
			 
			 *SA_burnin, *SA_interval, *min_MH_interval, *max_MH_interval, *MH_pval, *MH_interval_add,
			 *fVerbose);
  
  /* record the final network to pass back to R */

  if(*status==MCMCDyn_OK){
    newnetworktail[0]=newnetworkhead[0]=EdgeTree2EdgeList((Vertex*)newnetworktail+1,(Vertex*)newnetworkhead+1,nwp,*maxedges);
    *time = nwp->duration_info.time;
    if(nwp->duration_info.lasttoggle)
    memcpy(lasttoggle, nwp->duration_info.lasttoggle, sizeof(int)*DYADCOUNT(*n_nodes, *bipartite, *dflag));
  }

  MCMCDyn_finish_common(nwp, F_m, D_m, M_m, F_MH, D_MH);
  Free(difftime);
  Free(difftail);
  Free(diffhead);
}


/*********************
 void MCMCSampleDynPhase12
*********************/
MCMCDynStatus MCMCDynSArun(// Observed and discordant network.
			      Network *nwp,
			      // Formation terms and proposals.
			      Model *F_m, MHProposal *F_MH,
			      // Dissolution terms and proposals.
			      Model *D_m, MHProposal *D_MH,
			      // Model fitting.
			      double *eta, 
			      Model *M_m,
			      double *dev, // DEViation of the current network's targeted statistics from the target statistics.
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
  Edge nextdiffedge=1;
  unsigned int hist_pos=0, p=F_m->n_stats+D_m->n_stats, n, rowsize = p*2 + M_m->n_stats;
  double *meandev=(double*)R_alloc(M_m->n_stats,sizeof(double)), *last_jitter=(double*)R_alloc(p,sizeof(double)), *init_eta=(double*)R_alloc(p,sizeof(double));
  memcpy(init_eta, eta, (F_m->n_stats+D_m->n_stats)*sizeof(double));

  for(unsigned int i=0; i < runlength; i++){
    n = 0;
    for(unsigned int j=0; j<M_m->n_stats; j++){
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
      MCMCDynStatus status = MCMCDyn1Step(nwp,
					  F_m, F_MH, eta,
					  D_m, D_MH, eta+F_m->n_stats,
					  M_m,
					  0,
					  NULL, NULL, dev,
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
      MCMCDynStatus status = MCMCDyn1Step(nwp,
					  F_m, F_MH, eta,
					  D_m, D_MH, eta+F_m->n_stats,
					  M_m,
					  0,
					  NULL, NULL, dev,
					  maxchanges, &nextdiffedge,
					  difftime, difftail, diffhead, NULL,
					  min_MH_interval, max_MH_interval, MH_pval, MH_interval_add,
					  fVerbose);

      if(status==MCMCDyn_TOO_MANY_CHANGES)
        return MCMCDyn_TOO_MANY_CHANGES;
      
      if(nwp->nedges >= maxedges-1)
	return MCMCDyn_TOO_MANY_EDGES;
    
      for(unsigned int k=0;k<M_m->n_stats; k++){
	meandev[k]+=dev[k]*1;
	n+=1;
      }
      if (fVerbose>2){
	for(unsigned int k=0; k<p; k++){
	  Rprintf("eta[%d] = %f\n", k, eta[k]);
	}
	for(unsigned int k=0; k<M_m->n_stats; k++){
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
      for(unsigned int j=0; j<M_m->n_stats; j++){
	opt_history[hist_pos*rowsize+p+p+j] = dev[j];
      }
      hist_pos++;
    }
    
    if(fVerbose>1){
      for(unsigned int k=0; k<p; k++){
	Rprintf("eta[%d] = %f\n", k, eta[k]);
      }
      for(unsigned int k=0; k<M_m->n_stats; k++){
	Rprintf("meandev[%d] = %f\n", k, meandev[k]/n);
      }
      
      Rprintf("\n");
    }

    // Evaluate mean deviations.
    for(unsigned int j=0; j<M_m->n_stats; j++){
      meandev[j]/=n;
    }
    
    // If the statistics are getting worse by too much, stop updating eta and collect data for the rest of the run.
    for(unsigned int j=0; j<M_m->n_stats; j++){
      if(fabs(meandev[j]) > dev_guard[j]){
	memset(WinvGradient, 0, M_m->n_stats*p*sizeof(double));
	memset(dejitter, 0, p*p*sizeof(double));
      }
    }

    
    // Update formation and dissolution parameters, and cancel the effect of jitter.
    // eta[t+1] = eta[t] - a*(G^-1)*W*(d[t] - G*jit[t])
    //          = eta[t] - a*(G^-1)*W*d[t] + a*(G^-1)*W*G*jit[t]
    for(unsigned int k=0; k<M_m->n_stats; k++){
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

