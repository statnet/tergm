/*  File src/MCMCDyn.c in package tergm, part of the Statnet suite
 *  of packages for network analysis, http://statnet.org .
 *
 *  This software is distributed under the GPL-3 license.  It is free,
 *  open source, and has the attribution requirements (GPL Section 7) at
 *  http://statnet.org/attribution
 *
 *  Copyright 2008-2018 Statnet Commons
 */
#include "MCMCDyn.h"
#include <stdbool.h>

/*****************
 Note on undirected networks:  For j<k, edge {j,k} should be stored
 as (j,k) rather than (k,j).  In other words, only directed networks
 should have (k,j) with k>j.
*****************/

void MCMCDyn_init_common(int *tails, int *heads, int time, int *lasttoggle, int n_edges,
			 int n_nodes, int dflag, int bipartite, Network **nwp,
			 
			 int F_nterms, char *F_funnames, char *F_sonames, double *F_inputs, Model **F_m,
			 int D_nterms, char *D_funnames, char *D_sonames, double *D_inputs, Model **D_m,
			 int M_nterms, char *M_funnames, char *M_sonames, double *M_inputs, Model **M_m,
			 
			 int *attribs, int *maxout, int *maxin, int *minout,
			 int *minin, int condAllDegExact, int attriblength,
			 
			 char *F_MHProposaltype, char *F_MHProposalpackage, MHProposal **F_MH,
			 char *D_MHProposaltype, char *D_MHProposalpackage, MHProposal **D_MH,
			 int fVerbose){
  
  if(F_nterms) *F_m=ModelInitialize(F_funnames, F_sonames, &F_inputs, F_nterms); else if(F_m) *F_m=NULL;
  if(D_nterms) *D_m=ModelInitialize(D_funnames, D_sonames, &D_inputs, D_nterms); else if(D_m) *D_m=NULL;
  if(M_nterms) *M_m=ModelInitialize(M_funnames, M_sonames, &M_inputs, M_nterms); else if(M_m) *M_m=NULL;

  nwp[0]=NetworkInitialize(tails, heads, n_edges, 
                          n_nodes, dflag, bipartite, 1, time, lasttoggle);

  // This code initializes the discordance network. Note that both
  // F_MH and D_MH should expect MH->discord[0] to be the discordance
  // network relative to the last time step. However, they may want to
  // maintain other discordance networks, in which case they might
  // choose to allocate their own MH->discord vectors. If they do,
  // discordl won't be allocated.
  if(F_MH || D_MH){
    Network *discord=NetworkInitialize(NULL, NULL, 0,
				       n_nodes, dflag, bipartite, 0, 0, NULL);
    Network **discordl=NULL;
  
    if(F_MH){
      *F_MH = MHProposalInitialize(F_MHProposaltype, F_MHProposalpackage, F_inputs, fVerbose, *nwp, attribs, maxout, maxin, minout, minin,
				   condAllDegExact, attriblength);
      if((*F_MH)->discord==NULL){
	discordl = (Network**) Calloc(2,Network*); // A space for the sentinel NULL pointer.
	(*F_MH)->discord = discordl;
      }
      (*F_MH)->discord[0] = discord;
    }
    
    if(D_MH){
      *D_MH = MHProposalInitialize(D_MHProposaltype, D_MHProposalpackage, D_inputs, fVerbose, *nwp, attribs, maxout, maxin, minout, minin,
				   condAllDegExact, attriblength);
      if((*D_MH)->discord==NULL){
	if(discordl==NULL)
	  discordl = (Network**) Calloc(2,Network*); // A space for the sentinel NULL pointer.
	(*D_MH)->discord = discordl;
      }
      (*D_MH)->discord[0] = discord;
    }
  }

  GetRNGstate();  /* R function enabling uniform RNG. It needs to come after NetworkInitialize and MH_init, since they may call GetRNGstate as well. */  

}

void MCMCDyn_finish_common(Network *nwp,
			   Model *F_m,
			   Model *D_m,
			   Model *M_m,
			   MHProposal *F_MH,
			   MHProposal *D_MH){
  if(F_MH || D_MH){
    bool same_discord = F_MH->discord==D_MH->discord;
    for(Network **nwp=F_MH->discord; *nwp!=NULL; nwp++){
      NetworkDestroy(*nwp);
    }
    // First discordance network has already been destroyed.
    for(Network **nwp=D_MH->discord+1; *nwp!=NULL; nwp++){
      NetworkDestroy(*nwp);
    }
    
    Free(F_MH->discord);
    F_MH->discord=NULL;
    if(!same_discord) Free(D_MH->discord);
    D_MH->discord=NULL; // Always set to NULL, since it may have been freed with F_MH->discord.

    if(F_MH) MHProposalDestroy(F_MH);
    if(D_MH) MHProposalDestroy(D_MH);
  }
  if(F_m) ModelDestroy(F_m);
  if(D_m) ModelDestroy(D_m);
  if(M_m) ModelDestroy(M_m);
  NetworkDestroy(nwp);
  PutRNGstate();  /* Disable RNG before returning */

}

/*****************
 void MCMCDyn_wrapper

 Wrapper for a call from R.
*****************/
void MCMCDyn_wrapper(// Starting network.
		     int *tails, int *heads, int *time, int *lasttoggle, int *n_edges,
		     int *n_nodes, int *dflag, int *bipartite,
		     // Formation terms and proposals.
		     int *F_nterms, char **F_funnames, char **F_sonames, 
		     char **F_MHProposaltype, char **F_MHProposalpackage,
		     double *F_inputs, double *F_eta, 
		     // Dissolution terms and proposals.
		     int *D_nterms, char **D_funnames, char **D_sonames,
		     char **D_MHProposaltype, char **D_MHProposalpackage,
		     double *D_inputs, double *D_eta,
		     // Monitored terms.
		     int *M_nterms, char **M_funnames, char **M_sonames,  double *M_inputs,
		     // Degree bounds.
		     int *attribs, int *maxout, int *maxin, int *minout,
		     int *minin, int *condAllDegExact, int *attriblength, 
		     // MCMC settings.
		     int *nsteps,  int *min_MH_interval, int *max_MH_interval, double *MH_pval, double *MH_interval_add,
		     int *burnin, int *interval,  
		     // Space for output.
		     int *F_collect, double *F_sample, 
		     int *D_collect, double *D_sample, 
		     double *M_sample,
		     int *maxedges,
		     int *newnetworktails, int *newnetworkheads, 
		     int *maxchanges,
		     int *log_changes,
		     int *diffnetworktime, int *diffnetworktail, int *diffnetworkhead, int *diffnetworkto,
		     // Verbosity.
		     int *fVerbose,
		     int *status){
  Network *nwp;
  Model *F_m, *D_m, *M_m;
  MHProposal *F_MH, *D_MH;

  if(*lasttoggle == 0) lasttoggle = NULL;

  Vertex *difftime, *difftail, *diffhead;
  int *diffto;

  if(*log_changes){
    difftime = (Vertex *) diffnetworktime;
    difftail = (Vertex *) diffnetworktail;
    diffhead = (Vertex *) diffnetworkhead;
    diffto = (int *) diffnetworkto;
  }else{
    difftime = (Vertex *) Calloc(*maxchanges,Vertex);
    difftail = (Vertex *) Calloc(*maxchanges,Vertex);
    diffhead = (Vertex *) Calloc(*maxchanges,Vertex);
    diffto = NULL;
  }
  
  if(*log_changes){
    memset(difftime,0,*maxchanges*sizeof(Vertex));
    memset(difftail,0,*maxchanges*sizeof(Vertex));
    memset(diffhead,0,*maxchanges*sizeof(Vertex));
    memset(diffto,0,*maxchanges*sizeof(int));
  }

  MCMCDyn_init_common(tails, heads, *time, lasttoggle, *n_edges,
		      *n_nodes, *dflag, *bipartite, &nwp,
		      *F_nterms, *F_funnames, *F_sonames, F_inputs, &F_m,
		      *D_nterms, *D_funnames, *D_sonames, D_inputs, &D_m,
		      *M_nterms, *M_nterms ? *M_funnames : NULL, *M_nterms ? *M_sonames : NULL, M_inputs, &M_m,   // I am not sure  whether the ?:s are necessary here, but if it keeps us from dereferencing a NULL pointer...
		      attribs, maxout, maxin, minout,
		      minin, *condAllDegExact, *attriblength, 
		      *F_MHProposaltype, *F_MHProposalpackage, &F_MH,
		      *D_MHProposaltype, *D_MHProposalpackage, &D_MH,
		      *fVerbose);

  *status = MCMCSampleDyn(nwp,
			  F_m, F_MH, F_eta,
			  D_m, D_MH, D_eta,
			  M_m,
			  *F_collect?F_sample:NULL, *D_collect?D_sample:NULL, M_m?M_sample:NULL, *maxedges, *maxchanges, *log_changes, difftime, difftail, diffhead, diffto,
			  *nsteps, *min_MH_interval, *max_MH_interval, *MH_pval, *MH_interval_add, *burnin, *interval,
			  *fVerbose);
   
  /* record new generated network to pass back to R */

  if(*status == MCMCDyn_OK && *maxedges>0){
    newnetworktails[0]=newnetworkheads[0]=EdgeTree2EdgeList(newnetworktails+1,newnetworkheads+1,nwp,*maxedges-1);
    *time = nwp->duration_info.time;

    if(nwp->duration_info.lasttoggle)
    memcpy(lasttoggle, nwp->duration_info.lasttoggle, sizeof(int)*DYADCOUNT(*n_nodes, *bipartite, *dflag));
  }

  if(!*log_changes){
    Free(difftime);
    Free(difftail);
    Free(diffhead);
  }

  MCMCDyn_finish_common(nwp, F_m, D_m, M_m, F_MH, D_MH);
}

/*********************
 MCMCDynStatus MCMCSampleDyn

 Using the parameters contained in the array F_eta, obtain the
 network statistics for a sample of size nsteps.  burnin is the
 initial number of Markov chain steps before sampling anything
 and interval is the number of MC steps between successive 
 networks in the sample.  Put all the sampled statistics into
 the F_statistics array. 
*********************/
MCMCDynStatus MCMCSampleDyn(// Observed and discordant network.
			    Network *nwp,
			    // Formation terms and proposals.
			    Model *F_m, MHProposal *F_MH,
			    double *F_eta,
			    // Dissolution terms and proposals.
			    Model *D_m, MHProposal *D_MH,
			    double *D_eta,
			    // Monitored terms.
			    Model *M_m,
			    // Space for output.
			    double *F_stats, double *D_stats, double *M_stats,
			    Edge maxedges,
			    Edge maxchanges,
			    int log_changes,
			    Vertex *difftime, Vertex *difftail, Vertex *diffhead, int *diffto,		    
			    // MCMC settings.
			    unsigned int nsteps, unsigned int min_MH_interval, unsigned int max_MH_interval, double MH_pval, double MH_interval_add,
			    unsigned int burnin, unsigned int interval, 
			    // Verbosity.
			    int fVerbose){

  int i, j;
  Edge nextdiffedge=1;


  /*if (fVerbose)
    Rprintf("Total m->n_stats is %i; total nsteps is %d\n",
    F_m->n_stats,nsteps);*/
  
  
  /* Burn in step. */

  for(i=0;i<burnin;i++){
    MCMCDynStatus status = MCMCDyn1Step(nwp,
					F_m, F_MH, F_eta, D_m, D_MH, D_eta, M_m,
					log_changes, F_stats, D_stats, M_stats,
					maxchanges, &nextdiffedge, difftime, difftail, diffhead, diffto,
					min_MH_interval, max_MH_interval, MH_pval, MH_interval_add, fVerbose);
    // Check that we didn't run out of log space.
    if(status==MCMCDyn_TOO_MANY_CHANGES)
      return MCMCDyn_TOO_MANY_CHANGES;
  
    // If we need to return a network, then stop right there, since the network is too big to return, so stop early.
    if(maxedges!=0 && nwp->nedges >= maxedges-1)
      return MCMCDyn_TOO_MANY_EDGES;
  }
  
  //Rprintf("MCMCSampleDyn post burnin numdissolve %d\n", *numdissolve);
  
  if (fVerbose){
    Rprintf("Returned from STERGM burnin\n");
  }
  
  /* Now sample networks */
  for (i=0; i < nsteps; i++){
    /* Set current vector of stats equal to previous vector */
    if(F_stats){
      for (j=0; j<F_m->n_stats; j++){
	F_stats[j+F_m->n_stats] = F_stats[j];
      }
      F_stats += F_m->n_stats;
    }

    if(D_stats){
      for (j=0; j<D_m->n_stats; j++){
	D_stats[j+D_m->n_stats] = D_stats[j];
      }
      D_stats += D_m->n_stats;
    }
    
    if(M_stats){
      for (j=0; j<M_m->n_stats; j++){
	M_stats[j+M_m->n_stats] = M_stats[j];
      }
      M_stats += M_m->n_stats;
    }


    /* This then adds the change statistics to these values */
    for(j=0;j<interval;j++){
      MCMCDynStatus status = MCMCDyn1Step(nwp,
					  F_m, F_MH, F_eta, D_m, D_MH, D_eta, M_m,
					  log_changes, F_stats, D_stats, M_stats,
					  maxchanges, &nextdiffedge, difftime, difftail, diffhead, diffto,
					  min_MH_interval, max_MH_interval, MH_pval, MH_interval_add, fVerbose);
      
      // Check that we didn't run out of log space.
      if(status==MCMCDyn_TOO_MANY_CHANGES)
        return MCMCDyn_TOO_MANY_CHANGES;
      
      // If we need to return a network, then stop right there, since the network is too big to return, so stop early.
      if(maxedges!=0 && nwp->nedges >= maxedges-1)
	return MCMCDyn_TOO_MANY_EDGES;

    }
    
    //Rprintf("MCMCSampleDyn loop numdissolve %d\n", *numdissolve);
    if (fVerbose){
      if( ((3*i) % nsteps)<3 && nsteps > 500){
        Rprintf("Advanced %d time steps.\n", i);
      }
    }
  }

  if(log_changes) difftime[0]=difftail[0]=diffhead[0]=diffto[0]=nextdiffedge-1;
  return MCMCDyn_OK;
}

/* Helper function to run the formation or the dissolution side of the process, 
   depending on proposal, statistics, and parameters passed. 

   NOTE: Here, "stats" and "m" are for the process that "drives" the sampling 
   in this phase, while "O_stats" and "O_m" are for the other process (whose
   statistcs are still affected and need to be updated).
*/
 void MCMCDyn1Step_sample(MHProposal *MH,
			  double *par,
			  unsigned int min_MH_interval, unsigned int max_MH_interval, double MH_pval, double MH_interval_add,
			  Network *nwp,
			  Model *m, int fVerbose){
  double cutoff;
  double 
    si = 0, // sum of increments
    si2 = 0, // sum of squared increments
    s = 0, // sum of weights 
    s2 = 0 // sum of squared weights
    ;
  double sdecay = 1 - 1.0/min_MH_interval;

  MH->ntoggles = 0;
  (*(MH->func))(MH, nwp); /* Call MH proposal function to initialize */
  
  unsigned int step=0; // So that we could print out the number of steps later.
  for(unsigned int finished = 0, extrasteps = 0; step < max_MH_interval && finished < extrasteps+1; step++) {
    unsigned int prev_discord = MH->discord[0]->nedges;

    MH->logratio = 0;
    (*(MH->func))(MH, nwp); /* Call MH function to propose toggles */
    //      Rprintf("Back from proposal; step=%d\n",step);

    // Proposal failed.
    if(MH->toggletail[0]==MH_FAILED){
      switch(MH->togglehead[0]){
      case MH_UNRECOVERABLE:
	error("Something very bad happened during proposal. Memory has not been deallocated, so restart R soon.");
	
      case MH_IMPOSSIBLE:
	// Proceed to the next phase.
	return;
	
      case MH_UNSUCCESSFUL:
      case MH_CONSTRAINT:
	continue;
      }
    }

    ChangeStats(MH->ntoggles, MH->toggletail, MH->togglehead, nwp, m);
    
    //  Rprintf("change stats:"); 
    /* Calculate inner product */
    double ip = 0;
    for (unsigned int i=0; i<m->n_stats; i++){
      ip += par[i] * m->workspace[i];
      //  Rprintf("%f ", m->workspace[i]); 
    }
    //  Rprintf("\n ip %f dedges %f\n", ip, m->workspace[0]); 
    /* The logic is to set exp(cutoff) = exp(ip) * qratio ,
       then let the MH probability equal min{exp(cutoff), 1.0}.
       But we'll do it in log space instead.  */
    cutoff = ip + MH->logratio;
    
    /* if we accept the proposed network */
    if (cutoff >= 0.0 || log(unif_rand()) < cutoff) { 
      /* Hold off updating timesteamps until the changes are committed,
      which doesn't happen until later. */
      for (unsigned int i=0; i < MH->ntoggles; i++){
        ToggleEdge(MH->toggletail[i], MH->togglehead[i], nwp);
	
	if(MH->discord)
	  for(Network **nwd=MH->discord; *nwd!=NULL; nwd++){
	    ToggleEdge(MH->toggletail[i],  MH->togglehead[i], *nwd);
	  }
      }
      /* Do NOT record network statistics for posterity yet. */
    }

    int i = MH->discord[0]->nedges - prev_discord;
    s *= sdecay; si *= sdecay;
    s++; si += i;
    s2 *= sdecay*sdecay; si2 *= sdecay;
    s2++; si2 += i*i;
         
    if(step >= min_MH_interval && !finished) { 
      // Now, perform the test:
      double mi = (double)si / s, mi2 = (double)si2 / s;
      
      double vi = mi2 - mi*mi;
      double zi = mi / sqrt(vi * s2/(s*s)); // denom = sqrt(sum(w^2 * v)/sum(w)^2)
      double pi = pnorm(zi, 0, 1, FALSE, FALSE); // Pr(Z > zi)

      if(fVerbose>=5) Rprintf("%u: s=%2.2f s2=%2.2f d=%d i=%d si=%2.2f si2=%2.2f mi=%2.2f vi=%2.2f ni=%2.2f zi=%2.2f pi=%2.2f\n", step, s, s2, MH->discord[0]->nedges, i, si, si2, mi, vi, (s*s)/s2, zi, pi);
  
      if(pi > MH_pval){
	extrasteps = step*MH_interval_add+round(runif(0,1));
	finished++;
      }
    }

    if(finished) finished++;
  }
  if(fVerbose>=4){
    if(step>=max_MH_interval ) Rprintf("Convergence not achieved after %u M-H steps.\n",step);
    else Rprintf("Convergence achieved after %u M-H steps.\n",step);
  }

}

/*
  MCMCDyn1Step_advance
  Applies a list of toggles to a network, updating change statistics and time stamps.
*/
void MCMCDyn1Step_advance(unsigned int ntoggles,
			  Vertex *difftail, Vertex *diffhead,
			  Network *nwp,
			  Model *F_m, double *F_stats,
			  Model *D_m, double *D_stats,
			  Model *M_m, double *M_stats){

  if(F_stats){
    ChangeStats(ntoggles,difftail,diffhead,nwp,F_m);
    for (unsigned int i = 0; i < F_m->n_stats; i++)
      F_stats[i] += F_m->workspace[i];
  }
  
  if(D_stats){
    ChangeStats(ntoggles,difftail,diffhead,nwp,D_m);
    for (unsigned int i = 0; i < D_m->n_stats; i++)
      D_stats[i] += D_m->workspace[i];
  }
  
  if(M_stats){
    ChangeStats(ntoggles,difftail,diffhead,nwp,M_m);
    for (unsigned int i = 0; i < M_m->n_stats; i++)
      M_stats[i] += M_m->workspace[i];
  }

  for(Edge i=0;i<ntoggles;i++){
    ToggleEdge(difftail[i],diffhead[i],nwp);
  }

  for(Edge i=0;i<ntoggles;i++){
    TouchEdge(difftail[i],diffhead[i],nwp);
  }
  nwp->duration_info.time++;

}

/* 
   MCMCDyn1Step_record_reset
   Record a set of changes to the state of the sampler at the end of a phase:
   * record new generated network differences to pass back to R
   * undo the toggles to nwp[0]
   * empty the discordant network (nwp[1])
   * return the number of edges toggled, or -1 if ran out of space
*/
int MCMCDyn1Step_record_reset(Edge maxchanges,
			      Vertex *difftime, Vertex *difftail, Vertex *diffhead,
			      Network **nwp, 
			      Edge *nextdiffedge){
  Vertex tail, head;
  const unsigned int t=nwp[0]->duration_info.time+1; // Note that the toggle only takes effect on the next time step.
  Edge ntoggles = nwp[1]->nedges;
  
  for(unsigned int i=0; i<ntoggles; i++){
    FindithEdge(&tail, &head, 1, nwp[1]); // Grab the next edge that changed;
    ToggleEdge(tail, head, nwp[1]); // delete it from nwp[1];
    ToggleEdge(tail, head, nwp[0]); // undo the toggle in nwp[0];
    
    if(*nextdiffedge<maxchanges){
      // and record the toggle.
      difftime[*nextdiffedge] = t; // Effective toggle time is t+1.
      difftail[*nextdiffedge] = tail;
      diffhead[*nextdiffedge] = head;
      (*nextdiffedge)++;
    }else{
      return(-1);
    }
  }
  return(ntoggles);
}

/*********************
 void MCMCDyn1Step

 Simulate evolution of a dynamic network for nsteps steps.
*********************/
MCMCDynStatus MCMCDyn1Step(// Observed and discordant network.
		  Network *nwp,
		  // Formation terms and proposals.
		  Model *F_m, MHProposal *F_MH, double *F_eta,
		  // Dissolution terms and proposals.
		  Model *D_m, MHProposal *D_MH, double *D_eta,
		  // Monitored statistics.
		  Model *M_m,
		  // Space for output.
		  unsigned log_changes,
		  double *F_stats, double *D_stats, double *M_stats,
		  unsigned int maxchanges, Edge *nextdiffedge,
		  Vertex *difftime, Vertex *difftail, Vertex *diffhead, int *diffto,
		  // MCMC settings.
		  unsigned int min_MH_interval, unsigned int max_MH_interval, double MH_pval, double MH_interval_add,
		  // Verbosity.
		  int fVerbose){
  
  Edge ntoggles;
  int ntoggles_status;

  Edge nde=1;
  if(nextdiffedge) nde=*nextdiffedge;

  /* Run the dissolution process. */
  MCMCDyn1Step_sample(D_MH, D_eta, min_MH_interval, max_MH_interval, MH_pval, MH_interval_add, nwp, D_m, fVerbose);
  ntoggles_status = MCMCDyn1Step_record_reset(maxchanges, difftime, difftail, diffhead, (Network*[]){nwp,D_MH->discord[0]}, &nde);
  if(ntoggles_status<0) return MCMCDyn_TOO_MANY_CHANGES;
  if(diffto) for(unsigned int i = 1; i<=ntoggles_status; i++) diffto[nde-i] = 0;
  ntoggles = ntoggles_status;
  
  /* Run the formation process. */
  MCMCDyn1Step_sample(F_MH, F_eta, min_MH_interval, max_MH_interval, MH_pval, MH_interval_add, nwp, F_m, fVerbose);
  ntoggles_status = MCMCDyn1Step_record_reset(maxchanges, difftime, difftail, diffhead, (Network*[]){nwp,F_MH->discord[0]}, &nde);
  if(ntoggles_status<0) return MCMCDyn_TOO_MANY_CHANGES;
  if(diffto) for(unsigned int i = 1; i<=ntoggles_status; i++) diffto[nde-i] = 1;
  ntoggles += ntoggles_status;
  
  /* Commit both. */
  MCMCDyn1Step_advance(ntoggles, difftail+nde-ntoggles, diffhead+nde-ntoggles, nwp, F_m, F_stats, D_m, D_stats, M_m, M_stats);

  // If we don't keep a log of toggles, reset the position to save space.
  if(log_changes)
    *nextdiffedge=nde;
  return MCMCDyn_OK;
}
