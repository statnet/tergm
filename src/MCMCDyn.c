/*  File src/MCMCDyn.c in package tergm, part of the Statnet suite
 *  of packages for network analysis, https://statnet.org .
 *
 *  This software is distributed under the GPL-3 license.  It is free,
 *  open source, and has the attribution requirements (GPL Section 7) at
 *  https://statnet.org/attribution
 *
 *  Copyright 2008-2019 Statnet Commons
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
			 
			 int nterms, char *funnames, char *sonames, double *inputs, Model **m,
			 
			 int *attribs, int *maxout, int *maxin, int *minout,
			 int *minin, int condAllDegExact, int attriblength,
			 
			 char *MHProposaltype, char *MHProposalpackage, MHProposal **MHp,
                         StoreDyadMapInt **discord,
			 int fVerbose){
  GetRNGstate();  /* R function enabling uniform RNG. It needs to come after NetworkInitialize and MH_init, since they may call GetRNGstate as well. */  
  
  *m=ModelInitialize(funnames, sonames, &inputs, nterms);

  *nwp=NetworkInitialize((Vertex*)tails, (Vertex*)heads, n_edges, 
                          n_nodes, dflag, bipartite, 1, time, lasttoggle);
  
  /* Trigger initial storage update */
  InitStats(*nwp, *m);

  if(MHp){
    *MHp = MHProposalInitialize(MHProposaltype, MHProposalpackage, inputs, fVerbose, *nwp, attribs, maxout, maxin, minout, minin,
                                condAllDegExact, attriblength, (*m)->termarray->aux_storage);
  }

  *discord = kh_init(DyadMapInt); (*discord)->directed = dflag;
}

void MCMCDyn_finish_common(Network *nwp,
			   Model *m,
			   MHProposal *MHp,
                           StoreDyadMapInt *discord){
  if(MHp) MHProposalDestroy(MHp,nwp);
  ModelDestroy(nwp,m);
  NetworkDestroy(nwp);
  kh_destroy(DyadMapInt, discord);
  PutRNGstate();  /* Disable RNG before returning */

}

/*****************
 void MCMCDyn_wrapper

 Wrapper for a call from R.
*****************/
void MCMCDyn_wrapper(// Starting network.
		     int *tails, int *heads, int *time, int *lasttoggle, int *n_edges,
		     int *n_nodes, int *dflag, int *bipartite,
		     // Terms and proposals.
		     int *nterms, char **funnames, char **sonames, 
		     char **MHProposaltype, char **MHProposalpackage,
		     double *inputs, double *eta, 
		     // Degree bounds.
		     int *attribs, int *maxout, int *maxin, int *minout,
		     int *minin, int *condAllDegExact, int *attriblength, 
		     // MCMC settings.
		     int *nsteps,  int *min_MH_interval, int *max_MH_interval, double *MH_pval, double *MH_interval_add,
		     int *burnin, int *interval,  
		     // Space for output.
		     int *collect, double *sample, 
		     int *maxedges,
		     int *newnetworktails, int *newnetworkheads, 
		     int *maxchanges,
		     int *log_changes,
		     int *diffnetworktime, int *diffnetworktail, int *diffnetworkhead, int *diffnetworkto,
		     // Verbosity.
		     int *fVerbose,
		     int *status){
  Network *nwp;
  Model *m;
  MHProposal *MHp;
  StoreDyadMapInt *discord;

  Vertex *difftime=NULL, *difftail=NULL, *diffhead=NULL;
  int *diffto=NULL;

  if(*log_changes){
    difftime = (Vertex *) diffnetworktime;
    difftail = (Vertex *) diffnetworktail;
    diffhead = (Vertex *) diffnetworkhead;
    diffto = (int *) diffnetworkto;
    memset(difftime,0,*maxchanges*sizeof(Vertex));
    memset(difftail,0,*maxchanges*sizeof(Vertex));
    memset(diffhead,0,*maxchanges*sizeof(Vertex));
    memset(diffto,0,*maxchanges*sizeof(int));
  }

  MCMCDyn_init_common(tails, heads, *time, lasttoggle, *n_edges,
		      *n_nodes, *dflag, *bipartite, &nwp,
		      *nterms, *funnames, *sonames, inputs, &m,
		      attribs, maxout, maxin, minout,
		      minin, *condAllDegExact, *attriblength, 
		      *MHProposaltype, *MHProposalpackage, &MHp,
                      &discord,
		      *fVerbose);

  *status = MCMCSampleDyn(nwp, discord,
			  m, MHp, eta,
			  *collect?sample:NULL, *maxedges, *maxchanges, *log_changes, difftime, difftail, diffhead, diffto,
			  *nsteps, *min_MH_interval, *max_MH_interval, *MH_pval, *MH_interval_add, *burnin, *interval,
			  *fVerbose);
   
  /* record new generated network to pass back to R */

  if(*status == MCMCDyn_OK && *maxedges>0){
    newnetworktails[0]=newnetworkheads[0]=EdgeTree2EdgeList((Vertex*)newnetworktails+1,(Vertex*)newnetworkheads+1,nwp,*maxedges-1);
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

  MCMCDyn_finish_common(nwp, m, MHp, discord);
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
MCMCDynStatus MCMCSampleDyn(// Observed and discordant network.
			    Network *nwp, StoreDyadMapInt *discord,
			    // terms and proposals.
			    Model *m, MHProposal *MHp,
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
			    int fVerbose){

  int i, j;
  Edge nextdiffedge=1;


  /*if (fVerbose)
    Rprintf("Total m->n_stats is %i; total nsteps is %d\n",
    m->n_stats,nsteps);*/
  
  
  /* Burn in step. */

  for(i=0;i<burnin;i++){
    MCMCDynStatus status = MCMCDyn1Step(nwp, discord,
					m, MHp, eta,
					stats,
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
    if(stats){
      for (j=0; j<m->n_stats; j++){
	stats[j+m->n_stats] = stats[j];
      }
      stats += m->n_stats;
    }

    /* This then adds the change statistics to these values */
    for(j=0;j<interval;j++){
      MCMCDynStatus status = MCMCDyn1Step(nwp, discord,
					  m, MHp, eta,
					  stats,
					  maxchanges, &nextdiffedge, difftime, difftail, diffhead, diffto,
					  min_MH_interval, max_MH_interval, MH_pval, MH_interval_add, fVerbose);
      
      // Check that we didn't run out of log space.
      if(status==MCMCDyn_TOO_MANY_CHANGES)
        return MCMCDyn_TOO_MANY_CHANGES;
      
      // If we need to return a network, then stop right there, since the network is too big to return, so stop early.
      if(maxedges!=0 && nwp->nedges >= maxedges-1)
	return MCMCDyn_TOO_MANY_EDGES;

      // Periodically expire older timestamps.
      if(nwp->duration_info->time%TIMESTAMP_HORIZON_FREQ==0) ExpireTimestamps(TIMESTAMP_HORIZON_EDGE, TIMESTAMP_HORIZON_NONEDGE, nwp);
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


/* The following outlines what happens during each time step:

   1. All statistics with x_ functions are sent a TICK signal.
   2. Any updates to statistics are recorded.
   3. Timer is incremented.
   4. A toggle proposal is made.
   5. Change statistics are calculated.
   6. If proposal is rejected, continue to 10.
   7. Relevant dyad is toggled in nwp.
   8. If it hasn't been toggled during this time step, its lasttoggle information, if any, is saved in discord, and its lasttoggle information is updated.
   9. If it has been toggled during this time step, its lasttoggle information, if any, is restored from discord. Discord's backup is deleted.
   10. Continue from 4 until sufficiently long run.
   11. All statistics with x_ functions are sent a TOCK signal.
   12. Any updates to statistics are recorded.
   13. Toggles in discord are logged.
   14. discord is cleared.

 */


/*********************
 void MCMCDyn1Step

 Simulate evolution of a dynamic network for 1 step.
*********************/
MCMCDynStatus MCMCDyn1Step(// Observed and discordant network.
                           Network *nwp, StoreDyadMapInt *discord,
		  // terms and proposals.
		  Model *m, MHProposal *MHp, double *eta,
		  // Space for output.
		  double *stats,
		  unsigned int maxchanges, Edge *nextdiffedge,
		  Vertex *difftime, Vertex *difftail, Vertex *diffhead, int *diffto,
		  // MCMC settings.
		  unsigned int min_MH_interval, unsigned int max_MH_interval, double MH_pval, double MH_interval_add,
		  // Verbosity.
		  int fVerbose){

  /* If the term has an extension, send it a "TICK" signal. */
  memset(m->workspace, 0, m->n_stats*sizeof(double)); /* Zero all change stats. */
  EXEC_THROUGH_TERMS_INTO(m, m->workspace, {
      mtp->dstats = dstats; /* Stuck the change statistic here.*/
      if(mtp->x_func)
        (*(mtp->x_func))(TICK, NULL, mtp, nwp);
    });
  /* Record network statistics for posterity. */
  if(stats) {
    for (unsigned int i = 0; i < m->n_stats; i++)
      stats[i] += m->workspace[i];
  }

  nwp->duration_info->time++; // Advance the clock.

  /* Run the process. */
  
  double cutoff;
  double 
    si = 0, // sum of increments
    si2 = 0, // sum of squared increments
    s = 0, // sum of weights 
    s2 = 0 // sum of squared weights
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
    double ip = 0;
    for (unsigned int i=0; i<m->n_stats; i++){
      ip += eta[i] * m->workspace[i];
      //  Rprintf("%f ", m->workspace[i]); 
    }
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
        GET_EDGE_UPDATE_STORAGE_TOGGLE(MHp->toggletail[i], MHp->togglehead[i], nwp, m, MHp);

        {
          TailHead dyad = THKey(discord,MHp->toggletail[i], MHp->togglehead[i]);
          khint_t i = kh_get(DyadMapInt,discord,dyad);
          if(i == kh_none){
            // If the dyad has *not* been toggled in this time step, then save its last toggle info (if any) and make a provisional change in lasttoggle.
            // Here, current time is used as a placeholder for no last toggle info.
            kh_set(DyadMapInt, discord, dyad, kh_getval(DyadMapInt, nwp->duration_info->lasttoggle, dyad, nwp->duration_info->time));
            kh_set(DyadMapInt, nwp->duration_info->lasttoggle, dyad, nwp->duration_info->time);
          }else{
            // If the dyad *has* been toggled in this timestep, then untoggle it by restoring its change in lasttoggle.
            if(kh_value(discord, i) != nwp->duration_info->time){
              kh_set(DyadMapInt, nwp->duration_info->lasttoggle, dyad, kh_value(discord, i));
            }else{
              kh_unset(DyadMapInt, nwp->duration_info->lasttoggle, dyad);
            }

            kh_del(DyadMapInt, discord, i);
          }
        }
      }
      /* Record network statistics for posterity. */
      if(stats) {
        for (unsigned int i = 0; i < m->n_stats; i++)
          stats[i] += m->workspace[i];
      }
    }

    int i = kh_size(discord) - prev_discord;
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

      if(fVerbose>=5) Rprintf("%u: s=%2.2f s2=%2.2f d=%d i=%d si=%2.2f si2=%2.2f mi=%2.2f vi=%2.2f ni=%2.2f zi=%2.2f pi=%2.2f\n", step, s, s2, kh_size(discord), i, si, si2, mi, vi, (s*s)/s2, zi, pi);
  
      if(pi > MH_pval){
	extrasteps = step*MH_interval_add+round(runif(0,1));
	finished++;
      }
    }

    if(finished) finished++;
  }

  /* Step finished: record changes. */
  
  if(fVerbose>=4){
    if(step>=max_MH_interval ) Rprintf("Convergence not achieved after %u M-H steps.\n",step);
    else Rprintf("Convergence achieved after %u M-H steps.\n",step);
  }

  return MCMCDyn1Step_advance(nwp, discord, m, stats,
                              maxchanges, nextdiffedge, difftime, difftail, diffhead, diffto,
                              fVerbose);
}


MCMCDynStatus MCMCDyn1Step_advance(// Observed and discordant network.
                                   Network *nwp, StoreDyadMapInt *discord,
                                   // terms and proposals.
                                   Model *m,
                                   // Space for output.
                                   double *stats,
                                   unsigned int maxchanges, Edge *nextdiffedge,
                                   Vertex *difftime, Vertex *difftail, Vertex *diffhead, int *diffto,
                                   // Verbosity.
                                   int fVerbose){

  /* If the term has an extension, send it a "TOCK" signal and the set
     of dyads that changed. */
  memset(m->workspace, 0, m->n_stats*sizeof(double)); /* Zero all change stats. */
  EXEC_THROUGH_TERMS_INTO(m, m->workspace, {
      mtp->dstats = dstats; /* Stuck the change statistic here.*/
      if(mtp->x_func)
        (*(mtp->x_func))(TOCK, discord, mtp, nwp);
    });
  /* Record network statistics for posterity. */
  if(stats) {
    for (unsigned int i = 0; i < m->n_stats; i++)
      stats[i] += m->workspace[i];
  }

  const unsigned int t=nwp->duration_info->time;

  TailHead dyad;
  kh_foreach_key(discord, dyad,{    
      if(*nextdiffedge<maxchanges){
        // and record the toggle.
        if(difftime) difftime[*nextdiffedge] = t;
        if(difftail) difftail[*nextdiffedge] = dyad.tail;
        if(diffhead) diffhead[*nextdiffedge] = dyad.head;
        if(diffto) diffto[*nextdiffedge] = GetEdge(dyad.tail,dyad.head,nwp);
        (*nextdiffedge)++;
      }else{
        return(MCMCDyn_TOO_MANY_CHANGES);
      }
    });
  kh_clear(DyadMapInt, discord);

  return MCMCDyn_OK;
}
