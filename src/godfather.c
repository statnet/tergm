/*  File src/godfather.c in package tergm, part of the Statnet suite
 *  of packages for network analysis, https://statnet.org .
 *
 *  This software is distributed under the GPL-3 license.  It is free,
 *  open source, and has the attribution requirements (GPL Section 7) at
 *  https://statnet.org/attribution
 *
 *  Copyright 2008-2019 Statnet Commons
 */
#include "godfather.h"

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
void godfather_wrapper(int *tails, int *heads, int *time, int *lasttoggle, int *n_edges,
               int *n_nodes, int *directed_flag, int *bipartite, 
               int *nterms, char **funnames, char **sonames, double *inputs,
               int *total_toggles, int *toggletimes, 
               int *toggletails, int *toggleheads,
               int *start_time, int *end_time,
               double *changestats, 
               int *maxedges,
               int *newnetworktails, 
               int *newnetworkheads, 
               int *fVerbose, 
               int *status){
  ErgmState *s;
  StoreDyadMapInt *discord;

  MCMCDyn_init_common(tails, heads, *time, lasttoggle, *n_edges,
              *n_nodes, *directed_flag, *bipartite,
              *nterms, *funnames, *sonames, inputs,
              NULL, NULL, NULL, NULL,
              NULL, 0, 0,
              NULL, NULL,
              &s, &discord);
  
  Network *nwp = s->nwp;
  Model *m = s->m;  
  
  /*********************
  changestats are modified in groups of m->n_stats, and they
  reflect the CHANGE in the values of the statistics from the
  original network.  Thus, when we begin, the initial 
  values of the first group of m->n_stats changestats should 
  all be zero
  *********************/
  
  memset(changestats, 0, m->n_stats*sizeof(double));
  
  /* Now start obtaining change statistics */

  unsigned int pos = 0;
  // The reason it's = start_time but < end_time is that by the time t_stat arrives at end_time, the toggles for end_time will have already been applied.
  for(unsigned int t_stat = *start_time; t_stat < *end_time; t_stat++){
    changestats += m->n_stats;
    memcpy(changestats, changestats-m->n_stats, m->n_stats*sizeof(double));
    
    /* If the term has an extension, send it a "TICK" signal. */
    memset(m->workspace, 0, m->n_stats*sizeof(double)); /* Zero all change stats. */
    SIGNAL_TERMS_INTO(m, m->workspace, TICK, NULL);
    /* Record network statistics for posterity. */
    addonto(changestats, m->workspace, m->n_stats);

    while(pos < *total_toggles && toggletimes[pos]==t_stat+1){
      ChangeStats(1, (Vertex *)(toggletails+pos), (Vertex *)(toggleheads+pos), nwp, m);
    
      GET_EDGE_UPDATE_STORAGE_TOGGLE(toggletails[pos], toggleheads[pos], nwp, m, NULL);

    }
        
    addonto(changestats, m->workspace, m->n_stats);

    pos++;
  }
    
    MCMCDyn1Step_advance(s, discord, changestats,
                         0, NULL, NULL, NULL, NULL, NULL,
                         *fVerbose);
                         
  }

  if(*maxedges!=0 && nwp->nedges >= *maxedges-1){
    *status = MCMCDyn_TOO_MANY_EDGES;
  }

  if(*status == MCMCDyn_OK && *maxedges>0){
    newnetworktails[0]=newnetworkheads[0]=EdgeTree2EdgeList((Vertex*)newnetworktails+1,(Vertex*)newnetworkheads+1,nwp,*maxedges-1);
  }

  /* Clean up and return */
  MCMCDyn_finish_common(s, discord);
}
