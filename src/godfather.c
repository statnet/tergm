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
  Network *nwp;
  Model *m;
  StoreDyadMapInt *discord;

  MCMCDyn_init_common(tails, heads, *time, lasttoggle, *n_edges,
              *n_nodes, *directed_flag, *bipartite, &nwp,
              *nterms, *funnames, *sonames, inputs, &m,
              NULL, NULL, NULL, NULL,
              NULL, 0, 0,
              NULL, NULL, NULL,
              &discord,
              *fVerbose);
  
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
    for(unsigned int i = 0; i < m->n_stats; i++)
      changestats[i] += m->workspace[i];

    nwp->duration_info->time++; // Advance the clock.

    while(pos < *total_toggles && toggletimes[pos]==t_stat+1){
      ChangeStats(1, (Vertex *)(toggletails+pos), (Vertex *)(toggleheads+pos), nwp, m);
    
      GET_EDGE_UPDATE_STORAGE_TOGGLE(toggletails[pos], toggleheads[pos], nwp, m, NULL);

      {
        TailHead dyad = THKey(discord, toggletails[pos], toggleheads[pos]);
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
        
      for(unsigned int i = 0; i < m->n_stats; i++)
        changestats[i] += m->workspace[i];

      pos++;
    }
    
    MCMCDyn1Step_advance(nwp, discord, m, changestats,
                         0, NULL, NULL, NULL, NULL, NULL,
                         *fVerbose);
                         
    if(nwp->duration_info->time%TIMESTAMP_HORIZON_FREQ==0) ExpireTimestamps(TIMESTAMP_HORIZON_EDGE, TIMESTAMP_HORIZON_NONEDGE, nwp);
  }

  if(*maxedges!=0 && nwp->nedges >= *maxedges-1){
    *status = MCMCDyn_TOO_MANY_EDGES;
  }

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

  /* Clean up and return */
  MCMCDyn_finish_common(nwp, m, NULL, discord);
}
