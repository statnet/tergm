/*  File src/MHproposals_DynMLE.c in package tergm, part of the Statnet suite
 *  of packages for network analysis, http://statnet.org .
 *
 *  This software is distributed under the GPL-3 license.  It is free,
 *  open source, and has the attribution requirements (GPL Section 7) at
 *  http://statnet.org/attribution
 *
 *  Copyright 2003-2013 Statnet Commons
 */
#include "MHproposals_DynMLE.h"
#include "edgelist.h"

/* Shorthand. */
#define Mtail (MHp->toggletail)
#define Mhead (MHp->togglehead)

/********************
   void MH_FormationMLE
   Propose ONLY edges not in the reference graph
***********************/
void MH_FormationMLE (MHproposal *MHp, Network *nwp) 
{  
  static Vertex nnodes;
  static Edge ndyads, nedges0;

  if(MHp->ntoggles == 0) { /* Initialize */
    MHp->ntoggles=1;
    nnodes = nwp[0].nnodes;
    ndyads = DYADCOUNT(nnodes, 0, nwp[0].directed_flag);
    nedges0 = MHp->inputs[0];
    return;
  }
  
  if(nedges0==ndyads){ /* Attempting formation on a complete graph. */
    Mtail[0]=MH_FAILED;
    Mhead[0]=MH_IMPOSSIBLE;
    return;
  }

  BD_COND_LOOP({
    /* Keep trying dyads until a one that is not an edge in the reference network is found. */
    /* Generate. */
      GetRandDyad(Mtail, Mhead, nwp);
    }, !dEdgeListSearch(Mtail[0],Mhead[0],MHp->inputs), 2);
}

/********************
   void MH_FormationMLE
   Propose ONLY edges not in the reference graph
***********************/
void MH_FormationMLETNT(MHproposal *MHp, Network *nwp) 
{  
  static Vertex nnodes;
  Vertex tail,head;
  static Edge ndyads, nedges0;
  static double comp=0.5, odds;
  static Network discord;

  if(MHp->ntoggles == 0) { /* Initialize */
    MHp->ntoggles=1;
    nnodes = nwp->nnodes;
    odds = comp/(1.0-comp);
    ndyads = DYADCOUNT(nnodes, nwp->bipartite, nwp[0].directed_flag);

    nedges0 = MHp->inputs[0];
    MHp->discord = (Network**) calloc(2,sizeof(Network*)); // A space for the sentinel NULL pointer.
    MHp->discord[0] = &discord;
    discord = NetworkInitializeD(MHp->inputs+1, MHp->inputs+1+nedges0, nedges0, nnodes, nwp->directed_flag, nwp->bipartite, 0, 0, NULL);
   
    for(Edge i=0; i<nwp->nedges; i++){
      FindithEdge(&tail, &head, i+1, nwp);
      ToggleEdge(tail,head, &discord);
    }
    
    return;
  }
  

  Edge nedges = nwp->nedges, ndedges = discord.nedges, nempty = ndyads-nedges;

  if(nempty==0 && ndedges==0){ /* Attempting formation on a complete graph. */
    Mtail[0]=MH_FAILED;
    Mhead[0]=MH_IMPOSSIBLE;
    return;
  }

  double logratio=0;
  BD_LOOP({
      if(ndedges != 0 && (nempty == 0 || unif_rand() < comp)) { /* Select a discordant dyad at random */
	GetRandEdge(Mtail, Mhead, &discord);
	
	if(nempty==0){
	  logratio = log(ndedges*(1-comp));
	}else{
	  if(ndedges==1){
	    logratio = log(1.0 / (nempty + 1) /comp);
	  }else{
	    logratio = log((double) ndedges / (nempty + 1) / odds);
	  }
	}
      }else{ /* select an empty dyad in nwp[0] at random */
	GetRandNonedge(Mtail, Mhead, nwp);
	
	if(ndedges==0){
	  logratio = log(nempty*comp);
	}else{
	  if(nempty==1){
	    logratio = log(1.0 / ndyads / (1-comp));
	  }else{
	    logratio = log((double) nempty / (ndedges+1) * odds);
	  }
	}
      }
    });

  MHp->logratio += logratio;
}

/********************
   void MH_DissolutionMLE
   Propose ONLY edges in the reference graph
***********************/
void MH_DissolutionMLETNT(MHproposal *MHp, Network *nwp) 
{  
  static Vertex nnodes;
  Vertex tail,head;
  static double comp=0.5, odds;
  static Network discord;
  static Edge nedges0;

  if(MHp->ntoggles == 0) { /* Initialize */
    MHp->ntoggles=1;
    nnodes = nwp->nnodes;
    odds = comp/(1.0-comp);

    nedges0 = MHp->inputs[0];
    MHp->discord = (Network**) calloc(2,sizeof(Network*)); // A space for the sentinel NULL pointer.
    MHp->discord[0] = &discord;
    discord = NetworkInitializeD(MHp->inputs+1, MHp->inputs+1+nedges0, nedges0, nnodes, nwp->directed_flag, nwp->bipartite, 0, 0, NULL);
   
    for(Edge i=0; i<nwp->nedges; i++){
      FindithEdge(&tail, &head, i+1, nwp);
      ToggleEdge(tail,head, &discord);
    }
    
    return;
  }
  

  Edge nedges = nwp->nedges, ndedges = discord.nedges;

  if(nedges==0 && ndedges==0){ /* Attempting dissolution on an empty graph. */
    Mtail[0]=MH_FAILED;
    Mhead[0]=MH_IMPOSSIBLE;
    return;
  }
  
  double logratio=0;
  BD_LOOP({
      if(ndedges != 0 && (nedges == 0 || unif_rand() < comp)) { /* Select a discordant dyad at random */
	GetRandEdge(Mtail, Mhead, &discord);
	
	if(nedges==0){
	  logratio = log(ndedges*(1-comp));
	}else{
	  if(ndedges==1){
	    logratio = log((double) ndedges / nedges0 / comp);
	  }else{
	    logratio = log((double) ndedges / (nedges + 1) / odds);
	  }
	}
      }else{ /* select an edge in nwp[0] at random */
	GetRandEdge(Mtail, Mhead, nwp);
	
	if(ndedges==0){
	  logratio = log(nedges*comp);
	}else{
	  if(nedges==1){
	    logratio = log(1.0 / nedges0 / (1-comp));
	  }else{
	    logratio = log((double) nedges / (ndedges+1) * odds);
	  }
	}
      }
    });
  MHp->logratio += logratio;
}
