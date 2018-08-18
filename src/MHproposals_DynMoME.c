/*  File src/MHproposals_DynMoME.c in package tergm, part of the Statnet suite
 *  of packages for network analysis, http://statnet.org .
 *
 *  This software is distributed under the GPL-3 license.  It is free,
 *  open source, and has the attribution requirements (GPL Section 7) at
 *  http://statnet.org/attribution
 *
 *  Copyright 2008-2018 Statnet Commons
 */
#include "MHproposals_DynMoME.h"

/* Shorthand. */
#define Mtail (MHp->toggletail)
#define Mhead (MHp->togglehead)

/********************
   void MH_Formation
   Propose ONLY edges not in the reference graph
***********************/
void MH_Formation (MHproposal *MHp, Network *nwp) 
{  
  static Dyad ndyads;

  if(MHp->ntoggles == 0) { /* Initialize */
    MHp->ntoggles=1;
    ndyads = DYADCOUNT(nwp[0].nnodes, 0, nwp[0].directed_flag);
    return;
  }
  
  if(nwp[0].nedges==ndyads && nwp[1].nedges==0){ /* Attempting formation on a complete graph. */
    Mtail[0]=MH_FAILED;
    Mhead[0]=MH_IMPOSSIBLE;
    return;
  }
  
  BD_LOOP({
      /* A dyad eligible to be formed is either a nontie in nwp[0] or a tie in nwp[1]. */
      if(unif_rand() < ((double)nwp[1].nedges)/(ndyads-nwp[0].nedges + nwp[1].nedges)){
	// Tie in nwp[1]:
	GetRandEdge(Mtail, Mhead, nwp+1);
      }else{
	// Nontie in nwp[0]:
	GetRandNonedge(Mtail, Mhead, nwp);
      }
    });
}

/********************
   void MH_FormationTNT
   Tie/no tie:  Gives at least 50% chance of
   proposing a toggle of an existing (non-reference) edge, as opposed
   to simple random toggles that rarely do so in sparse 
   networks
   Propose ONLY edges not in the reference graph
***********************/
void MH_FormationTNT (MHproposal *MHp, Network *nwp) 
{  
  Edge nedges, ndedges;
  Dyad nempty;
  static double comp=0.5, odds;
  static Dyad ndyads;
  
  if(MHp->ntoggles == 0) { /* Initialize */
    MHp->ntoggles=1;
    odds = comp/(1.0-comp);
    ndyads = DYADCOUNT(nwp[0].nnodes, nwp->bipartite, nwp[0].directed_flag);
    return;
  }
  
  nedges  = nwp[0].nedges;
  ndedges = nwp[1].nedges;
  nempty = ndyads-nedges;

  if(nempty==0 && ndedges==0){  /* Attempting formation on a complete graph. */
    Mtail[0]=MH_FAILED;
    Mhead[0]=MH_IMPOSSIBLE;
    return;
  }

  double logratio=0;
  BD_LOOP({
      if(ndedges != 0 && (nempty == 0 || unif_rand() < comp)) { /* Select a discordant dyad at random */
	GetRandEdge(Mtail, Mhead, nwp+1);
	
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
	    logratio = log((double) nempty / (ndedges+1) / (1-comp));
	  }else{
	    logratio = log((double) nempty / (ndedges+1) * odds);
	  }
	}
      }
    });
  MHp->logratio += logratio;
  //   Rprintf("nedges %d reference nddyads %d h %d t %d MHp->ratio %f\n", 
  //	    nwp[0].nedges, nwp[1].nedges, tail, head, MHp->ratio); 
}

/********************
   void MH_Dissolution
   Propose ONLY edges in the reference graph
   Any candidate edge must be either in nwp[0] or in nwp[1]. This makes
   proposing easy.
***********************/
void MH_Dissolution (MHproposal *MHp, Network *nwp) 
{  
  Edge nedges=nwp[0].nedges, ndedges=nwp[1].nedges;
  
  if(MHp->ntoggles == 0) { /* Initialize */
    MHp->ntoggles=1;
    return;
  }
  
  if(nedges==0 && ndedges==0){ /* Attempting dissolution on an empty network. */
    Mtail[0]=MH_FAILED;
    Mhead[0]=MH_IMPOSSIBLE;
    return;
  }
  
  /* Make sure the edge is not in nwp[0] and nwp[1], which shouldn't
     happen anyway, but better safe than sorry. */
  BD_LOOP({
      Edge rane = 1 + unif_rand() * (nedges+ndedges);
      if(rane<=nedges){
	GetRandEdge(Mtail, Mhead, nwp);
	if(EdgetreeSearch(Mtail[0],Mhead[0],nwp[1].outedges)==0 && CheckTogglesValid(MHp, nwp)) break;
      }else{
	GetRandEdge(Mtail, Mhead, nwp+1);
	if(EdgetreeSearch(Mtail[0],Mhead[0],nwp[0].outedges)==0 && CheckTogglesValid(MHp, nwp)) break;
      }
    });
}

/********************
   void MH_DissolutionTNT
   Tie/no tie:  Gives at least 50% chance of
   proposing a toggle of an existing (non-reference) edge, as opposed
   to simple random toggles that rarely do so in sparse 
   networks
   Propose ONLY edges in the reference graph
***********************/
void MH_DissolutionTNT (MHproposal *MHp, Network *nwp) 
{  
  Edge nedges=nwp[0].nedges, ndedges=nwp[1].nedges;
  Edge nedges0 = nedges + ndedges;
  static double comp=0.5, odds;
  
  if(MHp->ntoggles == 0) { /* Initialize */
    MHp->ntoggles=1;
    odds = comp/(1.0-comp);
    return;
  }

  if(nedges==0 && ndedges==0){  /* Attempting dissolution on an empty graph. */
    Mtail[0]=MH_FAILED;
    Mhead[0]=MH_IMPOSSIBLE;
    return;
  }

  double logratio=0;
  BD_LOOP({
      if(ndedges != 0 && (nedges==0 || unif_rand() < comp)) { /* Select a discordant dyad at random */
	GetRandEdge(Mtail, Mhead, nwp+1);
	
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
  //   Rprintf("nedges %d reference nddyads %d h %d t %d MHp->ratio %f\n", 
  //	    nwp[0].nedges, nwp[1].nedges, tail, head, MHp->ratio); 
}
