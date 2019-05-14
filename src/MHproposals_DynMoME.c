/*  File src/MHproposals_DynMoME.c in package tergm, part of the Statnet suite
 *  of packages for network analysis, https://statnet.org .
 *
 *  This software is distributed under the GPL-3 license.  It is free,
 *  open source, and has the attribution requirements (GPL Section 7) at
 *  https://statnet.org/attribution
 *
 *  Copyright 2008-2019 Statnet Commons
 */
#include "MHproposals_DynMoME.h"

/* Shorthand. */
#define Mtail (MHp->toggletail)
#define Mhead (MHp->togglehead)

/********************
   void MH_Formation
   Propose ONLY edges not in the reference graph
***********************/
void MH_Formation (MHProposal *MHp, Network *nwp) 
{  
  static Dyad ndyads;

  if(MHp->ntoggles == 0) { /* Initialize */
    MHp->ntoggles=1;
    ndyads = DYADCOUNT(nwp->nnodes, 0, nwp->directed_flag);
    return;
  }
  
  if(nwp->nedges==ndyads && MHp->discord[0]->nedges==0){ /* Attempting formation on a complete graph. */
    Mtail[0]=MH_FAILED;
    Mhead[0]=MH_IMPOSSIBLE;
    return;
  }
  
  BD_LOOP({
      /* A dyad eligible to be formed is either a nontie in nwp[0] or a tie in MHp->discord-> */
      if(unif_rand() < ((double)MHp->discord[0]->nedges)/(ndyads-nwp->nedges + MHp->discord[0]->nedges)){
	// Tie in MHp->discord[0]:
	GetRandEdge(Mtail, Mhead, MHp->discord[0]);
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
void MH_FormationTNT (MHProposal *MHp, Network *nwp) 
{  
  static double comp=0.5, odds;
  static Dyad ndyads;
  
  if(MHp->ntoggles == 0) { /* Initialize */
    MHp->ntoggles=1;
    odds = comp/(1.0-comp);
    ndyads = DYADCOUNT(nwp->nnodes, nwp->bipartite, nwp->directed_flag);
    return;
  }
  
  Edge nedges  = nwp->nedges, ndedges = MHp->discord[0]->nedges;
  Dyad nempty = ndyads-nedges;

  if(nempty==0 && ndedges==0){  /* Attempting formation on a complete graph. */
    Mtail[0]=MH_FAILED;
    Mhead[0]=MH_IMPOSSIBLE;
    return;
  }

  double logratio=0;
  BD_LOOP({
      if(ndedges != 0 && (nempty == 0 || unif_rand() < comp)) { /* Select a discordant dyad at random */
	GetRandEdge(Mtail, Mhead, MHp->discord[0]);
	
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
  //	    nwp->nedges, MHp->discord[0]->nedges, tail, head, MHp->ratio); 
}

/********************
   void MH_Dissolution
   Propose ONLY edges in the reference graph
   Any candidate edge must be either in nwp[0] or in MHp->discord[0]-> This makes
   proposing easy.
***********************/
void MH_Dissolution (MHProposal *MHp, Network *nwp) 
{  
  if(MHp->ntoggles == 0) { /* Initialize */
    MHp->ntoggles=1;
    return;
  }

  Edge nedges=nwp->nedges, ndedges=MHp->discord[0]->nedges;
  
  if(nedges==0 && ndedges==0){ /* Attempting dissolution on an empty network. */
    Mtail[0]=MH_FAILED;
    Mhead[0]=MH_IMPOSSIBLE;
    return;
  }
  
  /* Make sure the edge is not in nwp[0] and MHp->discord[0], which shouldn't
     happen anyway, but better safe than sorry. */
  BD_LOOP({
      Edge rane = 1 + unif_rand() * (nedges+ndedges);
      if(rane<=nedges){
	GetRandEdge(Mtail, Mhead, nwp);
	if(EdgetreeSearch(Mtail[0],Mhead[0],MHp->discord[0]->outedges)==0 && CheckTogglesValid(MHp, nwp)) break;
      }else{
	GetRandEdge(Mtail, Mhead, MHp->discord[0]);
	if(EdgetreeSearch(Mtail[0],Mhead[0],nwp->outedges)==0 && CheckTogglesValid(MHp, nwp)) break;
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
void MH_DissolutionTNT (MHProposal *MHp, Network *nwp) 
{  
  static double comp=0.5, odds;
  
  if(MHp->ntoggles == 0) { /* Initialize */
    MHp->ntoggles=1;
    odds = comp/(1.0-comp);
    return;
  }

  Edge nedges=nwp->nedges, ndedges=MHp->discord[0]->nedges;
  Edge nedges0 = nedges + ndedges;
  
  if(nedges==0 && ndedges==0){  /* Attempting dissolution on an empty graph. */
    Mtail[0]=MH_FAILED;
    Mhead[0]=MH_IMPOSSIBLE;
    return;
  }

  double logratio=0;
  BD_LOOP({
      if(ndedges != 0 && (nedges==0 || unif_rand() < comp)) { /* Select a discordant dyad at random */
	GetRandEdge(Mtail, Mhead, MHp->discord[0]);
	
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
  //	    nwp->nedges, MHp->discord[0]->nedges, tail, head, MHp->ratio); 
}
