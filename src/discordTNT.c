/*  File src/discordTNT.c in package tergm, part of the Statnet suite
 *  of packages for network analysis, https://statnet.org .
 *
 *  This software is distributed under the GPL-3 license.  It is free,
 *  open source, and has the attribution requirements (GPL Section 7) at
 *  https://statnet.org/attribution
 *
 *  Copyright 2008-2020 Statnet Commons
 */
#include "ergm_MHproposal.h"
#include "ergm_edgelist.h"
#include "ergm_changestat.h"
#include "ergm_rlebdm.h"
#include "ergm_MHstorage.h"
#include "ergm_unsorted_edgelist.h"
#include "ergm_weighted_population.h"
#include "changestats_lasttoggle.h"
#include "tergm_model.h"
#include "ergm_Rutil.h"
#include "ergm_nodelist_dyad_sampler.h"
#include "ergm_BDStrat_proposals.h"

#define OUTVAL_NET(e,n) ((n)->outedges[(e)].value)
#define INVAL_NET(e,n) ((n)->inedges[(e)].value)
#define MIN_OUTEDGE_NET(a,n) (EdgetreeMinimum((n)->outedges, (a)))
#define MIN_INEDGE_NET(a,n) (EdgetreeMinimum((n)->inedges, (a)))
#define NEXT_OUTEDGE_NET(e,n) (EdgetreeSuccessor((n)->outedges,(e)))
#define NEXT_INEDGE_NET(e,n) (EdgetreeSuccessor((n)->inedges,(e)))
#define STEP_THROUGH_OUTEDGES_NET(a,e,v,n) for((e)=MIN_OUTEDGE_NET((a),(n));((v)=OUTVAL_NET((e),(n)))!=0;(e)=NEXT_OUTEDGE_NET((e),(n)))
#define STEP_THROUGH_INEDGES_NET(a,e,v,n) for((e)=MIN_INEDGE_NET((a),(n));((v)=INVAL_NET((e),(n)))!=0;(e)=NEXT_INEDGE_NET((e),(n)))


typedef struct {
  UnsrtEL *discordantEdges;
  UnsrtEL *discordantNonEdges;
  
  double discordance_fraction;
  
  int in_discord;
  
} discordTNTStorage; 

MH_I_FN(Mi_discordTNT) {
  MHp->ntoggles = 1;
  
  ALLOC_STORAGE(1, discordTNTStorage, sto);
  sto->discordantEdges = UnsrtELInitialize(0, NULL, NULL, FALSE);
  sto->discordantNonEdges = UnsrtELInitialize(0, NULL, NULL, FALSE);

  sto->discordance_fraction = asReal(getListElement(MHp->R, "discordance_fraction"));  
  // we ignore discord for this initialization (assuming a TICK will precede any proposals)
}

MH_X_FN(Mx_discordTNT) {
  GET_STORAGE(discordTNTStorage, sto);
  
  if(type == TICK) {
    // "clear" the discordant dyads
    sto->discordantEdges->nedges = 0;
    sto->discordantNonEdges->nedges = 0;
  }
}

MH_P_FN(MH_discordTNT) {
  GET_AUX_STORAGE(StoreTimeAndLasttoggle, dur_inf);
  
  GET_STORAGE(discordTNTStorage, sto);
  
  int in_discord;
  int in_network;
  
  int nedges = EDGECOUNT(nwp);
  int nddyads = kh_size(dur_inf->discord);
  
  if(nddyads == 0 || unif_rand() < 1 - sto->discordance_fraction) {
    // propose from network
    if(nedges == 0 || unif_rand() < 0.5) {
      // propose toggling a random dyad
      GetRandDyad(Mtail, Mhead, nwp);
      in_network = IS_OUTEDGE(Mtail[0], Mhead[0]);
      in_discord = kh_get(DyadMapInt, dur_inf->discord, THKey(dur_inf->discord, Mtail[0], Mhead[0])) != kh_none;
    } else {
      // propose toggling off an edge in network
      GetRandEdge(Mtail, Mhead, nwp);
      in_network = TRUE;
      in_discord = kh_get(DyadMapInt, dur_inf->discord, THKey(dur_inf->discord, Mtail[0], Mhead[0])) != kh_none;      
    }

    if(in_discord) {
      // need to resample to know index
      if(in_network) {
        UnsrtELGetRand(Mtail, Mhead, sto->discordantEdges);
      } else {
        UnsrtELGetRand(Mtail, Mhead, sto->discordantNonEdges);
      }
    }
  } else {
    // propose from discord
    if(unif_rand() < sto->discordantEdges->nedges/((double) nddyads)) {
      UnsrtELGetRand(Mtail, Mhead, sto->discordantEdges);
      in_network = TRUE;
    } else {
      UnsrtELGetRand(Mtail, Mhead, sto->discordantNonEdges);
      in_network = FALSE;
    }
    in_discord = TRUE;
  }
  
  // compute logratio
  
  Dyad ndyads = DYADCOUNT(nwp);
  
  double forward_discord = in_discord ? 1.0/nddyads : 0;
  double backward_discord = in_discord ? 0 : 1.0/(1 + nddyads);
  
  double forward_network = in_network ? (0.5/nedges + 0.5/ndyads) : (nedges == 0 ? 1.0/ndyads : 0.5/ndyads);
  double backward_network = in_network ? (nedges == 1 ? 1.0/ndyads : 0.5/ndyads) : (0.5/(nedges + 1) + 0.5/ndyads);
    
  double forward = (nddyads == 0) ? forward_network : (sto->discordance_fraction*forward_discord + (1 - sto->discordance_fraction)*forward_network);
  double backward = (nddyads == 1 && in_discord) ? backward_network : (sto->discordance_fraction*backward_discord + (1 - sto->discordance_fraction)*backward_network);

  MHp->logratio = log(backward/forward);

  sto->in_discord = in_discord;
}

MH_U_FN(Mu_discordTNT) {
  // add or remove the dyad from the appropriate discordance edgelist
  
  GET_STORAGE(discordTNTStorage, sto);
  
  if(sto->in_discord) {
    // currently in discord; take it out
    if(edgeflag) {        
      UnsrtELDelete(tail, head, sto->discordantEdges);
    } else {
      UnsrtELDelete(tail, head, sto->discordantNonEdges);  
    }
  } else {
    // not currently in discord; add it in
    if(edgeflag) {
      UnsrtELInsert(tail, head, sto->discordantNonEdges);
    } else {
      UnsrtELInsert(tail, head, sto->discordantEdges);        
    }
  }
}

MH_F_FN(Mf_discordTNT) {
  GET_STORAGE(discordTNTStorage, sto);
  
  UnsrtELDestroy(sto->discordantNonEdges);
  UnsrtELDestroy(sto->discordantEdges);
}



/********************
   discordStratTNT
********************/

typedef struct {
  UnsrtEL **nonDiscordantELs;
  UnsrtEL **discordantELs;
  UnsrtEL **discordantNonELs;
  
  Vertex **nodesbycode;
  
  WtPop *wtp;
  
  int *nodecountsbycode;
  Dyad *dyadcounts;
  
  int *tailtypes;
  int *headtypes;
  
  int mixingtype;
  
  int nmixtypes;
  
  int in_discord;

  StratTNTStorage *static_sto;
} discordStratTNTStorage; 


MH_I_FN(Mi_discordStratTNT) {
  // let StratTNT's I_FN do most of the work
  Mi_StratTNT(MHp, nwp);

  // clean up from static proposal initialization and handle purely temporal aspects
  StratTNTStorage *static_sto = MH_STORAGE;
  ALLOC_STORAGE(1, discordStratTNTStorage, sto);
  sto->nonDiscordantELs = static_sto->els;
  sto->nodesbycode = static_sto->nodesbycode;
  sto->wtp = static_sto->wtp;
  sto->nodecountsbycode = static_sto->nodecountsbycode;
  sto->dyadcounts = static_sto->ndyadstype;
  sto->tailtypes = static_sto->tailtypes;
  sto->headtypes = static_sto->headtypes;
  sto->nmixtypes = static_sto->nmixtypes;
  sto->static_sto = static_sto;

  sto->discordantELs = Calloc(sto->nmixtypes, UnsrtEL *);
  sto->discordantNonELs = Calloc(sto->nmixtypes, UnsrtEL *);
  for(int i = 0; i < sto->nmixtypes; i++) {
    sto->discordantELs[i] = UnsrtELInitialize(0, NULL, NULL, FALSE);
    sto->discordantNonELs[i] = UnsrtELInitialize(0, NULL, NULL, FALSE);
  }
} 

MH_X_FN(Mx_discordStratTNT) {
  GET_STORAGE(discordStratTNTStorage, sto);
    
  if(type == TICK) {
    // transfer discordant edges to nondiscordant edges
    // clear discordant edges and discordant nonedges
    for(int i = 0; i < sto->nmixtypes; i++) {
      // UnsrtELs start at index 1 here
      for(int j = 1; j <= sto->discordantELs[i]->nedges; j++) {
        UnsrtELInsert(sto->discordantELs[i]->tails[j], sto->discordantELs[i]->heads[j], sto->nonDiscordantELs[i]);
      }
      
      sto->discordantELs[i]->nedges = 0;
      sto->discordantNonELs[i]->nedges = 0;      
    }
  }
}

MH_P_FN(MH_discordStratTNT) {
  GET_AUX_STORAGE(StoreTimeAndLasttoggle, dur_inf);
  GET_STORAGE(discordStratTNTStorage, sto);
        
  // sample a strat mixing type on which to make a proposal
  int i = WtPopGetRand(sto->wtp);
  
  // record the mixing type of the toggle, in case it's needed in the U function later
  sto->mixingtype = i;    
  
  int tailtype = sto->tailtypes[i];
  int headtype = sto->headtypes[i];
    
  // number of edges of this mixing type
  int nedgestype = sto->nonDiscordantELs[i]->nedges + sto->discordantELs[i]->nedges;

  // number of dyads of this mixing type
  Dyad ndyadstype = sto->dyadcounts[i];
  
  // number of discordant dyads of this mixing type
  int nddyadstype = sto->discordantNonELs[i]->nedges + sto->discordantELs[i]->nedges;
  
  // flags
  int in_discord;
  int in_network;
  
  if(nddyadstype == 0 || unif_rand() < 0.5) {
    // propose from network
    if(nedgestype == 0 || unif_rand() < 0.5) {
      // propose toggling a random dyad of the specified mixing type
      GetRandDyadFromLists(Mtail, // tail
                           Mhead, // head
                           sto->nodesbycode, // tails
                           sto->nodesbycode, // heads
                           &tailtype, // tailattrs
                           &headtype, // headattrs
                           sto->nodecountsbycode, // tailcounts
                           sto->nodecountsbycode, // headcounts
                           1, // length; only one allowed pairing since we've already sampled the strat type
                           ndyadstype, // dyadcount
                           TRUE, // diagonal; no higher level types here, so always TRUE
                           DIRECTED); // directed
      
      in_network = IS_OUTEDGE(Mtail[0], Mhead[0]);
      in_discord = kh_get(DyadMapInt, dur_inf->discord, THKey(dur_inf->discord, Mtail[0], Mhead[0])) != kh_none;
      
      // if it resides in any of the edgelists we store, we need to resample
      if(in_network) {
        if(in_discord) {
          UnsrtELGetRand(Mtail, Mhead, sto->discordantELs[i]);
        } else {
          UnsrtELGetRand(Mtail, Mhead, sto->nonDiscordantELs[i]);
        }
      } else if(in_discord) {
        UnsrtELGetRand(Mtail, Mhead, sto->discordantNonELs[i]);
      }
    } else {
      // propose toggling off an edge of the specified mixing type
      if(unif_rand() < sto->nonDiscordantELs[i]->nedges/((double) nedgestype)) {
        UnsrtELGetRand(Mtail, Mhead, sto->nonDiscordantELs[i]);
        in_discord = FALSE;
      } else {
        UnsrtELGetRand(Mtail, Mhead, sto->discordantELs[i]);
        in_discord = TRUE;
      }
      in_network = TRUE;
    }
  } else {
    // propose from discord
    if(unif_rand() < sto->discordantELs[i]->nedges/((double) nddyadstype)) {
      UnsrtELGetRand(Mtail, Mhead, sto->discordantELs[i]);
      in_network = TRUE;
    } else {
      UnsrtELGetRand(Mtail, Mhead, sto->discordantNonELs[i]);
      in_network = FALSE;
    }
    in_discord = TRUE;
  }
  
  // compute logratio
  
  // these ignore overall factor of 1/2 which cancels out in ratio
  // and overall factor of prob(mixing type), which also cancels out in ratio
  double forward_discord = in_discord ? 1.0/nddyadstype : 0;
  double backward_discord = in_discord ? 0 : 1.0/(1 + nddyadstype);
  
  double forward_network = in_network ? (0.5/nedgestype + 0.5/ndyadstype) : (nedgestype == 0 ? 1.0/ndyadstype : 0.5/ndyadstype);
  double backward_network = in_network ? (nedgestype == 1 ? 1.0/ndyadstype : 0.5/ndyadstype) : (0.5/(nedgestype + 1) + 0.5/ndyadstype);
  
  if(nddyadstype == 0) forward_network *= 2;
  if(nddyadstype == 1 && in_discord) backward_network *= 2;
  
  double forward = forward_discord + forward_network;
  double backward = backward_discord + backward_network;

  MHp->logratio = log(backward/forward);

  sto->in_discord = in_discord;
}

MH_U_FN(Mu_discordStratTNT) {
  // add or remove edge from appropriate edgelist
  GET_STORAGE(discordStratTNTStorage, sto);
  
  if(edgeflag) {
    // we are removing an existing edge
    if(sto->in_discord) {
      UnsrtELDelete(tail, head, sto->discordantELs[sto->mixingtype]);
    } else {
      UnsrtELDelete(tail, head, sto->nonDiscordantELs[sto->mixingtype]);
      UnsrtELInsert(tail, head, sto->discordantNonELs[sto->mixingtype]);
    }
  } else {
    // we are adding a new edge
    if(sto->in_discord) {
      UnsrtELInsert(tail, head, sto->nonDiscordantELs[sto->mixingtype]);
      UnsrtELDelete(tail, head, sto->discordantNonELs[sto->mixingtype]);
    } else {
      UnsrtELInsert(tail, head, sto->discordantELs[sto->mixingtype]);
    }
  }
}

MH_F_FN(Mf_discordStratTNT) {
  GET_STORAGE(discordStratTNTStorage, sto);

  // let StratTNT's F_FN do most of the work
  MH_STORAGE = sto->static_sto;
  Mf_StratTNT(MHp, nwp);
  MH_STORAGE = sto;
  Free(sto->static_sto);

  for(int i = 0; i < sto->nmixtypes; i++) {
    UnsrtELDestroy(sto->discordantELs[i]);
    UnsrtELDestroy(sto->discordantNonELs[i]);    
  }

  Free(sto->discordantELs);
  Free(sto->discordantNonELs);
  // MHp->storage itself should be Freed by MHProposalDestroy
}



/********************
    MH_discordBDTNT
********************/

typedef struct {
  int *attrcounts;
  Vertex **nodesvec;
  int *nodepos;
    
  int tailtype;
  int tailmaxl;
  
  int headtype;
  int headmaxl;
  
  Dyad currentdyads;
  Dyad proposeddyads;
  
  int currentsubmaxledges;
  int proposedsubmaxledges;  
  
  int bound;
  int nmixtypes;
  int *vattr;
  
  int *tailtypes;
  int *headtypes;

  UnsrtEL *BDTDNE;
  UnsrtEL *nonBDTDNE;

  UnsrtEL *nonDiscordantEdges;
  UnsrtEL *discordantEdges;
  
  Network *combined_BDTDNE;
  Network *combined_nonBDTDNE;
  UnsrtEL *transferEL;
  
  int in_discord;
  
  BDTNTStorage *static_sto;
} discordBDTNTStorage;

MH_I_FN(Mi_discordBDTNT) {
  // let BDTNT's I_FN do most of the work
  Mi_BDTNT(MHp, nwp);

  // clean up from static proposal initialization and handle purely temporal aspects
  BDTNTStorage *static_sto = MH_STORAGE;
  ALLOC_STORAGE(1, discordBDTNTStorage, sto);
  sto->attrcounts = static_sto->attrcounts;
  sto->nodesvec = static_sto->nodesvec;
  sto->nodepos = static_sto->nodepos;
  sto->currentdyads = static_sto->currentdyads;
  sto->currentsubmaxledges = static_sto->currentsubmaxledges;
  sto->bound = static_sto->bound;
  sto->nmixtypes = static_sto->nmixtypes;
  sto->vattr = static_sto->vattr;
  sto->tailtypes = static_sto->tailtypes;
  sto->headtypes = static_sto->headtypes;
  sto->nonDiscordantEdges = static_sto->edgelist;
  sto->static_sto = static_sto;

  sto->BDTDNE = UnsrtELInitialize(0, NULL, NULL, FALSE);
  sto->nonBDTDNE = UnsrtELInitialize(0, NULL, NULL, FALSE);
  sto->discordantEdges = UnsrtELInitialize(0, NULL, NULL, FALSE);
  sto->combined_BDTDNE = NetworkInitialize(NULL, NULL, 0, N_NODES, DIRECTED, BIPARTITE, FALSE, 0, NULL);
  sto->combined_nonBDTDNE = NetworkInitialize(NULL, NULL, 0, N_NODES, DIRECTED, BIPARTITE, FALSE, 0, NULL);
  sto->transferEL = UnsrtELInitialize(0, NULL, NULL, FALSE);
}

MH_X_FN(Mx_discordBDTNT) {
  if(type == TICK) {    
    GET_STORAGE(discordBDTNTStorage, sto);
    
    // transfer discordant edges to nondiscordant edges    
    for(int i = 1; i <= sto->discordantEdges->nedges; i++) {
      UnsrtELInsert(sto->discordantEdges->tails[i], sto->discordantEdges->heads[i], sto->nonDiscordantEdges);
    }
    
    // clear all the discordance information    
    sto->BDTDNE->nedges = 0;
    sto->nonBDTDNE->nedges = 0;
    sto->discordantEdges->nedges = 0;
    
    // for now, destroy and recreate each time step (can we do this more efficiently?)
    NetworkDestroy(sto->combined_BDTDNE);
    NetworkDestroy(sto->combined_nonBDTDNE);
    sto->combined_BDTDNE = NetworkInitialize(NULL, NULL, 0, N_NODES, DIRECTED, BIPARTITE, FALSE, 0, NULL);
    sto->combined_nonBDTDNE = NetworkInitialize(NULL, NULL, 0, N_NODES, DIRECTED, BIPARTITE, FALSE, 0, NULL);
  }
}

MH_P_FN(MH_discordBDTNT) {    
  GET_STORAGE(discordBDTNTStorage, sto);
  GET_AUX_STORAGE(StoreTimeAndLasttoggle, dur_inf);

  int nedges = EDGECOUNT(nwp);
  
  int in_network;
  int in_discord;
  
  int nddyads = sto->discordantEdges->nedges + sto->BDTDNE->nedges;
  
  if(nddyads == 0 || unif_rand() < 0.5) {
    // propose from network
    
    // if currentdyads == 0, we *must* propose toggling off an existing edge;
    // the case nedges == 0 && currentdyads == 0 was excluded during initialization,
    // and we cannot end up in that case if we don't start in that case
    // (assuming the initial network is valid)  
    if((unif_rand() < 0.5 && nedges > 0) || (sto->currentdyads == 0)) {
      // select an existing edge at random, and propose toggling it off
      if(unif_rand() < ((double) sto->nonDiscordantEdges->nedges)/nedges) {
        UnsrtELGetRand(Mtail, Mhead, sto->nonDiscordantEdges);
        in_discord = FALSE;
      } else {
        UnsrtELGetRand(Mtail, Mhead, sto->discordantEdges);
        in_discord = TRUE;
      }
          
      in_network = TRUE;
    } else {
      // select a BD-toggleable dyad and propose toggling it
      GetRandDyadFromLists(Mtail, // tail
                           Mhead, // head
                           sto->nodesvec, // tails
                           sto->nodesvec, // heads
                           sto->tailtypes, // tailattrs
                           sto->headtypes, // headattrs
                           sto->attrcounts, // tailcounts
                           sto->attrcounts, // headcounts
                           sto->nmixtypes, // length
                           sto->currentdyads, // dyadcount
                           TRUE, // diagonal; no higher level types here, so always TRUE
                           DIRECTED); // directed; always FALSE in discordBDTNT
          
      in_network = IS_OUTEDGE(Mtail[0],Mhead[0]);
      in_discord = kh_get(DyadMapInt, dur_inf->discord, THKey(dur_inf->discord, Mtail[0], Mhead[0])) != kh_none;

      if(in_network) {
        if(unif_rand() < ((double) sto->nonDiscordantEdges->nedges)/nedges) {
          UnsrtELGetRand(Mtail, Mhead, sto->nonDiscordantEdges);
          in_discord = FALSE;
        } else {
          UnsrtELGetRand(Mtail, Mhead, sto->discordantEdges);
          in_discord = TRUE;
        }
      } else if(in_discord) {
        UnsrtELGetRand(Mtail, Mhead, sto->BDTDNE);
      }
    }
  } else {
    // propose from discord      
    if(unif_rand() < ((double) sto->discordantEdges->nedges)/nddyads) {
      UnsrtELGetRand(Mtail, Mhead, sto->discordantEdges);
      in_network = TRUE;
    } else {
      UnsrtELGetRand(Mtail, Mhead, sto->BDTDNE);
      in_network = FALSE;      
    }
    in_discord = TRUE;
  }
  
  sto->in_discord = in_discord;
  
  sto->tailtype = sto->vattr[Mtail[0]];
  sto->headtype = sto->vattr[Mhead[0]];    
  
  sto->tailmaxl = IN_DEG[Mtail[0]] + OUT_DEG[Mtail[0]] == sto->bound - 1 + in_network;
  sto->headmaxl = IN_DEG[Mhead[0]] + OUT_DEG[Mhead[0]] == sto->bound - 1 + in_network;   
        
  // the count of dyads that can be toggled in the "GetRandBDDyad" branch,
  // in the proposed network
  sto->proposeddyads = sto->currentdyads;
  
  int delta = in_network ? +1 : -1;
  
  for(int i = 0; i < sto->nmixtypes; i++) {
    int corr = 0;
    int ha = 0;

    if(sto->tailtype == sto->headtypes[i] && sto->tailmaxl) {
      ha += delta;
      corr += sto->attrcounts[sto->tailtypes[i]];
    }
      
    if(sto->headtype == sto->headtypes[i] && sto->headmaxl) {
      ha += delta;
      corr += sto->attrcounts[sto->tailtypes[i]];        
    }
      
    if(sto->tailtype == sto->tailtypes[i] && sto->tailmaxl) {
      corr += sto->attrcounts[sto->headtypes[i]] + ha;
    }
      
    if(sto->headtype == sto->tailtypes[i] && sto->headmaxl) {
      corr += sto->attrcounts[sto->headtypes[i]] + ha;
    }
      
    if(sto->tailtypes[i] == sto->headtypes[i]) {
      if(in_network) {
        corr -= ha;
      } else {
        corr += ha;
      }
      corr /= 2;
    }
    
    if(in_network) {
      sto->proposeddyads += corr;
    } else {
      sto->proposeddyads -= corr;
    }
  }

  sto->proposedsubmaxledges = sto->currentsubmaxledges;
  
  // if we are adding an edge that will be submaximal in the post-toggle 
  // network, then increment proposedsubmaxledges for this particular edge
  if(!in_network && !sto->tailmaxl && !sto->headmaxl) {
    sto->proposedsubmaxledges++;
  }
  
  // if we are removing an edge that is submaximal in the current
  // network, decrement proposedsubmaxledges for this particular edge
  if(in_network && !sto->tailmaxl && !sto->headmaxl) {
    sto->proposedsubmaxledges--;
  }

  Edge e;
  Vertex v;

  // if tail will change maximality on toggle, then adjust
  // proposedsubmaxledges for all edges between tail and
  // a submaximal neighbor v, taking care not to count head,
  // since that was handled separately above
  if(sto->tailmaxl) {
    STEP_THROUGH_OUTEDGES(Mtail[0], e, v) {
      if(v != Mhead[0] && IN_DEG[v] + OUT_DEG[v] < sto->bound) {
        sto->proposedsubmaxledges += delta;
      }
    }
    STEP_THROUGH_INEDGES(Mtail[0], e, v) {
      if(IN_DEG[v] + OUT_DEG[v] < sto->bound) {
        sto->proposedsubmaxledges += delta;
      }
    }
  }
  
  // ditto head
  if(sto->headmaxl) {
    STEP_THROUGH_OUTEDGES(Mhead[0], e, v) {
      if(IN_DEG[v] + OUT_DEG[v] < sto->bound) {
        sto->proposedsubmaxledges += delta;
      }
    }
    STEP_THROUGH_INEDGES(Mhead[0], e, v) {
      if(v != Mtail[0] && IN_DEG[v] + OUT_DEG[v] < sto->bound) {
        sto->proposedsubmaxledges += delta;
      }
    }
  }


  // how nddyads can change:
  // the dyad we toggle can add/remove one discordant dyad
  // discordant nonedges can change toggleability status based on the toggle we make,
  // but only if they have Mtail or Mhead as an endpoint
  // so check BDTDNE and/or nonBDTDNE neighbors of Mtail and Mhead, other than Mtail -> Mhead itself
  // which should be taken into account at the outset
  int propnddyads = nddyads;
  
  // direct effect of the toggle we are proposing
  if(in_discord) {
    propnddyads--;
  } else {
    propnddyads++;
  }
  
  // indirect effect of causing other discordant nonedges to change toggleability status
  if(!in_network) {
    // may reduce propnddyads; subtract (i.e. add) in_discord to avoid repeatedly counting the Mtail[0] -> Mhead[0] edge
    if(sto->tailmaxl) {
      propnddyads -= sto->combined_BDTDNE->indegree[Mtail[0]] + sto->combined_BDTDNE->outdegree[Mtail[0]] - in_discord;
    }
    
    if(sto->headmaxl) {
      propnddyads -= sto->combined_BDTDNE->indegree[Mhead[0]] + sto->combined_BDTDNE->outdegree[Mhead[0]] - in_discord;
    }
  } else {
    // may increase propnddyads, but only if other endpoint is also submaximal
    if(sto->tailmaxl) {
      STEP_THROUGH_OUTEDGES_NET(Mtail[0], e, v, sto->combined_nonBDTDNE) {
        if(IN_DEG[v] + OUT_DEG[v] < sto->bound) {
          propnddyads++;
        }
      }
      STEP_THROUGH_INEDGES_NET(Mtail[0], e, v, sto->combined_nonBDTDNE) {
        if(IN_DEG[v] + OUT_DEG[v] < sto->bound) {
          propnddyads++;
        }
      }
    }      
    
    if(sto->headmaxl) {
      STEP_THROUGH_OUTEDGES_NET(Mhead[0], e, v, sto->combined_nonBDTDNE) {
        if(IN_DEG[v] + OUT_DEG[v] < sto->bound) {
          propnddyads++;
        }
      }
      STEP_THROUGH_INEDGES_NET(Mhead[0], e, v, sto->combined_nonBDTDNE) {
        if(IN_DEG[v] + OUT_DEG[v] < sto->bound) {
          propnddyads++;
        }
      }
    }
  }
  
  double forward_network = in_network ? (sto->currentdyads == 0 ? 1.0/nedges : 0.5/nedges + (0.5/sto->currentdyads)*((double)sto->currentsubmaxledges/nedges)) : (nedges == 0 ? 1.0/sto->currentdyads : 0.5/sto->currentdyads);
  
  double forward_discord = in_discord ? 1.0/nddyads : 0;
  
  double backward_network = in_network ? (nedges == 1 ? 1.0/sto->proposeddyads : 0.5/sto->proposeddyads) : (sto->proposeddyads == 0 ? 1.0/(nedges + 1) : 0.5/(nedges + 1) + (0.5/sto->proposeddyads)*((double) sto->proposedsubmaxledges/(nedges + 1)));
  
  double backward_discord = in_discord ? 0 : 1.0/propnddyads;
  
  if(nddyads == 0) forward_network *= 2;
  if(propnddyads == 0) backward_network *= 2;
  
  double forward = forward_discord + forward_network;
  double backward = backward_discord + backward_network;

  MHp->logratio = log(backward/forward);
}

// this U_FN is called *before* the toggle is made in the network
MH_U_FN(Mu_discordBDTNT) {  
  GET_STORAGE(discordBDTNTStorage, sto);
  
  // add or remove the dyad being toggled from the relevant edge set
  if(edgeflag) {
    if(sto->in_discord) {
      // discordant edge
      UnsrtELDelete(tail, head, sto->discordantEdges);
    } else {
      // nondiscordant edge; will become BDTDNE
      UnsrtELDelete(tail, head, sto->nonDiscordantEdges);
      UnsrtELInsert(tail, head, sto->BDTDNE);
      ToggleKnownEdge(tail, head, sto->combined_BDTDNE, FALSE);
    }
  } else {
    if(sto->in_discord) {
      // discordant non-edge; evidently BD toggleable
      UnsrtELDelete(tail, head, sto->BDTDNE);
      ToggleKnownEdge(tail, head, sto->combined_BDTDNE, TRUE);
      UnsrtELInsert(tail, head, sto->nonDiscordantEdges);
    } else {
      // nondiscordant nonedge; will become discordantEdge
      UnsrtELInsert(tail, head, sto->discordantEdges);        
    }
  }
  
  // update current dyad count
  sto->currentdyads = sto->proposeddyads;  
  
  // update the current submaximal edge count
  sto->currentsubmaxledges = sto->proposedsubmaxledges;  
  
  Edge e;
  Vertex v;  
  
  // if (sub)maximality/BD toggleability has changed for any nodes/dyads, update storage accordingly
  if(edgeflag) {
    if(sto->tailmaxl) {
      // tail will be newly submaxl after toggle, so add it to the appropriate node list
      sto->nodesvec[sto->tailtype][sto->attrcounts[sto->tailtype]] = tail;
      sto->nodepos[tail] = sto->attrcounts[sto->tailtype];
      sto->attrcounts[sto->tailtype]++;
      
      // iterate over all nonBDTDNEs with tail as an endpoint; if other endpoint is also
      // submaxl, move this edge from nonBDTDNE to BDTDNE
      
      sto->transferEL->nedges = 0;      
      
      STEP_THROUGH_OUTEDGES_NET(tail, e, v, sto->combined_nonBDTDNE) {
        if(IN_DEG[v] + OUT_DEG[v] < sto->bound) {
          UnsrtELDelete(tail, v, sto->nonBDTDNE);
          UnsrtELInsert(tail, v, sto->BDTDNE);
          UnsrtELInsert(tail, v, sto->transferEL);
        }
      }
      
      STEP_THROUGH_INEDGES_NET(tail, e, v, sto->combined_nonBDTDNE) {
        if(IN_DEG[v] + OUT_DEG[v] < sto->bound) {
          UnsrtELDelete(v, tail, sto->nonBDTDNE);
          UnsrtELInsert(v, tail, sto->BDTDNE);
          UnsrtELInsert(v, tail, sto->transferEL);
        }
      }

      for(int i = 1; i <= sto->transferEL->nedges; i++) {
        ToggleKnownEdge(sto->transferEL->tails[i], sto->transferEL->heads[i], sto->combined_nonBDTDNE, TRUE);
        ToggleKnownEdge(sto->transferEL->tails[i], sto->transferEL->heads[i], sto->combined_BDTDNE, FALSE);        
      }
    }
    
    if(sto->headmaxl) {
      // head will be newly submaxl after toggle, so add it to the appropriate node list    
      sto->nodesvec[sto->headtype][sto->attrcounts[sto->headtype]] = head;
      sto->nodepos[head] = sto->attrcounts[sto->headtype];
      sto->attrcounts[sto->headtype]++;

      // iterate over all nonBDTDNEs with head as an endpoint; if other endpoint is also
      // submaxl, move this edge from nonBDTDNE to BDTDNE
      sto->transferEL->nedges = 0;      
      
      STEP_THROUGH_OUTEDGES_NET(head, e, v, sto->combined_nonBDTDNE) {
        if(IN_DEG[v] + OUT_DEG[v] < sto->bound) {
          UnsrtELDelete(head, v, sto->nonBDTDNE);
          UnsrtELInsert(head, v, sto->BDTDNE);
          UnsrtELInsert(head, v, sto->transferEL);
        }
      }
      
      STEP_THROUGH_INEDGES_NET(head, e, v, sto->combined_nonBDTDNE) {
        if(IN_DEG[v] + OUT_DEG[v] < sto->bound) {
          UnsrtELDelete(v, head, sto->nonBDTDNE);
          UnsrtELInsert(v, head, sto->BDTDNE);
          UnsrtELInsert(v, head, sto->transferEL);
        }
      }

      for(int i = 1; i <= sto->transferEL->nedges; i++) {
        ToggleKnownEdge(sto->transferEL->tails[i], sto->transferEL->heads[i], sto->combined_nonBDTDNE, TRUE);
        ToggleKnownEdge(sto->transferEL->tails[i], sto->transferEL->heads[i], sto->combined_BDTDNE, FALSE);        
      }
    }
  } else {
    if(sto->tailmaxl) {
      // tail will be newly maxl after toggle, so remove it from the appropriate node list, updating nodepos
      // for whatever node is currently at the end of sto->nodesvec[sto->tailtype], since that will be
      // moved into tail's position
      sto->nodesvec[sto->tailtype][sto->nodepos[tail]] = sto->nodesvec[sto->tailtype][sto->attrcounts[sto->tailtype] - 1];
      sto->nodepos[sto->nodesvec[sto->tailtype][sto->nodepos[tail]]] = sto->nodepos[tail];
      sto->attrcounts[sto->tailtype]--;

      // transfer all BDTDNEs with tail as an endpoint to nonBDTDNE
      sto->transferEL->nedges = 0;

      STEP_THROUGH_OUTEDGES_NET(tail, e, v, sto->combined_BDTDNE) {
        UnsrtELDelete(tail, v, sto->BDTDNE);
        UnsrtELInsert(tail, v, sto->nonBDTDNE);
        UnsrtELInsert(tail, v, sto->transferEL);        
      }

      STEP_THROUGH_INEDGES_NET(tail, e, v, sto->combined_BDTDNE) {
        UnsrtELDelete(v, tail, sto->BDTDNE);
        UnsrtELInsert(v, tail, sto->nonBDTDNE);
        UnsrtELInsert(v, tail, sto->transferEL);
      }

      for(int i = 1; i <= sto->transferEL->nedges; i++) {
        ToggleKnownEdge(sto->transferEL->tails[i], sto->transferEL->heads[i], sto->combined_BDTDNE, TRUE);
        ToggleKnownEdge(sto->transferEL->tails[i], sto->transferEL->heads[i], sto->combined_nonBDTDNE, FALSE);        
      }
    }
    
    if(sto->headmaxl) {
      // head will be newly maxl after toggle, so remove it from the appropriate node list, updating nodepos
      // for whatever node is currently at the end of sto->nodesvec[sto->headtype], since that will be
      // moved into head's position
      sto->nodesvec[sto->headtype][sto->nodepos[head]] = sto->nodesvec[sto->headtype][sto->attrcounts[sto->headtype] - 1];
      sto->nodepos[sto->nodesvec[sto->headtype][sto->nodepos[head]]] = sto->nodepos[head];
      sto->attrcounts[sto->headtype]--;

      sto->transferEL->nedges = 0;

      STEP_THROUGH_OUTEDGES_NET(head, e, v, sto->combined_BDTDNE) {
        UnsrtELDelete(head, v, sto->BDTDNE);
        UnsrtELInsert(head, v, sto->nonBDTDNE);
        UnsrtELInsert(head, v, sto->transferEL);        
      }

      STEP_THROUGH_INEDGES_NET(head, e, v, sto->combined_BDTDNE) {
        UnsrtELDelete(v, head, sto->BDTDNE);
        UnsrtELInsert(v, head, sto->nonBDTDNE);
        UnsrtELInsert(v, head, sto->transferEL);
      }

      for(int i = 1; i <= sto->transferEL->nedges; i++) {
        ToggleKnownEdge(sto->transferEL->tails[i], sto->transferEL->heads[i], sto->combined_BDTDNE, TRUE);
        ToggleKnownEdge(sto->transferEL->tails[i], sto->transferEL->heads[i], sto->combined_nonBDTDNE, FALSE);        
      }
    }      
  }
}

MH_F_FN(Mf_discordBDTNT) {
  GET_STORAGE(discordBDTNTStorage, sto);

  // let BDTNT's F_FN do most of the work
  MH_STORAGE = sto->static_sto;
  Mf_BDTNT(MHp, nwp);
  MH_STORAGE = sto;
  Free(sto->static_sto);
  
  UnsrtELDestroy(sto->BDTDNE);
  UnsrtELDestroy(sto->nonBDTDNE);
  UnsrtELDestroy(sto->discordantEdges);  
  NetworkDestroy(sto->combined_BDTDNE);
  NetworkDestroy(sto->combined_nonBDTDNE);
  UnsrtELDestroy(sto->transferEL);
  // MHp->storage itself should be Freed by MHProposalDestroy
}



/********************
    MH_discordBDStratTNT
********************/

typedef struct {
  Vertex ***nodesvec;
  int **attrcounts;
  
  int *nodepos;

  int strattailtype;
  int bdtailtype;
  int tailmaxl;
  
  int stratheadtype;
  int bdheadtype;  
  int headmaxl;
  
  int stratmixingtype;
  
  double currentcumprob;
  double proposedcumprob;
  
  double *originalprobvec;
  
  WtPop *wtp;

  int bound;
  int nmixtypes;
  
  int *strat_vattr;
  int *bd_vattr;
  
  int *BDtypesbyStrattype;
  int **BDtailsbyStrattype;
  int **BDheadsbyStrattype;
  
  int *strattailtypes;
  int *stratheadtypes;

  Network *combined_BDTDNE;
  Network *combined_nonBDTDNE;

  UnsrtEL *transferEL;

  UnsrtEL **BDTDNE;
  UnsrtEL **nonBDTDNE;
  
  UnsrtEL **discordantEdges;
  UnsrtEL **nonDiscordantEdges;
  
  int in_discord;
  
  int nstratlevels;
  
  int *currentsubmaxledgestype;
  int **indmat;  
  
  int nmixtypestoupdate;
  int *mixtypestoupdate;
  
  BDStratTNTStorage *static_sto;
} discordBDStratTNTStorage;

MH_I_FN(Mi_discordBDStratTNT) {
  // let BDStratTNT's I_FN do most of the work
  Mi_BDStratTNT(MHp, nwp);
  
  // clean up from static proposal initialization and handle purely temporal aspects
  BDStratTNTStorage *static_sto = MH_STORAGE;
  ALLOC_STORAGE(1, discordBDStratTNTStorage, sto);
  sto->nodesvec = static_sto->nodesvec;
  sto->attrcounts = static_sto->attrcounts;
  sto->nodepos = static_sto->nodepos;
  sto->currentcumprob = static_sto->currentcumprob;
  sto->originalprobvec = static_sto->originalprobvec;
  sto->wtp = static_sto->wtp;
  sto->bound = static_sto->bound;
  sto->nmixtypes = static_sto->nmixtypes;
  sto->strat_vattr = static_sto->strat_vattr;
  sto->bd_vattr = static_sto->bd_vattr;
  sto->BDtypesbyStrattype = static_sto->BDtypesbyStrattype;
  sto->BDtailsbyStrattype = static_sto->BDtailsbyStrattype;
  sto->BDheadsbyStrattype = static_sto->BDheadsbyStrattype;
  sto->strattailtypes = static_sto->strattailtypes;
  sto->stratheadtypes = static_sto->stratheadtypes;
  sto->nonDiscordantEdges = static_sto->els;
  sto->nstratlevels = static_sto->nstratlevels;
  sto->currentsubmaxledgestype = static_sto->currentsubmaxledgestype;
  sto->indmat = static_sto->indmat;
  sto->mixtypestoupdate = static_sto->mixtypestoupdate;
  sto->static_sto = static_sto;
  
  sto->BDTDNE = Calloc(sto->nmixtypes, UnsrtEL *);
  sto->nonBDTDNE = Calloc(sto->nmixtypes, UnsrtEL *);
  sto->discordantEdges = Calloc(sto->nmixtypes, UnsrtEL *);  
  for(int i = 0; i < sto->nmixtypes; i++) {
    sto->BDTDNE[i] = UnsrtELInitialize(0, NULL, NULL, FALSE);
    sto->nonBDTDNE[i] = UnsrtELInitialize(0, NULL, NULL, FALSE);
    sto->discordantEdges[i] = UnsrtELInitialize(0, NULL, NULL, FALSE);
  }
  sto->combined_BDTDNE = NetworkInitialize(NULL, NULL, 0, N_NODES, DIRECTED, BIPARTITE, FALSE, 0, NULL);
  sto->combined_nonBDTDNE = NetworkInitialize(NULL, NULL, 0, N_NODES, DIRECTED, BIPARTITE, FALSE, 0, NULL);
  sto->transferEL = UnsrtELInitialize(0, NULL, NULL, FALSE);
}

MH_X_FN(Mx_discordBDStratTNT) {
  if(type == TICK) {
    GET_STORAGE(discordBDStratTNTStorage, sto);

    for(int i = 0; i < sto->nmixtypes; i++) {
      // transfer discordantEdges to nonDiscordantEdges
      for(int j = 1; j <= sto->discordantEdges[i]->nedges; j++) {
        UnsrtELInsert(sto->discordantEdges[i]->tails[j], sto->discordantEdges[i]->heads[j], sto->nonDiscordantEdges[i]);
      }

      // clear all the discordance information
      sto->BDTDNE[i]->nedges = 0;
      sto->nonBDTDNE[i]->nedges = 0;
      sto->discordantEdges[i]->nedges = 0;      
    }
    
    // for now, destroy and recreate each time step (can we do this more efficiently?)
    NetworkDestroy(sto->combined_BDTDNE);
    NetworkDestroy(sto->combined_nonBDTDNE);
    sto->combined_BDTDNE = NetworkInitialize(NULL, NULL, 0, N_NODES, DIRECTED, BIPARTITE, FALSE, 0, NULL);
    sto->combined_nonBDTDNE = NetworkInitialize(NULL, NULL, 0, N_NODES, DIRECTED, BIPARTITE, FALSE, 0, NULL);
  }
}

MH_P_FN(MH_discordBDStratTNT) {
  GET_STORAGE(discordBDStratTNTStorage, sto);
  GET_AUX_STORAGE(StoreTimeAndLasttoggle, dur_inf);

  // sample a toggleable strat mixing type on which to make a proposal
  int strat_i = WtPopGetRand(sto->wtp);
  
  // record the mixing type of the toggle, in case it's needed in the U function later
  sto->stratmixingtype = strat_i;    

  int strattailtype = sto->strattailtypes[strat_i];
  int stratheadtype = sto->stratheadtypes[strat_i];
    
  // number of edges of this mixing type
  int nedgestype = sto->nonDiscordantEdges[strat_i]->nedges + sto->discordantEdges[strat_i]->nedges;
  
  Dyad ndyadstype = 0;
  for(int j = 0; j < sto->BDtypesbyStrattype[strat_i]; j++) {
    if(strattailtype == stratheadtype && sto->BDtailsbyStrattype[strat_i][j] == sto->BDheadsbyStrattype[strat_i][j]) {
      ndyadstype += (Dyad)sto->attrcounts[strattailtype][sto->BDtailsbyStrattype[strat_i][j]]*(sto->attrcounts[stratheadtype][sto->BDheadsbyStrattype[strat_i][j]] - 1)/2;
    } else {
      ndyadstype += (Dyad)sto->attrcounts[strattailtype][sto->BDtailsbyStrattype[strat_i][j]]*sto->attrcounts[stratheadtype][sto->BDheadsbyStrattype[strat_i][j]];
    }
  }
  
  int nddyadstype = sto->discordantEdges[strat_i]->nedges + sto->BDTDNE[strat_i]->nedges;
  
  int in_network;
  int in_discord;
  
  if(nddyadstype == 0 || unif_rand() < 0.5) {
    // propose from network
    if((unif_rand() < 0.5 && nedgestype > 0) || ndyadstype == 0) {
      // propose toggling off an existing edge of strat mixing type strat_i
      if(unif_rand() < ((double) sto->nonDiscordantEdges[strat_i]->nedges)/nedgestype) {
        UnsrtELGetRand(Mtail, Mhead, sto->nonDiscordantEdges[strat_i]);
        in_discord = FALSE;
      } else {
        UnsrtELGetRand(Mtail, Mhead, sto->discordantEdges[strat_i]);
        in_discord = TRUE;
      }
      in_network = TRUE;
    } else {
      // select a random BD toggleable dyad of strat mixing type strat_i and propose toggling it
      GetRandDyadFromLists(Mtail, // tail
                           Mhead, // head
                           sto->nodesvec[strattailtype], // tails
                           sto->nodesvec[stratheadtype], // heads
                           sto->BDtailsbyStrattype[strat_i], // tailattrs
                           sto->BDheadsbyStrattype[strat_i], // headattrs
                           sto->attrcounts[strattailtype], // tailcounts
                           sto->attrcounts[stratheadtype], // headcounts
                           sto->BDtypesbyStrattype[strat_i], // length
                           ndyadstype, // dyadcount
                           strattailtype == stratheadtype, // diagonal
                           DIRECTED); // directed; always FALSE in discordBDStratTNT
         
      in_network = IS_OUTEDGE(Mtail[0],Mhead[0]);
      in_discord = kh_get(DyadMapInt, dur_inf->discord, THKey(dur_inf->discord, Mtail[0], Mhead[0])) != kh_none;

      if(in_network) {
        if(unif_rand() < ((double) sto->nonDiscordantEdges[strat_i]->nedges/nedgestype)) {
          UnsrtELGetRand(Mtail, Mhead, sto->nonDiscordantEdges[strat_i]);          
          in_discord = FALSE;          
        } else {
          UnsrtELGetRand(Mtail, Mhead, sto->discordantEdges[strat_i]);
          in_discord = TRUE;
        }
      } else if(in_discord) {
        UnsrtELGetRand(Mtail, Mhead, sto->BDTDNE[strat_i]);          
      }
    }
  } else {
    // propose from discord
    if(unif_rand() < ((double) sto->discordantEdges[strat_i]->nedges)/nddyadstype) {
      UnsrtELGetRand(Mtail, Mhead, sto->discordantEdges[strat_i]);
      in_network = TRUE;
    } else {
      UnsrtELGetRand(Mtail, Mhead, sto->BDTDNE[strat_i]);
      in_network = FALSE;
    }
    in_discord = TRUE;
  }

  sto->in_discord = in_discord;

  sto->strattailtype = sto->strat_vattr[Mtail[0]];
  sto->stratheadtype = sto->strat_vattr[Mhead[0]];
    
  sto->bdtailtype = sto->bd_vattr[Mtail[0]];
  sto->bdheadtype = sto->bd_vattr[Mhead[0]];

  sto->tailmaxl = IN_DEG[Mtail[0]] + OUT_DEG[Mtail[0]] == sto->bound - 1 + in_network;
  sto->headmaxl = IN_DEG[Mhead[0]] + OUT_DEG[Mhead[0]] == sto->bound - 1 + in_network;
  
  // here we compute the proposedcumprob, checking only those
  // mixing types that can be influenced by toggles made on 
  // the current mixing type
  sto->proposedcumprob = sto->currentcumprob;
  sto->nmixtypestoupdate = 0; // reset counter  
  // avoid these somewhat expensive checks in the typical case
  // where you have enough submaximal nodes that you cannot
  // be exhausting any mixing types of toggleable dyads
  if(sto->attrcounts[sto->strattailtype][sto->bdtailtype] <= 2 || sto->attrcounts[sto->stratheadtype][sto->bdheadtype] <= 2) {  
    if(sto->proposedcumprob != 1 && in_network) {
      // we are proposing removing an edge of the current mixing type and need to make sure 
      // that any influenced mixing type that couldn't be toggled before this toggle
      // but could be toggled after this toggle gets its prob added to proposedcumprob;
      // note that if proposedcumprob = 1 in the test above then all mixing types were
      // toggleable before this off-toggle and that will remain the case afterward, so
      // there is nothing to do
      int ntocheck = (2 - ((sto->strattailtype == sto->stratheadtype) || BIPARTITE))*sto->nstratlevels;
      for(int i = 0; i < ntocheck; i++) {
        int infl_i;
        if(i < sto->nstratlevels) {
          infl_i = sto->indmat[sto->strattailtype][i];
        } else {
          infl_i = sto->indmat[i - sto->nstratlevels][sto->stratheadtype];  
        }
        if(infl_i < 0 || infl_i == sto->stratmixingtype) {
          continue;
        }

        if(WtPopGetWt(infl_i, sto->wtp) > 0) {
          continue;
        }
        // else there are no toggleable dyads of type infl_i in the current network;
        // we need to check if that changes after hypothetically removing the proposed edge
        int anytoggleable = FALSE;
        
        for(int j = 0; j < sto->BDtypesbyStrattype[infl_i]; j++) {
          // adjustments
          int proposedtailadjustment = (sto->strattailtype == sto->strattailtypes[infl_i] && sto->bdtailtype == sto->BDtailsbyStrattype[infl_i][j] && sto->tailmaxl) + (sto->stratheadtype == sto->strattailtypes[infl_i] && sto->bdheadtype == sto->BDtailsbyStrattype[infl_i][j] && sto->headmaxl);
          int proposedheadadjustment = (sto->strattailtype == sto->stratheadtypes[infl_i] && sto->bdtailtype == sto->BDheadsbyStrattype[infl_i][j] && sto->tailmaxl) + (sto->stratheadtype == sto->stratheadtypes[infl_i] && sto->bdheadtype == sto->BDheadsbyStrattype[infl_i][j] && sto->headmaxl);
          
          int tailcounts = sto->attrcounts[sto->strattailtypes[infl_i]][sto->BDtailsbyStrattype[infl_i][j]];
          int headcounts = sto->attrcounts[sto->stratheadtypes[infl_i]][sto->BDheadsbyStrattype[infl_i][j]];
          
          proposedtailadjustment = -proposedtailadjustment;
          proposedheadadjustment = -proposedheadadjustment;
          
          if(tailcounts > proposedtailadjustment && headcounts > proposedheadadjustment + (sto->strattailtypes[infl_i] == sto->stratheadtypes[infl_i] && sto->BDtailsbyStrattype[infl_i][j] == sto->BDheadsbyStrattype[infl_i][j])) {
            anytoggleable = TRUE;
            break;
          }      
        }
        
        if(anytoggleable) {
          sto->proposedcumprob += sto->originalprobvec[infl_i];
          sto->mixtypestoupdate[sto->nmixtypestoupdate] = infl_i;
          sto->nmixtypestoupdate++;          
        }
      }
    } else if(!in_network) {
      // we are proposing toggling on an edge of the current mixing type, and need to make
      // sure that any influenced mixing type that could be toggled before this toggle
      // can also be toggled after this toggle, or else has its prob subtracted accordingly
      int ntocheck = (2 - ((sto->strattailtype == sto->stratheadtype) || BIPARTITE))*sto->nstratlevels;
      for(int i = 0; i < ntocheck; i++) {
        int infl_i;
        if(i < sto->nstratlevels) {
          infl_i = sto->indmat[sto->strattailtype][i];
        } else {
          infl_i = sto->indmat[i - sto->nstratlevels][sto->stratheadtype];  
        }
        if(infl_i < 0 || infl_i == sto->stratmixingtype) {
          continue;
        }

        if(WtPopGetWt(infl_i, sto->wtp) == 0 || sto->nonDiscordantEdges[infl_i]->nedges > 0 || sto->discordantEdges[infl_i]->nedges > 0) {
          continue;
        }
        // else there are no toggleable dyads of type infl_i in the current network;
        // we need to check if that changes after hypothetically removing the proposed edge
        int anytoggleable = FALSE;
        
        for(int j = 0; j < sto->BDtypesbyStrattype[infl_i]; j++) {
          // adjustments
          int proposedtailadjustment = (sto->strattailtype == sto->strattailtypes[infl_i] && sto->bdtailtype == sto->BDtailsbyStrattype[infl_i][j] && sto->tailmaxl) + (sto->stratheadtype == sto->strattailtypes[infl_i] && sto->bdheadtype == sto->BDtailsbyStrattype[infl_i][j] && sto->headmaxl);
          int proposedheadadjustment = (sto->strattailtype == sto->stratheadtypes[infl_i] && sto->bdtailtype == sto->BDheadsbyStrattype[infl_i][j] && sto->tailmaxl) + (sto->stratheadtype == sto->stratheadtypes[infl_i] && sto->bdheadtype == sto->BDheadsbyStrattype[infl_i][j] && sto->headmaxl);
          
          int tailcounts = sto->attrcounts[sto->strattailtypes[infl_i]][sto->BDtailsbyStrattype[infl_i][j]];
          int headcounts = sto->attrcounts[sto->stratheadtypes[infl_i]][sto->BDheadsbyStrattype[infl_i][j]];
                
          if(tailcounts > proposedtailadjustment && headcounts > proposedheadadjustment + (sto->strattailtypes[infl_i] == sto->stratheadtypes[infl_i] && sto->BDtailsbyStrattype[infl_i][j] == sto->BDheadsbyStrattype[infl_i][j])) {
            anytoggleable = TRUE;
            break;
          }      
        }
        
        if(!anytoggleable) {
          sto->proposedcumprob -= sto->originalprobvec[infl_i];
          sto->mixtypestoupdate[sto->nmixtypestoupdate] = infl_i;
          sto->nmixtypestoupdate++;          
        }
      }    
    }
  }
  
  // need to compute proposed dyad count for current mixing type (only)
  Dyad proposeddyadstype = ndyadstype;

  int delta = in_network ? +1 : -1;

  for(int j = 0; j < sto->BDtypesbyStrattype[strat_i]; j++) {
    int corr = 0;
    int ha = 0;

    if(sto->strattailtype == sto->stratheadtypes[strat_i] && sto->bdtailtype == sto->BDheadsbyStrattype[strat_i][j] && sto->tailmaxl) {
      ha += delta;
      corr += sto->attrcounts[sto->strattailtypes[strat_i]][sto->BDtailsbyStrattype[strat_i][j]];
    }
      
    if(sto->stratheadtype == sto->stratheadtypes[strat_i] && sto->bdheadtype == sto->BDheadsbyStrattype[strat_i][j] && sto->headmaxl) {
      ha += delta;
      corr += sto->attrcounts[sto->strattailtypes[strat_i]][sto->BDtailsbyStrattype[strat_i][j]];        
    }
      
    if(sto->strattailtype == sto->strattailtypes[strat_i] && sto->bdtailtype == sto->BDtailsbyStrattype[strat_i][j] && sto->tailmaxl) {
      corr += sto->attrcounts[sto->stratheadtypes[strat_i]][sto->BDheadsbyStrattype[strat_i][j]] + ha;
    }
      
    if(sto->stratheadtype == sto->strattailtypes[strat_i] && sto->bdheadtype == sto->BDtailsbyStrattype[strat_i][j] && sto->headmaxl) {
      corr += sto->attrcounts[sto->stratheadtypes[strat_i]][sto->BDheadsbyStrattype[strat_i][j]] + ha;
    }
      
    if(sto->strattailtypes[strat_i] == sto->stratheadtypes[strat_i] && sto->BDtailsbyStrattype[strat_i][j] == sto->BDheadsbyStrattype[strat_i][j]) {
      if(in_network) {
        corr -= ha;
      } else {
        corr += ha;
      }
      corr /= 2;
    }
    
    if(in_network) {
      proposeddyadstype += corr;
    } else {
      proposeddyadstype -= corr;
    }      
  }
    
  
  Edge e;
  Vertex v;

  // need propnddyads count
  
  int propnddyadstype = nddyadstype;
  
  // direct effect of the toggle we are proposing
  if(in_discord) {
    propnddyadstype--;
  } else {
    propnddyadstype++;
  }
  
  // indirect effect of causing other discordant nonedges to change toggleability status
  if(!in_network) {
    // may reduce propnddyads; add in_discord to avoid repeatedly counting the Mtail[0] -> Mhead[0] edge
    if(sto->tailmaxl) {
      propnddyadstype += in_discord;
      STEP_THROUGH_OUTEDGES_NET(Mtail[0], e, v, sto->combined_BDTDNE) {
        if(sto->indmat[sto->strattailtype][sto->strat_vattr[v]] == strat_i) {
          propnddyadstype--;
        }
      }
      STEP_THROUGH_INEDGES_NET(Mtail[0], e, v, sto->combined_BDTDNE) {
        if(sto->indmat[sto->strattailtype][sto->strat_vattr[v]] == strat_i) {
          propnddyadstype--;
        }
      }
    }
    
    if(sto->headmaxl) {
      propnddyadstype += in_discord;
      STEP_THROUGH_OUTEDGES_NET(Mhead[0], e, v, sto->combined_BDTDNE) {
        if(sto->indmat[sto->stratheadtype][sto->strat_vattr[v]] == strat_i) {
          propnddyadstype--;
        }
      }
      STEP_THROUGH_INEDGES_NET(Mhead[0], e, v, sto->combined_BDTDNE) {
        if(sto->indmat[sto->stratheadtype][sto->strat_vattr[v]] == strat_i) {
          propnddyadstype--;
        }
      }
    }
  } else {
    // may increase propnddyads, but only if other endpoint is also submaximal
    if(sto->tailmaxl) {
      STEP_THROUGH_OUTEDGES_NET(Mtail[0], e, v, sto->combined_nonBDTDNE) {
        if(sto->indmat[sto->strattailtype][sto->strat_vattr[v]] == strat_i && IN_DEG[v] + OUT_DEG[v] < sto->bound) {
          propnddyadstype++;
        }
      }
      STEP_THROUGH_INEDGES_NET(Mtail[0], e, v, sto->combined_nonBDTDNE) {
        if(sto->indmat[sto->strattailtype][sto->strat_vattr[v]] == strat_i && IN_DEG[v] + OUT_DEG[v] < sto->bound) {
          propnddyadstype++;
        }
      }
    }      
    
    if(sto->headmaxl) {
      STEP_THROUGH_OUTEDGES_NET(Mhead[0], e, v, sto->combined_nonBDTDNE) {
        if(sto->indmat[sto->stratheadtype][sto->strat_vattr[v]] == strat_i && IN_DEG[v] + OUT_DEG[v] < sto->bound) {
          propnddyadstype++;
        }
      }
      STEP_THROUGH_INEDGES_NET(Mhead[0], e, v, sto->combined_nonBDTDNE) {
        if(sto->indmat[sto->stratheadtype][sto->strat_vattr[v]] == strat_i && IN_DEG[v] + OUT_DEG[v] < sto->bound) {
          propnddyadstype++;
        }
      }
    }
  }

  int proposedsubmaxledgestype = sto->currentsubmaxledgestype[strat_i];

  // if we are adding an edge that will be submaximal in the post-toggle 
  // network, then increment proposedsubmaxledgestype for this particular edge
  if(!in_network && !sto->tailmaxl && !sto->headmaxl) {
    proposedsubmaxledgestype++;
  }

  // if we are removing an edge that is submaximal in the current
  // network, decrement proposedsubmaxledgestype for this particular edge
  if(in_network && !sto->tailmaxl && !sto->headmaxl) {
    proposedsubmaxledgestype--;
  }

  // if tail will change maximality on toggle, then adjust
  // proposedsubmaxledgestype for all edges between tail and
  // a submaximal neighbor v with the edge between tail and v
  // having the mixing type strat_i, taking care not to count head,
  // since that was handled separately above
  if(sto->tailmaxl) {
    STEP_THROUGH_OUTEDGES(Mtail[0], e, v) {
      if(v != Mhead[0] && IN_DEG[v] + OUT_DEG[v] < sto->bound && sto->indmat[sto->strattailtype][sto->strat_vattr[v]] == strat_i) {
        proposedsubmaxledgestype += delta;
      }
    }
    STEP_THROUGH_INEDGES(Mtail[0], e, v) {
      if(IN_DEG[v] + OUT_DEG[v] < sto->bound && sto->indmat[sto->strattailtype][sto->strat_vattr[v]] == strat_i) {
        proposedsubmaxledgestype += delta;
      }
    }
  }

  // ditto head
  if(sto->headmaxl) {
    STEP_THROUGH_OUTEDGES(Mhead[0], e, v) {
      if(IN_DEG[v] + OUT_DEG[v] < sto->bound && sto->indmat[sto->stratheadtype][sto->strat_vattr[v]] == strat_i) {
        proposedsubmaxledgestype += delta;
      }
    }
    STEP_THROUGH_INEDGES(Mhead[0], e, v) {
      if(v != Mtail[0] && IN_DEG[v] + OUT_DEG[v] < sto->bound && sto->indmat[sto->stratheadtype][sto->strat_vattr[v]] == strat_i) {
        proposedsubmaxledgestype += delta;
      }
    }
  }
    
  double prob_weight = sto->currentcumprob/sto->proposedcumprob;
  
  double forward_network = in_network ? (ndyadstype == 0 ? 1.0/nedgestype : 0.5/nedgestype + (0.5/ndyadstype)*((double)sto->currentsubmaxledgestype[strat_i]/nedgestype)) : (nedgestype == 0 ? 1.0/ndyadstype : 0.5/ndyadstype);
  
  double forward_discord = in_discord ? 1.0/nddyadstype : 0;
  
  double backward_network = in_network ? (nedgestype == 1 ? 1.0/proposeddyadstype : 0.5/proposeddyadstype) : (proposeddyadstype == 0 ? 1.0/(nedgestype + 1) : 0.5/(nedgestype + 1) + (0.5/proposeddyadstype)*((double) proposedsubmaxledgestype/(nedgestype + 1)));
  
  double backward_discord = in_discord ? 0 : 1.0/propnddyadstype;
  
  if(nddyadstype == 0) forward_network *= 2;
  if(propnddyadstype == 0) backward_network *= 2;
  
  double forward = forward_discord + forward_network;
  double backward = backward_discord + backward_network;

  MHp->logratio = log(prob_weight*backward/forward);
  
}

MH_U_FN(Mu_discordBDStratTNT) {   
  GET_STORAGE(discordBDStratTNTStorage, sto);

  // if any strat mixing types have changed toggleability status, update prob info accordingly
  if(sto->nmixtypestoupdate > 0) {
    sto->currentcumprob = sto->proposedcumprob;

    if(edgeflag) {
      for(int i = 0; i < sto->nmixtypestoupdate; i++) {
        WtPopSetWt(sto->mixtypestoupdate[i], sto->originalprobvec[sto->mixtypestoupdate[i]], sto->wtp);          
      }
    } else {
      for(int i = 0; i < sto->nmixtypestoupdate; i++) {
        WtPopSetWt(sto->mixtypestoupdate[i], 0, sto->wtp);          
      }
    }
  }

  // if we are adding an edge that will be submaximal in the post-toggle 
  // network, then increment currentsubmaxledgestype for this particular edge
  if(!edgeflag && !sto->tailmaxl && !sto->headmaxl) {
    sto->currentsubmaxledgestype[sto->stratmixingtype]++;
  }

  // if we are removing an edge that is submaximal in the current
  // network, decrement currentsubmaxledgestype for this particular edge
  if(edgeflag && !sto->tailmaxl && !sto->headmaxl) {
    sto->currentsubmaxledgestype[sto->stratmixingtype]--;
  }

  int delta = edgeflag ? +1 : -1;

  Edge e;
  Vertex v;

  // if tail will change maximality on toggle, then adjust
  // currentsubmaxledgestype for all edges between tail and
  // a submaximal neighbor v, taking care not to count head,
  // since that was handled separately above
  if(sto->tailmaxl) {
    STEP_THROUGH_OUTEDGES(tail, e, v) {
      if(v != head && IN_DEG[v] + OUT_DEG[v] < sto->bound && sto->indmat[sto->strattailtype][sto->strat_vattr[v]] >= 0) {
        sto->currentsubmaxledgestype[sto->indmat[sto->strattailtype][sto->strat_vattr[v]]] += delta;
      }
    }
    STEP_THROUGH_INEDGES(tail, e, v) {
      if(IN_DEG[v] + OUT_DEG[v] < sto->bound && sto->indmat[sto->strattailtype][sto->strat_vattr[v]] >= 0) {
        sto->currentsubmaxledgestype[sto->indmat[sto->strattailtype][sto->strat_vattr[v]]] += delta;
      }
    }
  }

  // ditto head
  if(sto->headmaxl) {
    STEP_THROUGH_OUTEDGES(head, e, v) {
      if(IN_DEG[v] + OUT_DEG[v] < sto->bound && sto->indmat[sto->stratheadtype][sto->strat_vattr[v]] >= 0) {
        sto->currentsubmaxledgestype[sto->indmat[sto->stratheadtype][sto->strat_vattr[v]]] += delta;
      }
    }
    STEP_THROUGH_INEDGES(head, e, v) {
      if(v != tail && IN_DEG[v] + OUT_DEG[v] < sto->bound && sto->indmat[sto->stratheadtype][sto->strat_vattr[v]] >= 0) {
        sto->currentsubmaxledgestype[sto->indmat[sto->stratheadtype][sto->strat_vattr[v]]] += delta;
      }
    }
  }  
    
  // add or remove the dyad being toggled from the relevant edge set
  if(edgeflag) {
    if(sto->in_discord) {
      // discordant edge
      UnsrtELDelete(tail, head, sto->discordantEdges[sto->stratmixingtype]);
    } else {
      // nondiscordant edge; will become BDTDNE
      UnsrtELDelete(tail, head, sto->nonDiscordantEdges[sto->stratmixingtype]);      
      UnsrtELInsert(tail, head, sto->BDTDNE[sto->stratmixingtype]);
      ToggleKnownEdge(tail, head, sto->combined_BDTDNE, FALSE);
    }
  } else {
    if(sto->in_discord) {
      // discordant non-edge; evidently BD toggleable
      UnsrtELDelete(tail, head, sto->BDTDNE[sto->stratmixingtype]);
      ToggleKnownEdge(tail, head, sto->combined_BDTDNE, TRUE);
      UnsrtELInsert(tail, head, sto->nonDiscordantEdges[sto->stratmixingtype]);      
    } else {
      // nondiscordant nonedge; will become discordantEdge
      UnsrtELInsert(tail, head, sto->discordantEdges[sto->stratmixingtype]);        
    }
  }
  
  // if (sub)maximality/BD toggleability has changed for any nodes/dyads, update storage accordingly
  if(edgeflag) {
    if(sto->tailmaxl) {
      // tail will be newly submaxl after toggle, so add it to the appropriate node list
      sto->nodesvec[sto->strattailtype][sto->bdtailtype][sto->attrcounts[sto->strattailtype][sto->bdtailtype]] = tail;
      sto->nodepos[tail] = sto->attrcounts[sto->strattailtype][sto->bdtailtype];
      sto->attrcounts[sto->strattailtype][sto->bdtailtype]++;
      
      // iterate over all nonBDTDNEs with tail as an endpoint; if other endpoint is also
      // submaxl, move this edge from nonBDTDNE to BDTDNE
      
      sto->transferEL->nedges = 0;      
      
      STEP_THROUGH_OUTEDGES_NET(tail, e, v, sto->combined_nonBDTDNE) {
        if(IN_DEG[v] + OUT_DEG[v] < sto->bound) {
          int stratmixingtype = sto->indmat[sto->strattailtype][sto->strat_vattr[v]];
          UnsrtELDelete(tail, v, sto->nonBDTDNE[stratmixingtype]);
          UnsrtELInsert(tail, v, sto->BDTDNE[stratmixingtype]);
          UnsrtELInsert(tail, v, sto->transferEL);
        }
      }
      
      STEP_THROUGH_INEDGES_NET(tail, e, v, sto->combined_nonBDTDNE) {
        if(IN_DEG[v] + OUT_DEG[v] < sto->bound) {
          int stratmixingtype = sto->indmat[sto->strattailtype][sto->strat_vattr[v]];
          UnsrtELDelete(v, tail, sto->nonBDTDNE[stratmixingtype]);
          UnsrtELInsert(v, tail, sto->BDTDNE[stratmixingtype]);
          UnsrtELInsert(v, tail, sto->transferEL);
        }
      }

      for(int i = 1; i <= sto->transferEL->nedges; i++) {
        ToggleKnownEdge(sto->transferEL->tails[i], sto->transferEL->heads[i], sto->combined_nonBDTDNE, TRUE);
        ToggleKnownEdge(sto->transferEL->tails[i], sto->transferEL->heads[i], sto->combined_BDTDNE, FALSE);        
      }
    }
    
    if(sto->headmaxl) {
      // head will be newly submaxl after toggle, so add it to the appropriate node list    
      sto->nodesvec[sto->stratheadtype][sto->bdheadtype][sto->attrcounts[sto->stratheadtype][sto->bdheadtype]] = head;
      sto->nodepos[head] = sto->attrcounts[sto->stratheadtype][sto->bdheadtype];
      sto->attrcounts[sto->stratheadtype][sto->bdheadtype]++;

      sto->transferEL->nedges = 0;      
      
      STEP_THROUGH_OUTEDGES_NET(head, e, v, sto->combined_nonBDTDNE) {
        if(IN_DEG[v] + OUT_DEG[v] < sto->bound) {
          int stratmixingtype = sto->indmat[sto->stratheadtype][sto->strat_vattr[v]];
          UnsrtELDelete(head, v, sto->nonBDTDNE[stratmixingtype]);
          UnsrtELInsert(head, v, sto->BDTDNE[stratmixingtype]);
          UnsrtELInsert(head, v, sto->transferEL);
        }
      }
      
      STEP_THROUGH_INEDGES_NET(head, e, v, sto->combined_nonBDTDNE) {
        if(IN_DEG[v] + OUT_DEG[v] < sto->bound) {
          int stratmixingtype = sto->indmat[sto->stratheadtype][sto->strat_vattr[v]];
          UnsrtELDelete(v, head, sto->nonBDTDNE[stratmixingtype]);
          UnsrtELInsert(v, head, sto->BDTDNE[stratmixingtype]);
          UnsrtELInsert(v, head, sto->transferEL);
        }
      }

      for(int i = 1; i <= sto->transferEL->nedges; i++) {
        ToggleKnownEdge(sto->transferEL->tails[i], sto->transferEL->heads[i], sto->combined_nonBDTDNE, TRUE);
        ToggleKnownEdge(sto->transferEL->tails[i], sto->transferEL->heads[i], sto->combined_BDTDNE, FALSE);        
      }
    }
  } else {
    if(sto->tailmaxl) {
      // tail will be newly maxl after toggle, so remove it from the appropriate node list
      sto->nodesvec[sto->strattailtype][sto->bdtailtype][sto->nodepos[tail]] = sto->nodesvec[sto->strattailtype][sto->bdtailtype][sto->attrcounts[sto->strattailtype][sto->bdtailtype] - 1];
      sto->nodepos[sto->nodesvec[sto->strattailtype][sto->bdtailtype][sto->nodepos[tail]]] = sto->nodepos[tail];
      sto->attrcounts[sto->strattailtype][sto->bdtailtype]--;

      sto->transferEL->nedges = 0;

      STEP_THROUGH_OUTEDGES_NET(tail, e, v, sto->combined_BDTDNE) {
        int stratmixingtype = sto->indmat[sto->strattailtype][sto->strat_vattr[v]];
        UnsrtELDelete(tail, v, sto->BDTDNE[stratmixingtype]);
        UnsrtELInsert(tail, v, sto->nonBDTDNE[stratmixingtype]);
        UnsrtELInsert(tail, v, sto->transferEL);        
      }

      STEP_THROUGH_INEDGES_NET(tail, e, v, sto->combined_BDTDNE) {
        int stratmixingtype = sto->indmat[sto->strattailtype][sto->strat_vattr[v]];
        UnsrtELDelete(v, tail, sto->BDTDNE[stratmixingtype]);
        UnsrtELInsert(v, tail, sto->nonBDTDNE[stratmixingtype]);
        UnsrtELInsert(v, tail, sto->transferEL);
      }

      for(int i = 1; i <= sto->transferEL->nedges; i++) {
        ToggleKnownEdge(sto->transferEL->tails[i], sto->transferEL->heads[i], sto->combined_BDTDNE, TRUE);
        ToggleKnownEdge(sto->transferEL->tails[i], sto->transferEL->heads[i], sto->combined_nonBDTDNE, FALSE);        
      }
    }
    
    if(sto->headmaxl) {
      // head will be newly maxl after toggle, so remove it from the appropriate node list
      sto->nodesvec[sto->stratheadtype][sto->bdheadtype][sto->nodepos[head]] = sto->nodesvec[sto->stratheadtype][sto->bdheadtype][sto->attrcounts[sto->stratheadtype][sto->bdheadtype] - 1];
      sto->nodepos[sto->nodesvec[sto->stratheadtype][sto->bdheadtype][sto->nodepos[head]]] = sto->nodepos[head];
      sto->attrcounts[sto->stratheadtype][sto->bdheadtype]--;

      sto->transferEL->nedges = 0;

      STEP_THROUGH_OUTEDGES_NET(head, e, v, sto->combined_BDTDNE) {
        int stratmixingtype = sto->indmat[sto->stratheadtype][sto->strat_vattr[v]];
        UnsrtELDelete(head, v, sto->BDTDNE[stratmixingtype]);
        UnsrtELInsert(head, v, sto->nonBDTDNE[stratmixingtype]);
        UnsrtELInsert(head, v, sto->transferEL);        
      }

      STEP_THROUGH_INEDGES_NET(head, e, v, sto->combined_BDTDNE) {
        int stratmixingtype = sto->indmat[sto->stratheadtype][sto->strat_vattr[v]];
        UnsrtELDelete(v, head, sto->BDTDNE[stratmixingtype]);
        UnsrtELInsert(v, head, sto->nonBDTDNE[stratmixingtype]);
        UnsrtELInsert(v, head, sto->transferEL);
      }

      for(int i = 1; i <= sto->transferEL->nedges; i++) {
        ToggleKnownEdge(sto->transferEL->tails[i], sto->transferEL->heads[i], sto->combined_BDTDNE, TRUE);
        ToggleKnownEdge(sto->transferEL->tails[i], sto->transferEL->heads[i], sto->combined_nonBDTDNE, FALSE);        
      }
    }
  }
}

MH_F_FN(Mf_discordBDStratTNT) {
  GET_STORAGE(discordBDStratTNTStorage, sto);

  // let BDStratTNT's F_FN do most of the work
  MH_STORAGE = sto->static_sto;
  Mf_BDStratTNT(MHp, nwp);
  MH_STORAGE = sto;
  Free(sto->static_sto);

  for(int i = 0; i < sto->nmixtypes; i++) {
    UnsrtELDestroy(sto->BDTDNE[i]);
    UnsrtELDestroy(sto->nonBDTDNE[i]);    
    UnsrtELDestroy(sto->discordantEdges[i]);
  }
  Free(sto->BDTDNE);
  Free(sto->nonBDTDNE);
  Free(sto->discordantEdges);
  
  NetworkDestroy(sto->combined_BDTDNE);
  NetworkDestroy(sto->combined_nonBDTDNE);
  UnsrtELDestroy(sto->transferEL);
  // MHp->storage itself should be Freed by MHProposalDestroy
}


#undef OUTVAL_NET
#undef INVAL_NET
#undef MIN_OUTEDGE_NET
#undef MIN_INEDGE_NET
#undef NEXT_OUTEDGE_NET
#undef NEXT_INEDGE_NET
#undef STEP_THROUGH_OUTEDGES_NET
#undef STEP_THROUGH_INEDGES_NET
