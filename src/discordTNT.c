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
#include "ergm_nodelist.h"
#include "ergm_BDStrat_proposals.h"
#include "ergm_hash_edgelist.h"

typedef struct {
  UnsrtEL *nonDiscordantEdges;
  UnsrtEL *discordantEdges;
  UnsrtEL *discordantNonEdges;
  
  double discordance_fraction;
  
  int in_discord;
  
} discordTNTStorage; 

MH_I_FN(Mi_discordTNT) {
  MHp->ntoggles = 1;
  
  ALLOC_STORAGE(1, discordTNTStorage, sto);
  sto->nonDiscordantEdges = UnsrtELInitialize(0, NULL, NULL, FALSE);
  sto->discordantEdges = UnsrtELInitialize(0, NULL, NULL, FALSE);
  sto->discordantNonEdges = UnsrtELInitialize(0, NULL, NULL, FALSE);

  sto->discordance_fraction = asReal(getListElement(MHp->R, "discordance_fraction"));  

  EXEC_THROUGH_NET_EDGES(tail, head, e, {
    UnsrtELInsert(tail, head, sto->nonDiscordantEdges);
  });
}

MH_X_FN(Mx_discordTNT) {
  GET_STORAGE(discordTNTStorage, sto);
  
  if(type == TICK) {
    // transfer discordant edges to non-discordant edges
    for(int i = 1; i <= sto->discordantEdges->nedges; i++) {
      UnsrtELInsert(sto->discordantEdges->tails[i], sto->discordantEdges->heads[i], sto->nonDiscordantEdges);
    }      
      
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

      if(in_discord) {
        // need to resample to know index
        if(in_network) {
          UnsrtELGetRand(Mtail, Mhead, sto->discordantEdges);
        } else {
          UnsrtELGetRand(Mtail, Mhead, sto->discordantNonEdges);
        }
      } else if(in_network) {
        UnsrtELGetRand(Mtail, Mhead, sto->nonDiscordantEdges);
      }
    } else {
      // propose toggling off an edge in network
      if(unif_rand() < sto->nonDiscordantEdges->nedges/((double) nedges)) {
        UnsrtELGetRand(Mtail, Mhead, sto->nonDiscordantEdges);
        in_discord = FALSE;      
      } else {
        UnsrtELGetRand(Mtail, Mhead, sto->discordantEdges);
        in_discord = TRUE;
      }
      
      in_network = TRUE;
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
  GET_STORAGE(discordTNTStorage, sto);
  
  // add or remove the dyad from the appropriate discordance edgelist  
  if(sto->in_discord == edgeflag) {
    UnsrtELToggleKnown(tail, head, sto->discordantEdges, edgeflag);
  } else {
    UnsrtELToggleKnown(tail, head, sto->nonDiscordantEdges, edgeflag);
    UnsrtELToggleKnown(tail, head, sto->discordantNonEdges, !edgeflag);
  }
}

MH_F_FN(Mf_discordTNT) {
  GET_STORAGE(discordTNTStorage, sto);
  
  UnsrtELDestroy(sto->nonDiscordantEdges);
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
    
  int in_discord;

  StratTNTStorage *static_sto;
} discordStratTNTStorage; 


MH_I_FN(Mi_discordStratTNT) {
  // let StratTNT's I_FN do most of the work
  Mi_StratTNT(MHp, nwp);

  // do a little copying and handle purely temporal aspects
  StratTNTStorage *static_sto = MH_STORAGE;
  ALLOC_STORAGE(1, discordStratTNTStorage, sto);
  sto->static_sto = static_sto;
  sto->nonDiscordantELs = static_sto->els;

  sto->discordantELs = Calloc(static_sto->nmixtypes, UnsrtEL *);
  sto->discordantNonELs = Calloc(static_sto->nmixtypes, UnsrtEL *);
  for(int i = 0; i < static_sto->nmixtypes; i++) {
    sto->discordantELs[i] = UnsrtELInitialize(0, NULL, NULL, FALSE);
    sto->discordantNonELs[i] = UnsrtELInitialize(0, NULL, NULL, FALSE);
  }
} 

MH_X_FN(Mx_discordStratTNT) {
  GET_STORAGE(discordStratTNTStorage, sto);
  
  if(type == TICK) {
    // transfer discordant edges to nondiscordant edges
    // clear discordant edges and discordant nonedges
    for(int i = 0; i < sto->static_sto->nmixtypes; i++) {
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
  sto->static_sto->currentmixingtype = WtPopGetRand(sto->static_sto->wtp);
        
  // number of edges of this mixing type
  int nedgestype = sto->nonDiscordantELs[sto->static_sto->currentmixingtype]->nedges + sto->discordantELs[sto->static_sto->currentmixingtype]->nedges;

  // number of dyads of this mixing type
  Dyad ndyadstype = sto->static_sto->ndyadstype[sto->static_sto->currentmixingtype];
  
  // number of discordant dyads of this mixing type
  int nddyadstype = sto->discordantNonELs[sto->static_sto->currentmixingtype]->nedges + sto->discordantELs[sto->static_sto->currentmixingtype]->nedges;
  
  // flags
  int in_discord;
  int in_network;
  
  if(nddyadstype == 0 || unif_rand() < 0.5) {
    // propose from network
    if(nedgestype == 0 || unif_rand() < 0.5) {
      // propose toggling a random dyad of the specified mixing type
      NodeListGetRandWithCount(Mtail, Mhead, sto->static_sto->nodelist, sto->static_sto->currentmixingtype, ndyadstype);
      
      in_network = IS_OUTEDGE(Mtail[0], Mhead[0]);
      in_discord = kh_get(DyadMapInt, dur_inf->discord, THKey(dur_inf->discord, Mtail[0], Mhead[0])) != kh_none;
      
      // if it resides in any of the edgelists we store, we need to resample
      if(in_network) {
        if(in_discord) {
          UnsrtELGetRand(Mtail, Mhead, sto->discordantELs[sto->static_sto->currentmixingtype]);
        } else {
          UnsrtELGetRand(Mtail, Mhead, sto->nonDiscordantELs[sto->static_sto->currentmixingtype]);
        }
      } else if(in_discord) {
        UnsrtELGetRand(Mtail, Mhead, sto->discordantNonELs[sto->static_sto->currentmixingtype]);
      }
    } else {
      // propose toggling off an edge of the specified mixing type
      if(unif_rand() < sto->nonDiscordantELs[sto->static_sto->currentmixingtype]->nedges/((double) nedgestype)) {
        UnsrtELGetRand(Mtail, Mhead, sto->nonDiscordantELs[sto->static_sto->currentmixingtype]);
        in_discord = FALSE;
      } else {
        UnsrtELGetRand(Mtail, Mhead, sto->discordantELs[sto->static_sto->currentmixingtype]);
        in_discord = TRUE;
      }
      in_network = TRUE;
    }
  } else {
    // propose from discord
    if(unif_rand() < sto->discordantELs[sto->static_sto->currentmixingtype]->nedges/((double) nddyadstype)) {
      UnsrtELGetRand(Mtail, Mhead, sto->discordantELs[sto->static_sto->currentmixingtype]);
      in_network = TRUE;
    } else {
      UnsrtELGetRand(Mtail, Mhead, sto->discordantNonELs[sto->static_sto->currentmixingtype]);
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
  GET_STORAGE(discordStratTNTStorage, sto);

  // add or remove edge from appropriate edgelist(s)
  if(sto->in_discord == edgeflag) {
    UnsrtELToggleKnown(tail, head, sto->discordantELs[sto->static_sto->currentmixingtype], edgeflag);      
  } else {
    UnsrtELToggleKnown(tail, head, sto->nonDiscordantELs[sto->static_sto->currentmixingtype], edgeflag);
    UnsrtELToggleKnown(tail, head, sto->discordantNonELs[sto->static_sto->currentmixingtype], !edgeflag);      
  }
}

MH_F_FN(Mf_discordStratTNT) {
  GET_STORAGE(discordStratTNTStorage, sto);
  // free things used only in the dynamic proposal
  for(int i = 0; i < sto->static_sto->nmixtypes; i++) {
    UnsrtELDestroy(sto->discordantELs[i]);
    UnsrtELDestroy(sto->discordantNonELs[i]);    
  }

  Free(sto->discordantELs);
  Free(sto->discordantNonELs);

  // let StratTNT's F_FN do most of the work
  MH_STORAGE = sto->static_sto;
  Mf_StratTNT(MHp, nwp);
  Free(sto->static_sto);
  MH_STORAGE = sto;
  // MHp->storage itself should be Freed by MHProposalDestroy
}

/********************
    MH_discordBDStratTNT
********************/

typedef struct {
  Network *combined_BDTDNE;
  Network *combined_nonBDTDNE;

  UnsrtEL *transferEL;

  HashEL **BDTDNE;
  
  HashEL **discordantEdges;
  HashEL **nonDiscordantEdges;
  
  Vertex *nodes;
  int *maxl;
  
  int in_discord;
    
  BDStratTNTStorage *static_sto;
} discordBDStratTNTStorage;

MH_I_FN(Mi_discordBDStratTNT) {
  // let BDStratTNT's I_FN do most of the work
  Mi_BDStratTNT(MHp, nwp);
  
  // copy a few things and handle purely temporal aspects
  BDStratTNTStorage *static_sto = MH_STORAGE;
  ALLOC_STORAGE(1, discordBDStratTNTStorage, sto);
  sto->nonDiscordantEdges = static_sto->hash;
  sto->static_sto = static_sto;
  
  sto->nodes = Calloc(2, Vertex);
  sto->maxl = Calloc(2, int);
  
  sto->BDTDNE = Calloc(sto->static_sto->nmixtypes, HashEL *);
  sto->discordantEdges = Calloc(sto->static_sto->nmixtypes, HashEL *);  
  for(int i = 0; i < sto->static_sto->nmixtypes; i++) {
    sto->BDTDNE[i] = HashELInitialize(0, NULL, NULL, FALSE, DIRECTED);
    sto->discordantEdges[i] = HashELInitialize(0, NULL, NULL, FALSE, DIRECTED);
  }
  sto->combined_BDTDNE = NetworkInitialize(NULL, NULL, 0, N_NODES, DIRECTED, BIPARTITE, FALSE, 0, NULL);
  sto->combined_nonBDTDNE = NetworkInitialize(NULL, NULL, 0, N_NODES, DIRECTED, BIPARTITE, FALSE, 0, NULL);
  sto->transferEL = UnsrtELInitialize(0, NULL, NULL, FALSE);
}

MH_X_FN(Mx_discordBDStratTNT) {
  if(type == TICK) {
    GET_STORAGE(discordBDStratTNTStorage, sto);

    for(int i = 0; i < sto->static_sto->nmixtypes; i++) {
      // transfer discordantEdges to nonDiscordantEdges
      for(int j = 1; j <= sto->discordantEdges[i]->list->nedges; j++) {
        HashELInsert(sto->discordantEdges[i]->list->tails[j], sto->discordantEdges[i]->list->heads[j], sto->nonDiscordantEdges[i]);
      }

      // clear all the discordance information
      if(sto->BDTDNE[i]->list->nedges > 0) HashELClear(sto->BDTDNE[i]);
      if(sto->discordantEdges[i]->list->nedges > 0) HashELClear(sto->discordantEdges[i]);      
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
  int strat_i = WtPopGetRand(sto->static_sto->wtp);
  sto->static_sto->stratmixingtype = strat_i;
  
  // number of edges of this mixing type
  int nedgestype = sto->nonDiscordantEdges[strat_i]->list->nedges + sto->discordantEdges[strat_i]->list->nedges;
  
  Dyad ndyadstype = NodeListDyadCount(sto->static_sto->nodelist, strat_i);
  
  int nddyadstype = sto->discordantEdges[strat_i]->list->nedges + sto->BDTDNE[strat_i]->list->nedges;
  
  int in_network;
  int in_discord;
  
  if(nddyadstype == 0 || unif_rand() < 0.5) {
    // propose from network
    if((unif_rand() < 0.5 && nedgestype > 0) || ndyadstype == 0) {
      // propose toggling off an existing edge of strat mixing type strat_i
      if(unif_rand() < ((double) sto->nonDiscordantEdges[strat_i]->list->nedges)/nedgestype) {
        HashELGetRand(Mtail, Mhead, sto->nonDiscordantEdges[strat_i]);
        in_discord = FALSE;
      } else {
        HashELGetRand(Mtail, Mhead, sto->discordantEdges[strat_i]);
        in_discord = TRUE;
      }
      in_network = TRUE;
    } else {
      // select a random BD toggleable dyad of strat mixing type strat_i and propose toggling it
      NodeListGetRandWithCount(Mtail, Mhead, sto->static_sto->nodelist, strat_i, ndyadstype);
         
      in_network = IS_OUTEDGE(Mtail[0],Mhead[0]);
      in_discord = kh_get(DyadMapInt, dur_inf->discord, THKey(dur_inf->discord, Mtail[0], Mhead[0])) != kh_none;
    }
  } else {
    // propose from discord
    if(unif_rand() < ((double) sto->discordantEdges[strat_i]->list->nedges)/nddyadstype) {
      HashELGetRand(Mtail, Mhead, sto->discordantEdges[strat_i]);
      in_network = TRUE;
    } else {
      HashELGetRand(Mtail, Mhead, sto->BDTDNE[strat_i]);
      in_network = FALSE;
    }
    in_discord = TRUE;
  }

  sto->in_discord = in_discord;

  sto->static_sto->tailmaxl = IN_DEG[Mtail[0]] + OUT_DEG[Mtail[0]] == sto->static_sto->bound - 1 + in_network;
  sto->static_sto->headmaxl = IN_DEG[Mhead[0]] + OUT_DEG[Mhead[0]] == sto->static_sto->bound - 1 + in_network;

  sto->nodes[0] = *Mtail;  
  sto->nodes[1] = *Mhead;
        
  sto->maxl[0] = sto->static_sto->tailmaxl;  
  sto->maxl[1] = sto->static_sto->headmaxl;  

  // compute proposed dyad count for current mixing type (only)
  Dyad proposeddyadstype = NodeListDyadCountOnToggle(*Mtail, *Mhead, sto->static_sto->nodelist, strat_i, in_network ? 1 : -1, sto->static_sto->tailmaxl, sto->static_sto->headmaxl);

  ComputeChangesToToggleability(Mtail, Mhead, in_network, sto->static_sto);
  
  // need to calculate number of BD-toggleable discordant dyads in the proposed network; 
  // this can involve both the proposal dyad itself and other dyads containing either of
  // the nodes in the proposal dyad
  int delta = in_network ? +1 : -1;

  int propnddyadstype = nddyadstype;  
  propnddyadstype += in_discord ? -1 : 1;
  
  Network *relevant_net = in_network ? sto->combined_nonBDTDNE : sto->combined_BDTDNE;
  for(int i = 0; i < 2; i++) {
    if(sto->maxl[i]) {
      EXEC_THROUGH_EDGES_EATH_NET_DECL(sto->nodes[i], ego, alter, tail, head, edge, relevant_net, {
        if(alter != sto->nodes[1 - i] && sto->static_sto->indmat[sto->static_sto->strat_vattr[tail]][sto->static_sto->strat_vattr[head]] == strat_i && (!in_network || IN_DEG[alter] + OUT_DEG[alter] < sto->static_sto->bound)) {
          propnddyadstype += delta;
        }
      });
    }
  }

  // calculate logratio
  double prob_weight = sto->static_sto->currentcumprob/sto->static_sto->proposedcumprob;
  
  double forward_network = in_network ? (ndyadstype == 0 ? 1.0/nedgestype : 0.5/nedgestype + (sto->static_sto->tailmaxl || sto->static_sto->headmaxl ? 0.0 : 0.5/ndyadstype)) : (nedgestype == 0 ? 1.0/ndyadstype : 0.5/ndyadstype);
  
  double forward_discord = in_discord ? 1.0/nddyadstype : 0;
  
  double backward_network = in_network ? (nedgestype == 1 ? 1.0/proposeddyadstype : 0.5/proposeddyadstype) : (proposeddyadstype == 0 ? 1.0/(nedgestype + 1) : 0.5/(nedgestype + 1) + (sto->static_sto->tailmaxl || sto->static_sto->headmaxl ? 0.0 : 0.5/proposeddyadstype));
  
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
  if(sto->static_sto->nmixtypestoupdate > 0) {
    sto->static_sto->currentcumprob = sto->static_sto->proposedcumprob;
    for(int i = 0; i < sto->static_sto->nmixtypestoupdate; i++) {
      WtPopSetWt(sto->static_sto->mixtypestoupdate[i], edgeflag ? sto->static_sto->originalprobvec[sto->static_sto->mixtypestoupdate[i]] : 0, sto->static_sto->wtp);          
    }
  }

  // add or remove the dyad being toggled from the relevant edge set(s)/network
  if(sto->in_discord == edgeflag) {
    HashELToggleKnown(tail, head, sto->discordantEdges[sto->static_sto->stratmixingtype], edgeflag);      
  } else {
    HashELToggleKnown(tail, head, sto->nonDiscordantEdges[sto->static_sto->stratmixingtype], edgeflag);
    HashELToggleKnown(tail, head, sto->BDTDNE[sto->static_sto->stratmixingtype], !edgeflag);
    ToggleKnownEdge(tail, head, sto->combined_BDTDNE, !edgeflag);      
  }

  NodeListToggleKnownIf(tail, head, sto->static_sto->nodelist, !edgeflag, sto->static_sto->tailmaxl, sto->static_sto->headmaxl);
  
  // update dyad toggleability statuses, as appropriate  
  Network *relevant_net = edgeflag ? sto->combined_nonBDTDNE : sto->combined_BDTDNE;
  sto->transferEL->nedges = 0; // reset transferEL
  for(int i = 0; i < 2; i++) {
    if(sto->maxl[i]) {
      EXEC_THROUGH_EDGES_EATH_NET_DECL(sto->nodes[i], ego, alter, _tail, _head, edge, relevant_net, {
        if(!edgeflag || IN_DEG[alter] + OUT_DEG[alter] < sto->static_sto->bound) {
          int stratmixingtype = sto->static_sto->indmat[sto->static_sto->strat_vattr[_tail]][sto->static_sto->strat_vattr[_head]];
          HashELToggleKnown(_tail, _head, sto->BDTDNE[stratmixingtype], !edgeflag);
          UnsrtELInsert(_tail, _head, sto->transferEL);
        }
      });
    }
  } 

  // apply changes in transferEL to the Network objects
  for(int i = 1; i <= sto->transferEL->nedges; i++) {
    ToggleKnownEdge(sto->transferEL->tails[i], sto->transferEL->heads[i], sto->combined_nonBDTDNE, edgeflag);
    ToggleKnownEdge(sto->transferEL->tails[i], sto->transferEL->heads[i], sto->combined_BDTDNE, !edgeflag);        
  }
}

MH_F_FN(Mf_discordBDStratTNT) {
  GET_STORAGE(discordBDStratTNTStorage, sto);
  // free things used only in the dynamic proposal
  for(int i = 0; i < sto->static_sto->nmixtypes; i++) {
    HashELDestroy(sto->BDTDNE[i]);
    HashELDestroy(sto->discordantEdges[i]);
  }
  Free(sto->BDTDNE);
  Free(sto->discordantEdges);
  
  Free(sto->nodes);
  Free(sto->maxl);
  
  NetworkDestroy(sto->combined_BDTDNE);
  NetworkDestroy(sto->combined_nonBDTDNE);
  UnsrtELDestroy(sto->transferEL);
  
  // let BDStratTNT's F_FN do most of the work
  MH_STORAGE = sto->static_sto;
  Mf_BDStratTNT(MHp, nwp);
  Free(sto->static_sto);
  MH_STORAGE = sto;
  // MHp->storage itself should be Freed by MHProposalDestroy
}
