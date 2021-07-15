/*  File src/discordTNT.c in package tergm, part of the
 *  Statnet suite of packages for network analysis, https://statnet.org .
 *
 *  This software is distributed under the GPL-3 license.  It is free,
 *  open source, and has the attribution requirements (GPL Section 7) at
 *  https://statnet.org/attribution .
 *
 *  Copyright 2008-2021 Statnet Commons
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
#include "ergm_BDStratBlocks.h"
#include "ergm_BDStrat_proposals.h"
#include "ergm_hash_edgelist.h"



/********************
    MH_discordTNT
********************/

typedef struct {
  UnsrtEL *nonDiscordantEdges;
  UnsrtEL *discordantEdges;
  UnsrtEL *discordantNonEdges;
  
  double discordance_fraction;  
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
      in_discord = kh_get(DyadMapInt, dur_inf->discord, TH(Mtail[0], Mhead[0])) != kh_none;

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
}

MH_U_FN(Mu_discordTNT) {  
  GET_AUX_STORAGE(StoreTimeAndLasttoggle, dur_inf);
  GET_STORAGE(discordTNTStorage, sto);
  int in_discord = kh_get(DyadMapInt, dur_inf->discord, TH(tail, head)) != kh_none;
  
  // add or remove the dyad from the appropriate discordance edgelist  
  if(in_discord == edgestate) {
    UnsrtELToggleKnown(tail, head, sto->discordantEdges, edgestate);
  } else {
    UnsrtELToggleKnown(tail, head, sto->nonDiscordantEdges, edgestate);
    UnsrtELToggleKnown(tail, head, sto->discordantNonEdges, !edgestate);
  }
}

MH_F_FN(Mf_discordTNT) {
  GET_STORAGE(discordTNTStorage, sto);
  
  UnsrtELDestroy(sto->nonDiscordantEdges);
  UnsrtELDestroy(sto->discordantNonEdges);
  UnsrtELDestroy(sto->discordantEdges);
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

  double discordance_fraction;    
    
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

  sto->discordance_fraction = asReal(getListElement(MHp->R, "discordance_fraction"));  
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
  
  Dyad ndyadstype = BDStratBlocksDyadCount(sto->static_sto->blocks, strat_i);
  
  int nddyadstype = sto->discordantEdges[strat_i]->list->nedges + sto->BDTDNE[strat_i]->list->nedges;
  
  int in_network;
  int in_discord;
  
  if(nddyadstype == 0 || unif_rand() < 1 - sto->discordance_fraction) {
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
      BDStratBlocksGetRandWithCount(Mtail, Mhead, sto->static_sto->blocks, strat_i, ndyadstype);
         
      in_network = IS_OUTEDGE(Mtail[0],Mhead[0]);
      in_discord = kh_get(DyadMapInt, dur_inf->discord, TH(Mtail[0], Mhead[0])) != kh_none;
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

  int tailattr = sto->static_sto->bd_vattr[*Mtail];
  int headattr = sto->static_sto->bd_vattr[*Mhead];
  
  sto->static_sto->tailmaxl = (DIRECTED ? sto->static_sto->outdegree[headattr][*Mtail] : sto->static_sto->indegree[headattr][*Mtail] + sto->static_sto->outdegree[headattr][*Mtail]) == sto->static_sto->maxout[headattr][*Mtail] - 1 + in_network;
  sto->static_sto->headmaxl = (DIRECTED ? sto->static_sto->indegree[tailattr][*Mhead] : sto->static_sto->indegree[tailattr][*Mhead] + sto->static_sto->outdegree[tailattr][*Mhead]) == sto->static_sto->maxin[tailattr][*Mhead] - 1 + in_network;

  sto->nodes[0] = *Mtail;  
  sto->nodes[1] = *Mhead;
        
  sto->maxl[0] = sto->static_sto->tailmaxl;  
  sto->maxl[1] = sto->static_sto->headmaxl;  

  BDStratBlocksSetLast(*Mtail, *Mhead, in_network, sto->static_sto->blocks);

  // compute proposed dyad count for current mixing type (only)
  Dyad proposeddyadstype = BDStratBlocksDyadCountOnToggle(*Mtail, *Mhead, sto->static_sto->blocks, strat_i, in_network ? +1 : -1, sto->static_sto->tailmaxl, sto->static_sto->headmaxl);

  ComputeChangesToToggleability(Mtail, Mhead, sto->static_sto);
  
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
        if(// not the proposal dyad, and
           (tail != *Mtail || head != *Mhead) &&
           // same strat mixing type as the proposal dyad, and
           sto->static_sto->indmat[sto->static_sto->strat_vattr[tail]][sto->static_sto->strat_vattr[head]] == strat_i && 
           // same bd mixing type as the proposal dyad, and
           ((tailattr == sto->static_sto->bd_vattr[tail] && headattr == sto->static_sto->bd_vattr[head]) || 
            (!DIRECTED && tailattr == sto->static_sto->bd_vattr[head] && headattr == sto->static_sto->bd_vattr[tail])) &&
           // respects directedness, and
           (!DIRECTED || tail == *Mtail || head == *Mhead) &&
           // will change toggleability status if we accept the proposed toggle
           (!in_network || (DIRECTED ? (alter == tail ? sto->static_sto->outdegree[sto->static_sto->bd_vattr[ego]][alter] < sto->static_sto->maxout[sto->static_sto->bd_vattr[ego]][alter] 
                                                      : sto->static_sto->indegree[sto->static_sto->bd_vattr[ego]][alter] < sto->static_sto->maxin[sto->static_sto->bd_vattr[ego]][alter]) 
                                     : (sto->static_sto->indegree[sto->static_sto->bd_vattr[ego]][alter] + sto->static_sto->outdegree[sto->static_sto->bd_vattr[ego]][alter] < sto->static_sto->maxout[sto->static_sto->bd_vattr[ego]][alter])))) {
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
  
  double forward = (nddyadstype == 0) ? forward_network : (sto->discordance_fraction*forward_discord + (1 - sto->discordance_fraction)*forward_network);
  double backward = (propnddyadstype == 0) ? backward_network : (sto->discordance_fraction*backward_discord + (1 - sto->discordance_fraction)*backward_network);

  MHp->logratio = log(prob_weight*backward/forward);
}

MH_U_FN(Mu_discordBDStratTNT) {   
  GET_STORAGE(discordBDStratTNTStorage, sto);

  // if any strat mixing types have changed toggleability status, update prob info accordingly
  if(sto->static_sto->nmixtypestoupdate > 0) {
    sto->static_sto->currentcumprob = sto->static_sto->proposedcumprob;
    for(int i = 0; i < sto->static_sto->nmixtypestoupdate; i++) {
      WtPopSetWt(sto->static_sto->mixtypestoupdate[i], edgestate ? sto->static_sto->originalprobvec[sto->static_sto->mixtypestoupdate[i]] : 0, sto->static_sto->wtp);          
    }
  }

  // add or remove the dyad being toggled from the relevant edge set(s)/network
  if(sto->in_discord == edgestate) {
    HashELToggleKnown(tail, head, sto->discordantEdges[sto->static_sto->stratmixingtype], edgestate);      
  } else {
    HashELToggleKnown(tail, head, sto->nonDiscordantEdges[sto->static_sto->stratmixingtype], edgestate);
    HashELToggleKnown(tail, head, sto->BDTDNE[sto->static_sto->stratmixingtype], !edgestate);
    ToggleKnownEdge(tail, head, sto->combined_BDTDNE, !edgestate);      
  }

  BDStratBlocksToggleIf(tail, head, sto->static_sto->blocks, sto->static_sto->tailmaxl, sto->static_sto->headmaxl);

  int tailattr = sto->static_sto->bd_vattr[tail];
  int headattr = sto->static_sto->bd_vattr[head];
  
  // update dyad toggleability statuses, as appropriate  
  Network *relevant_net = edgestate ? sto->combined_nonBDTDNE : sto->combined_BDTDNE;
  sto->transferEL->nedges = 0; // reset transferEL
  for(int i = 0; i < 2; i++) {
    if(sto->maxl[i]) {
      EXEC_THROUGH_EDGES_EATH_NET_DECL(sto->nodes[i], ego, alter, _tail, _head, edge, relevant_net, {
        if(// same bd mixing type as the proposal dyad, and
           ((tailattr == sto->static_sto->bd_vattr[_tail] && headattr == sto->static_sto->bd_vattr[_head]) || 
            (!DIRECTED && tailattr == sto->static_sto->bd_vattr[_head] && headattr == sto->static_sto->bd_vattr[_tail])) &&
           // respects directedness, and
           (!DIRECTED || _tail == tail || _head == head) &&            
           // will change toggleability status if we accept the proposed toggle
           (!edgestate || (DIRECTED ? (alter == _tail ? sto->static_sto->outdegree[sto->static_sto->bd_vattr[ego]][alter] < sto->static_sto->maxout[sto->static_sto->bd_vattr[ego]][alter] 
                                                      : sto->static_sto->indegree[sto->static_sto->bd_vattr[ego]][alter] < sto->static_sto->maxin[sto->static_sto->bd_vattr[ego]][alter]) 
                                    : (sto->static_sto->indegree[sto->static_sto->bd_vattr[ego]][alter] + sto->static_sto->outdegree[sto->static_sto->bd_vattr[ego]][alter] < sto->static_sto->maxout[sto->static_sto->bd_vattr[ego]][alter])))) {
          int stratmixingtype = sto->static_sto->indmat[sto->static_sto->strat_vattr[_tail]][sto->static_sto->strat_vattr[_head]];
          HashELToggleKnown(_tail, _head, sto->BDTDNE[stratmixingtype], !edgestate);
          UnsrtELInsert(_tail, _head, sto->transferEL);
        }
      });
    }
  } 

  // apply changes in transferEL to the Network objects
  for(int i = 1; i <= sto->transferEL->nedges; i++) {
    ToggleKnownEdge(sto->transferEL->tails[i], sto->transferEL->heads[i], sto->combined_nonBDTDNE, edgestate);
    ToggleKnownEdge(sto->transferEL->tails[i], sto->transferEL->heads[i], sto->combined_BDTDNE, !edgestate);        
  }

  sto->static_sto->indegree[tailattr][head] += edgestate ? -1 : 1;
  sto->static_sto->outdegree[headattr][tail] += edgestate ? -1 : 1;  
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
