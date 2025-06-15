/*  File src/discordTNT.c in package tergm, part of the Statnet suite of
 *  packages for network analysis, https://statnet.org .
 *
 *  This software is distributed under the GPL-3 license.  It is free, open
 *  source, and has the attribution requirements (GPL Section 7) at
 *  https://statnet.org/attribution .
 *
 *  Copyright 2008-2025 Statnet Commons
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
    for(int i = 1; i <= UnsrtELSize(sto->discordantEdges); i++) {
      UnsrtELInsert(sto->discordantEdges->tails[i], sto->discordantEdges->heads[i], sto->nonDiscordantEdges);
    }      
      
    // "clear" the discordant dyads
    UnsrtELSize(sto->discordantEdges) = 0;
    UnsrtELSize(sto->discordantNonEdges) = 0;
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
      if(unif_rand() < UnsrtELSize(sto->nonDiscordantEdges)/((double) nedges)) {
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
    if(unif_rand() < UnsrtELSize(sto->discordantEdges)/((double) nddyads)) {
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

#define OUTVAL_NET(e,n) ((n)->outedges[(e)].value)
#define INVAL_NET(e,n) ((n)->inedges[(e)].value)
#define MIN_OUTEDGE_NET(a,n) (EdgetreeMinimum((n)->outedges, (a)))
#define MIN_INEDGE_NET(a,n) (EdgetreeMinimum((n)->inedges, (a)))
#define NEXT_OUTEDGE_NET(e,n) (EdgetreeSuccessor((n)->outedges,(e)))
#define NEXT_INEDGE_NET(e,n) (EdgetreeSuccessor((n)->inedges,(e)))

#define EXEC_THROUGH_EDGES_EATH_NET_DECL(node, ego, alter, tail, head, edge, net, subroutine) { \
  Vertex (ego) = (node); \
  Vertex (alter), (tail), (head); \
  Edge (edge); \
  (tail) = (ego); \
  for((edge) = MIN_OUTEDGE_NET((ego),(net)); ((head) = (alter) = OUTVAL_NET((edge),(net))) != 0; (edge) = NEXT_OUTEDGE_NET((edge),(net))) { \
    subroutine \
  } \
  (head) = (ego); \
  for((edge) = MIN_INEDGE_NET((ego),(net)); ((tail) = (alter) = INVAL_NET((edge),(net))) != 0; (edge) = NEXT_INEDGE_NET((edge),(net))) { \
    subroutine \
  } \
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
  HashEL **edges;
  
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
  sto->edges = static_sto->hash;
  sto->static_sto = static_sto;
  
  sto->nodes = R_Calloc(2, Vertex);
  sto->maxl = R_Calloc(2, int);
  
  sto->BDTDNE = R_Calloc(sto->static_sto->strat_nmixtypes, HashEL *);
  sto->discordantEdges = R_Calloc(sto->static_sto->strat_nmixtypes, HashEL *);  
  for(int i = 0; i < sto->static_sto->strat_nmixtypes; i++) {
    sto->BDTDNE[i] = HashELInitialize(0, NULL, NULL, FALSE);
    sto->discordantEdges[i] = HashELInitialize(0, NULL, NULL, FALSE);
  }
  sto->combined_BDTDNE = NetworkInitialize(NULL, NULL, 0, N_NODES, DIRECTED, BIPARTITE);
  sto->combined_nonBDTDNE = NetworkInitialize(NULL, NULL, 0, N_NODES, DIRECTED, BIPARTITE);
  sto->transferEL = UnsrtELInitialize(0, NULL, NULL, FALSE);

  sto->discordance_fraction = asReal(getListElement(MHp->R, "discordance_fraction"));  
}

MH_X_FN(Mx_discordBDStratTNT) {
  if(type == TICK) {
    GET_STORAGE(discordBDStratTNTStorage, sto);

    for(int i = 0; i < sto->static_sto->strat_nmixtypes; i++) {
      // clear all the discordance information
      if(HashELSize(sto->BDTDNE[i]) > 0) HashELClear(sto->BDTDNE[i]);
      if(HashELSize(sto->discordantEdges[i]) > 0) HashELClear(sto->discordantEdges[i]);      
    }
    
    // for now, destroy and recreate each time step (can we do this more efficiently?)
    NetworkDestroy(sto->combined_BDTDNE);
    NetworkDestroy(sto->combined_nonBDTDNE);
    sto->combined_BDTDNE = NetworkInitialize(NULL, NULL, 0, N_NODES, DIRECTED, BIPARTITE);
    sto->combined_nonBDTDNE = NetworkInitialize(NULL, NULL, 0, N_NODES, DIRECTED, BIPARTITE);
  }
}

MH_P_FN(MH_discordBDStratTNT) {
  GET_STORAGE(discordBDStratTNTStorage, sto);
  GET_AUX_STORAGE(StoreTimeAndLasttoggle, dur_inf);

  // sample a toggleable strat mixing type on which to make a proposal
  int strat_i = WtPopGetRand(sto->static_sto->wtp);
  sto->static_sto->stratmixingtype = strat_i;
  
  // number of edges of this mixing type
  int nedgestype = HashELSize(sto->edges[strat_i]);
  
  Dyad ndyadstype = BDStratBlocksDyadCount(sto->static_sto->blocks, strat_i);
  
  int nddyadstype = HashELSize(sto->discordantEdges[strat_i]) + HashELSize(sto->BDTDNE[strat_i]);
  
  int in_network;
  int in_discord;
  
  if(nddyadstype == 0 || unif_rand() < 1 - sto->discordance_fraction) {
    // propose from network
    if((unif_rand() < 0.5 && nedgestype > 0) || ndyadstype == 0) {
      // propose toggling off an existing edge of strat mixing type strat_i
      HashELGetRand(Mtail, Mhead, sto->edges[strat_i]);

      in_network = TRUE;
      in_discord = kh_get(DyadMapInt, dur_inf->discord, TH(Mtail[0], Mhead[0])) != kh_none;
    } else {
      // select a random BD toggleable dyad of strat mixing type strat_i and propose toggling it
      BDStratBlocksGetRandWithCount(Mtail, Mhead, sto->static_sto->blocks, strat_i, ndyadstype);
         
      in_network = IS_OUTEDGE(Mtail[0],Mhead[0]);
      in_discord = kh_get(DyadMapInt, dur_inf->discord, TH(Mtail[0], Mhead[0])) != kh_none;
    }
  } else {
    // propose from discord
    if(unif_rand() < ((double) HashELSize(sto->discordantEdges[strat_i]))/nddyadstype) {
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

  // compute proposed dyad count for current mixing type (only)
  Dyad proposeddyadstype = BDStratBlocksDyadCountOnToggle(*Mtail, *Mhead, sto->static_sto->blocks, strat_i, sto->static_sto->tailmaxl, sto->static_sto->headmaxl);

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
  double prob_weight = sto->static_sto->current_total_weight/sto->static_sto->proposed_total_weight;
  
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
  if(sto->static_sto->strat_nmixtypestoupdate > 0) {
    sto->static_sto->current_total_weight = sto->static_sto->proposed_total_weight;
    for(int i = 0; i < sto->static_sto->strat_nmixtypestoupdate; i++) {
      WtPopSetWt(sto->static_sto->strat_mixtypestoupdate[i], edgestate ? sto->static_sto->original_weights[sto->static_sto->strat_mixtypestoupdate[i]] : 0, sto->static_sto->wtp);          
    }
  }

  // add or remove the dyad being toggled from the relevant edge set(s)/network
  HashELToggleKnown(tail, head, sto->edges[sto->static_sto->stratmixingtype], edgestate);
  if(sto->in_discord == edgestate) {
    HashELToggleKnown(tail, head, sto->discordantEdges[sto->static_sto->stratmixingtype], edgestate);      
  } else {
    HashELToggleKnown(tail, head, sto->BDTDNE[sto->static_sto->stratmixingtype], !edgestate);
    ToggleKnownEdge(tail, head, sto->combined_BDTDNE, !edgestate);      
  }

  BDNodeListsToggleIf(tail, head, sto->static_sto->lists, sto->static_sto->tailmaxl, sto->static_sto->headmaxl);

  int tailattr = sto->static_sto->bd_vattr[tail];
  int headattr = sto->static_sto->bd_vattr[head];
  
  // update dyad toggleability statuses, as appropriate  
  Network *relevant_net = edgestate ? sto->combined_nonBDTDNE : sto->combined_BDTDNE;
  UnsrtELSize(sto->transferEL) = 0; // reset transferEL
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
  for(int i = 1; i <= UnsrtELSize(sto->transferEL); i++) {
    ToggleKnownEdge(sto->transferEL->tails[i], sto->transferEL->heads[i], sto->combined_nonBDTDNE, edgestate);
    ToggleKnownEdge(sto->transferEL->tails[i], sto->transferEL->heads[i], sto->combined_BDTDNE, !edgestate);        
  }

  sto->static_sto->indegree[tailattr][head] += edgestate ? -1 : 1;
  sto->static_sto->outdegree[headattr][tail] += edgestate ? -1 : 1;  
}

MH_F_FN(Mf_discordBDStratTNT) {
  GET_STORAGE(discordBDStratTNTStorage, sto);
  // free things used only in the dynamic proposal
  for(int i = 0; i < sto->static_sto->strat_nmixtypes; i++) {
    HashELDestroy(sto->BDTDNE[i]);
    HashELDestroy(sto->discordantEdges[i]);
  }
  R_Free(sto->BDTDNE);
  R_Free(sto->discordantEdges);
  
  R_Free(sto->nodes);
  R_Free(sto->maxl);
  
  NetworkDestroy(sto->combined_BDTDNE);
  NetworkDestroy(sto->combined_nonBDTDNE);
  UnsrtELDestroy(sto->transferEL);
  
  // let BDStratTNT's F_FN do most of the work
  MH_STORAGE = sto->static_sto;
  Mf_BDStratTNT(MHp, nwp);
  R_Free(sto->static_sto);
  MH_STORAGE = sto;
  // MHp->storage itself should be R_Freed by MHProposalDestroy
}
