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
  GET_STORAGE(discordTNTStorage, sto);
  
  // add or remove the dyad from the appropriate discordance edgelist  
  if(sto->in_discord == edgeflag) {
    UnsrtELToggleKnown(tail, head, sto->discordantEdges, edgeflag);
  } else {
    UnsrtELToggleKnown(tail, head, sto->discordantNonEdges, !edgeflag);
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
  int i = WtPopGetRand(sto->static_sto->wtp);
  
  // record the mixing type of the toggle, in case it's needed in the U function later
  sto->static_sto->currentmixingtype = i;    
      
  // number of edges of this mixing type
  int nedgestype = sto->nonDiscordantELs[i]->nedges + sto->discordantELs[i]->nedges;

  // number of dyads of this mixing type
  Dyad ndyadstype = sto->static_sto->ndyadstype[i];
  
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
                           sto->static_sto->nodesbycode, // tails
                           sto->static_sto->nodesbycode, // heads
                           sto->static_sto->tailtypes + i, // tailattrs
                           sto->static_sto->headtypes + i, // headattrs
                           sto->static_sto->nodecountsbycode, // tailcounts
                           sto->static_sto->nodecountsbycode, // headcounts
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
    MH_discordBDTNT
********************/

typedef struct {
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

  // copy a few things and handle purely temporal aspects
  BDTNTStorage *static_sto = MH_STORAGE;
  ALLOC_STORAGE(1, discordBDTNTStorage, sto);
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
    if((unif_rand() < 0.5 && nedges > 0) || (sto->static_sto->currentdyads == 0)) {
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
                           sto->static_sto->nodesvec, // tails
                           sto->static_sto->nodesvec, // heads
                           sto->static_sto->tailtypes, // tailattrs
                           sto->static_sto->headtypes, // headattrs
                           sto->static_sto->attrcounts, // tailcounts
                           sto->static_sto->attrcounts, // headcounts
                           sto->static_sto->nmixtypes, // length
                           sto->static_sto->currentdyads, // dyadcount
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
  
  sto->static_sto->tailtype = sto->static_sto->vattr[Mtail[0]];
  sto->static_sto->headtype = sto->static_sto->vattr[Mhead[0]];    
  
  sto->static_sto->tailmaxl = IN_DEG[Mtail[0]] + OUT_DEG[Mtail[0]] == sto->static_sto->bound - 1 + in_network;
  sto->static_sto->headmaxl = IN_DEG[Mhead[0]] + OUT_DEG[Mhead[0]] == sto->static_sto->bound - 1 + in_network;   

  sto->static_sto->nodes[0] = *Mtail;  
  sto->static_sto->nodes[1] = *Mhead;
        
  sto->static_sto->maxl[0] = sto->static_sto->tailmaxl;  
  sto->static_sto->maxl[1] = sto->static_sto->headmaxl;

  // obtain the count of dyads that can be toggled in the "GetRandBDDyad" branch, in the proposed network
  NodeListToggleKnownIf(*Mtail, sto->static_sto->nodesvec[sto->static_sto->tailtype], sto->static_sto->nodepos, sto->static_sto->attrcounts + sto->static_sto->tailtype, !in_network, sto->static_sto->tailmaxl);
  NodeListToggleKnownIf(*Mhead, sto->static_sto->nodesvec[sto->static_sto->headtype], sto->static_sto->nodepos, sto->static_sto->attrcounts + sto->static_sto->headtype, !in_network, sto->static_sto->headmaxl);      
  sto->static_sto->proposeddyads = NodeListDyadCount(sto->static_sto->attrcounts, sto->static_sto->attrcounts, sto->static_sto->tailtypes, sto->static_sto->headtypes, sto->static_sto->nmixtypes, TRUE, DIRECTED);
  NodeListToggleKnownIf(*Mtail, sto->static_sto->nodesvec[sto->static_sto->tailtype], sto->static_sto->nodepos, sto->static_sto->attrcounts + sto->static_sto->tailtype, in_network, sto->static_sto->tailmaxl);
  NodeListToggleKnownIf(*Mhead, sto->static_sto->nodesvec[sto->static_sto->headtype], sto->static_sto->nodepos, sto->static_sto->attrcounts + sto->static_sto->headtype, in_network, sto->static_sto->headmaxl);      
  
  // need to calculate number of submaxl edges and number of BD-toggleable discordant dyads
  // in the proposed network; this can involve both the proposal dyad itself and other dyads
  // containing either of the nodes in the proposal dyad
  int delta = in_network ? +1 : -1;

  sto->static_sto->proposedsubmaxledges = sto->static_sto->currentsubmaxledges;  
  if(!sto->static_sto->tailmaxl && !sto->static_sto->headmaxl) {
    sto->static_sto->proposedsubmaxledges -= delta;
  }
    
  int propnddyads = nddyads;  
  propnddyads += in_discord ? -1 : 1;
  
  Network *relevant_net = in_network ? sto->combined_nonBDTDNE : sto->combined_BDTDNE;
  for(int i = 0; i < 2; i++) {
    if(sto->static_sto->maxl[i]) {
      EXEC_THROUGH_EDGES_EA_NET_DECL(sto->static_sto->nodes[i], ego, alter, edge, nwp, {
        if(alter != sto->static_sto->nodes[1 - i] && IN_DEG[alter] + OUT_DEG[alter] < sto->static_sto->bound) {
          sto->static_sto->proposedsubmaxledges += delta;
        }
      });
      EXEC_THROUGH_EDGES_EA_NET_DECL(sto->static_sto->nodes[i], ego, alter, edge, relevant_net, {
        if(alter != sto->static_sto->nodes[1 - i] && (!in_network || IN_DEG[alter] + OUT_DEG[alter] < sto->static_sto->bound)) {
          propnddyads += delta;
        }
      });
    }
  } 
  
  double forward_network = in_network ? (sto->static_sto->currentdyads == 0 ? 1.0/nedges : 0.5/nedges + (0.5/sto->static_sto->currentdyads)*((double)sto->static_sto->currentsubmaxledges/nedges)) : (nedges == 0 ? 1.0/sto->static_sto->currentdyads : 0.5/sto->static_sto->currentdyads);
  
  double forward_discord = in_discord ? 1.0/nddyads : 0;
  
  double backward_network = in_network ? (nedges == 1 ? 1.0/sto->static_sto->proposeddyads : 0.5/sto->static_sto->proposeddyads) : (sto->static_sto->proposeddyads == 0 ? 1.0/(nedges + 1) : 0.5/(nedges + 1) + (0.5/sto->static_sto->proposeddyads)*((double) sto->static_sto->proposedsubmaxledges/(nedges + 1)));
  
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
  
  // add or remove the dyad being toggled from the relevant edge set(s)/network
  if(sto->in_discord == edgeflag) {
    UnsrtELToggleKnown(tail, head, sto->discordantEdges, edgeflag);      
  } else {
    UnsrtELToggleKnown(tail, head, sto->nonDiscordantEdges, edgeflag);
    UnsrtELToggleKnown(tail, head, sto->BDTDNE, !edgeflag);
    ToggleKnownEdge(tail, head, sto->combined_BDTDNE, !edgeflag);      
  }

  // update node lists as needed
  NodeListToggleKnownIf(tail, sto->static_sto->nodesvec[sto->static_sto->tailtype], sto->static_sto->nodepos, sto->static_sto->attrcounts + sto->static_sto->tailtype, !edgeflag, sto->static_sto->tailmaxl);
  NodeListToggleKnownIf(head, sto->static_sto->nodesvec[sto->static_sto->headtype], sto->static_sto->nodepos, sto->static_sto->attrcounts + sto->static_sto->headtype, !edgeflag, sto->static_sto->headmaxl);      
  
  // update dyad toggleability statuses as appropriate
  Network *relevant_net = edgeflag ? sto->combined_nonBDTDNE : sto->combined_BDTDNE;
  sto->transferEL->nedges = 0; // reset transferEL
  for(int i = 0; i < 2; i++) {
    if(sto->static_sto->maxl[i]) {
      EXEC_THROUGH_EDGES_EATH_NET_DECL(sto->static_sto->nodes[i], ego, alter, _tail, _head, edge, relevant_net, {
        if(!edgeflag || IN_DEG[alter] + OUT_DEG[alter] < sto->static_sto->bound) {
          UnsrtELToggleKnown(_tail, _head, sto->nonBDTDNE, edgeflag);
          UnsrtELToggleKnown(_tail, _head, sto->BDTDNE, !edgeflag);
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

  // update current dyad count
  sto->static_sto->currentdyads = sto->static_sto->proposeddyads;  
  
  // update the current submaximal edge count
  sto->static_sto->currentsubmaxledges = sto->static_sto->proposedsubmaxledges;    
}

MH_F_FN(Mf_discordBDTNT) {
  GET_STORAGE(discordBDTNTStorage, sto);
  // free things used only by the dynamic proposal
  UnsrtELDestroy(sto->BDTDNE);
  UnsrtELDestroy(sto->nonBDTDNE);
  UnsrtELDestroy(sto->discordantEdges);  
  NetworkDestroy(sto->combined_BDTDNE);
  NetworkDestroy(sto->combined_nonBDTDNE);
  UnsrtELDestroy(sto->transferEL);
  
  // let BDTNT's F_FN do most of the work
  MH_STORAGE = sto->static_sto;
  Mf_BDTNT(MHp, nwp);
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

  UnsrtEL **BDTDNE;
  UnsrtEL **nonBDTDNE;
  
  UnsrtEL **discordantEdges;
  UnsrtEL **nonDiscordantEdges;
  
  int in_discord;
    
  BDStratTNTStorage *static_sto;
} discordBDStratTNTStorage;

MH_I_FN(Mi_discordBDStratTNT) {
  // let BDStratTNT's I_FN do most of the work
  Mi_BDStratTNT(MHp, nwp);
  
  // copy a few things and handle purely temporal aspects
  BDStratTNTStorage *static_sto = MH_STORAGE;
  ALLOC_STORAGE(1, discordBDStratTNTStorage, sto);
  sto->nonDiscordantEdges = static_sto->els;
  sto->static_sto = static_sto;
  
  sto->BDTDNE = Calloc(sto->static_sto->nmixtypes, UnsrtEL *);
  sto->nonBDTDNE = Calloc(sto->static_sto->nmixtypes, UnsrtEL *);
  sto->discordantEdges = Calloc(sto->static_sto->nmixtypes, UnsrtEL *);  
  for(int i = 0; i < sto->static_sto->nmixtypes; i++) {
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

    for(int i = 0; i < sto->static_sto->nmixtypes; i++) {
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
  int strat_i = WtPopGetRand(sto->static_sto->wtp);
  
  // record the mixing type of the toggle, in case it's needed in the U function later
  sto->static_sto->stratmixingtype = strat_i;    

  int strattailtype = sto->static_sto->strattailtypes[strat_i];
  int stratheadtype = sto->static_sto->stratheadtypes[strat_i];
  int strat_diag = strattailtype == stratheadtype;
  
  // number of edges of this mixing type
  int nedgestype = sto->nonDiscordantEdges[strat_i]->nedges + sto->discordantEdges[strat_i]->nedges;
  
  Dyad ndyadstype = NodeListDyadCount(sto->static_sto->attrcounts[strattailtype], sto->static_sto->attrcounts[stratheadtype], sto->static_sto->bd_tails, sto->static_sto->bd_heads, sto->static_sto->bd_mixtypes[strat_diag], strat_diag, DIRECTED);
  
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
                           sto->static_sto->nodesvec[strattailtype], // tails
                           sto->static_sto->nodesvec[stratheadtype], // heads
                           sto->static_sto->bd_tails, // tailattrs
                           sto->static_sto->bd_heads, // headattrs
                           sto->static_sto->attrcounts[strattailtype], // tailcounts
                           sto->static_sto->attrcounts[stratheadtype], // headcounts
                           sto->static_sto->bd_mixtypes[strat_diag], // length
                           ndyadstype, // dyadcount
                           strat_diag, // diagonal
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

  sto->static_sto->strattailtype = sto->static_sto->strat_vattr[Mtail[0]];
  sto->static_sto->stratheadtype = sto->static_sto->strat_vattr[Mhead[0]];
    
  sto->static_sto->bdtailtype = sto->static_sto->bd_vattr[Mtail[0]];
  sto->static_sto->bdheadtype = sto->static_sto->bd_vattr[Mhead[0]];

  sto->static_sto->tailmaxl = IN_DEG[Mtail[0]] + OUT_DEG[Mtail[0]] == sto->static_sto->bound - 1 + in_network;
  sto->static_sto->headmaxl = IN_DEG[Mhead[0]] + OUT_DEG[Mhead[0]] == sto->static_sto->bound - 1 + in_network;

  sto->static_sto->nodes[0] = *Mtail;  
  sto->static_sto->nodes[1] = *Mhead;
        
  sto->static_sto->maxl[0] = sto->static_sto->tailmaxl;  
  sto->static_sto->maxl[1] = sto->static_sto->headmaxl;  

  // temporarily set tail and head toggleability to what it would be in the proposed network
  NodeListToggleKnownIf(*Mtail, sto->static_sto->nodesvec[sto->static_sto->strattailtype][sto->static_sto->bdtailtype], sto->static_sto->nodepos, sto->static_sto->attrcounts[sto->static_sto->strattailtype] + sto->static_sto->bdtailtype, !in_network, sto->static_sto->tailmaxl);
  NodeListToggleKnownIf(*Mhead, sto->static_sto->nodesvec[sto->static_sto->stratheadtype][sto->static_sto->bdheadtype], sto->static_sto->nodepos, sto->static_sto->attrcounts[sto->static_sto->stratheadtype] + sto->static_sto->bdheadtype, !in_network, sto->static_sto->headmaxl);

  // compute proposed dyad count for current mixing type (only)
  Dyad proposeddyadstype = NodeListDyadCount(sto->static_sto->attrcounts[strattailtype], sto->static_sto->attrcounts[stratheadtype], sto->static_sto->bd_tails, sto->static_sto->bd_heads, sto->static_sto->bd_mixtypes[strat_diag], strat_diag, DIRECTED);
    
  // here we compute the proposedcumprob, checking only those
  // mixing types that can be influenced by toggles made on 
  // the current mixing type
  sto->static_sto->proposedcumprob = sto->static_sto->currentcumprob;
  sto->static_sto->nmixtypestoupdate = 0; // reset counter
  // avoid these somewhat expensive checks in the typical case
  // where you have enough submaximal nodes that you cannot
  // be exhausting any mixing types of toggleable dyads
  if(sto->static_sto->attrcounts[sto->static_sto->strattailtype][sto->static_sto->bdtailtype] <= 2 || sto->static_sto->attrcounts[sto->static_sto->stratheadtype][sto->static_sto->bdheadtype] <= 2) {
    
    // how many strat types do we need to check?
    int ntocheck = strat_diag ? sto->static_sto->nstratlevels : 2*sto->static_sto->nstratlevels;

    for(int i = 0; i < ntocheck; i++) {
      // find the index of the i'th strat type we need to check, by looking it up in the indmat
      int infl_i = sto->static_sto->indmat[i < sto->static_sto->nstratlevels ? sto->static_sto->strattailtype : i - sto->static_sto->nstratlevels][i < sto->static_sto->nstratlevels ? i : sto->static_sto->stratheadtype];

      // if this strat type is not included in the proposal, or is the same as the strat type of the proposed toggle,
      // then it cannot change toggleability status, so skip it
      if(infl_i < 0 || infl_i == strat_i) {
        continue;
      }
      
      // can we toggle this mixing type in the current network?
      int toggle_curr = WtPopGetWt(infl_i, sto->static_sto->wtp) > 0;
      
      // will we be able to toggle this mixing type in the proposed network? 
      int toggle_prop = sto->nonDiscordantEdges[infl_i]->nedges > 0 || sto->discordantEdges[infl_i]->nedges > 0 || NodeListDyadCountPositive(sto->static_sto->attrcounts[sto->static_sto->strattailtypes[infl_i]], sto->static_sto->attrcounts[sto->static_sto->stratheadtypes[infl_i]], sto->static_sto->bd_tails, sto->static_sto->bd_heads, sto->static_sto->bd_mixtypes[sto->static_sto->strattailtypes[infl_i] == sto->static_sto->stratheadtypes[infl_i]], sto->static_sto->strattailtypes[infl_i] == sto->static_sto->stratheadtypes[infl_i]);
      
      // will there be a change in toggleability status?
      int change = toggle_curr - toggle_prop;

      // if so, take this into account      
      if(change) {
        sto->static_sto->proposedcumprob -= change*sto->static_sto->originalprobvec[infl_i];
        sto->static_sto->mixtypestoupdate[sto->static_sto->nmixtypestoupdate] = infl_i;
        sto->static_sto->nmixtypestoupdate++;        
      }
    }
  }
  
  // restore tail and head toggleability to their current status
  NodeListToggleKnownIf(*Mtail, sto->static_sto->nodesvec[sto->static_sto->strattailtype][sto->static_sto->bdtailtype], sto->static_sto->nodepos, sto->static_sto->attrcounts[sto->static_sto->strattailtype] + sto->static_sto->bdtailtype, in_network, sto->static_sto->tailmaxl);
  NodeListToggleKnownIf(*Mhead, sto->static_sto->nodesvec[sto->static_sto->stratheadtype][sto->static_sto->bdheadtype], sto->static_sto->nodepos, sto->static_sto->attrcounts[sto->static_sto->stratheadtype] + sto->static_sto->bdheadtype, in_network, sto->static_sto->headmaxl);


  // need to calculate number of submaxl edges and number of BD-toggleable discordant dyads
  // in the proposed network; this can involve both the proposal dyad itself and other dyads
  // containing either of the nodes in the proposal dyad
  int delta = in_network ? +1 : -1;

  int proposedsubmaxledgestype = sto->static_sto->currentsubmaxledgestype[strat_i];
  if(!sto->static_sto->tailmaxl && !sto->static_sto->headmaxl) {
    proposedsubmaxledgestype -= delta;
  }
    
  int propnddyadstype = nddyadstype;  
  propnddyadstype += in_discord ? -1 : 1;
  
  Network *relevant_net = in_network ? sto->combined_nonBDTDNE : sto->combined_BDTDNE;
  for(int i = 0; i < 2; i++) {
    if(sto->static_sto->maxl[i]) {
      EXEC_THROUGH_EDGES_EATH_NET_DECL(sto->static_sto->nodes[i], ego, alter, tail, head, edge, nwp, {
        if(alter != sto->static_sto->nodes[1 - i] && sto->static_sto->indmat[sto->static_sto->strat_vattr[tail]][sto->static_sto->strat_vattr[head]] == strat_i && IN_DEG[alter] + OUT_DEG[alter] < sto->static_sto->bound) {
          proposedsubmaxledgestype += delta;
        }
      });
      EXEC_THROUGH_EDGES_EATH_NET_DECL(sto->static_sto->nodes[i], ego, alter, tail, head, edge, relevant_net, {
        if(alter != sto->static_sto->nodes[1 - i] && sto->static_sto->indmat[sto->static_sto->strat_vattr[tail]][sto->static_sto->strat_vattr[head]] == strat_i && (!in_network || IN_DEG[alter] + OUT_DEG[alter] < sto->static_sto->bound)) {
          propnddyadstype += delta;
        }
      });
    }
  } 

  // calculate logratio
  double prob_weight = sto->static_sto->currentcumprob/sto->static_sto->proposedcumprob;
  
  double forward_network = in_network ? (ndyadstype == 0 ? 1.0/nedgestype : 0.5/nedgestype + (0.5/ndyadstype)*((double)sto->static_sto->currentsubmaxledgestype[strat_i]/nedgestype)) : (nedgestype == 0 ? 1.0/ndyadstype : 0.5/ndyadstype);
  
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
  if(sto->static_sto->nmixtypestoupdate > 0) {
    sto->static_sto->currentcumprob = sto->static_sto->proposedcumprob;

    if(edgeflag) {
      sto->static_sto->nmixtypes_toggleable += sto->static_sto->nmixtypestoupdate;
      for(int i = 0; i < sto->static_sto->nmixtypestoupdate; i++) {
        WtPopSetWt(sto->static_sto->mixtypestoupdate[i], sto->static_sto->originalprobvec[sto->static_sto->mixtypestoupdate[i]], sto->static_sto->wtp);          
      }
    } else {
      sto->static_sto->nmixtypes_toggleable -= sto->static_sto->nmixtypestoupdate;
      for(int i = 0; i < sto->static_sto->nmixtypestoupdate; i++) {
        WtPopSetWt(sto->static_sto->mixtypestoupdate[i], 0, sto->static_sto->wtp);          
      }
    }
  }

  // add or remove the dyad being toggled from the relevant edge set(s)/network
  if(sto->in_discord == edgeflag) {
    UnsrtELToggleKnown(tail, head, sto->discordantEdges[sto->static_sto->stratmixingtype], edgeflag);      
  } else {
    UnsrtELToggleKnown(tail, head, sto->nonDiscordantEdges[sto->static_sto->stratmixingtype], edgeflag);
    UnsrtELToggleKnown(tail, head, sto->BDTDNE[sto->static_sto->stratmixingtype], !edgeflag);
    ToggleKnownEdge(tail, head, sto->combined_BDTDNE, !edgeflag);      
  }

  NodeListToggleKnownIf(tail, sto->static_sto->nodesvec[sto->static_sto->strattailtype][sto->static_sto->bdtailtype], sto->static_sto->nodepos, sto->static_sto->attrcounts[sto->static_sto->strattailtype] + sto->static_sto->bdtailtype, !edgeflag, sto->static_sto->tailmaxl);
  NodeListToggleKnownIf(head, sto->static_sto->nodesvec[sto->static_sto->stratheadtype][sto->static_sto->bdheadtype], sto->static_sto->nodepos, sto->static_sto->attrcounts[sto->static_sto->stratheadtype] + sto->static_sto->bdheadtype, !edgeflag, sto->static_sto->headmaxl);      
  
  // update submaxl edge counts and dyad toggleability statuses, as appropriate
  int delta = edgeflag ? +1 : -1;
  
  if(!sto->static_sto->tailmaxl && !sto->static_sto->headmaxl) {
    sto->static_sto->currentsubmaxledgestype[sto->static_sto->stratmixingtype] -= delta;
  }
  
  Network *relevant_net = edgeflag ? sto->combined_nonBDTDNE : sto->combined_BDTDNE;
  sto->transferEL->nedges = 0; // reset transferEL
  for(int i = 0; i < 2; i++) {
    if(sto->static_sto->maxl[i]) {
      EXEC_THROUGH_EDGES_EATH_NET_DECL(sto->static_sto->nodes[i], ego, alter, _tail, _head, edge, nwp, {
        if(alter != sto->static_sto->nodes[1 - i] && IN_DEG[alter] + OUT_DEG[alter] < sto->static_sto->bound) {
          int stratmixingtype = sto->static_sto->indmat[sto->static_sto->strat_vattr[_tail]][sto->static_sto->strat_vattr[_head]];
          if(stratmixingtype >= 0) {
            sto->static_sto->currentsubmaxledgestype[stratmixingtype] += delta;
          }
        }
      });
      
      EXEC_THROUGH_EDGES_EATH_NET_DECL(sto->static_sto->nodes[i], ego, alter, _tail, _head, edge, relevant_net, {
        if(!edgeflag || IN_DEG[alter] + OUT_DEG[alter] < sto->static_sto->bound) {
          int stratmixingtype = sto->static_sto->indmat[sto->static_sto->strat_vattr[_tail]][sto->static_sto->strat_vattr[_head]];
          UnsrtELToggleKnown(_tail, _head, sto->nonBDTDNE[stratmixingtype], edgeflag);
          UnsrtELToggleKnown(_tail, _head, sto->BDTDNE[stratmixingtype], !edgeflag);
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
  
  // let BDStratTNT's F_FN do most of the work
  MH_STORAGE = sto->static_sto;
  Mf_BDStratTNT(MHp, nwp);
  Free(sto->static_sto);
  MH_STORAGE = sto;
  // MHp->storage itself should be Freed by MHProposalDestroy
}
