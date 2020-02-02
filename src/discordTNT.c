#include "ergm_MHproposal.h"
#include "ergm_edgelist.h"
#include "ergm_changestat.h"
#include "ergm_rlebdm.h"
#include "ergm_MHstorage.h"
#include "ergm_unsorted_edgelist.h"
#include "changestats_lasttoggle.h"
#include "tergm_model.h"

MH_I_FN(Mi_discordTNT) {
  MHp->ntoggles = 1;
  
  ALLOC_STORAGE(3, void *, sto);
  sto[0] = UnsrtELInitialize(0, NULL, NULL, FALSE);
  sto[1] = UnsrtELInitialize(0, NULL, NULL, FALSE);
  sto[2] = Calloc(1, int);
  
  // we shouldn't actually use these for anything, but
  // for consistency I'm initializing the discordant edgelists
  // to what's in dur_inf->discord
  // it would probably be faster for some purposes to skip this though
  GET_AUX_STORAGE(StoreTimeAndLasttoggle, dur_inf);
  StoreDyadMapInt *discord = dur_inf->discord;
  TailHead dyad;
  kh_foreach_key(discord, dyad, {
    if(IS_OUTEDGE(dyad.tail, dyad.head)) {
      UnsrtELInsert(dyad.tail, dyad.head, sto[0]);
    } else {
      UnsrtELInsert(dyad.tail, dyad.head, sto[1]);
    }
  });
}

MH_X_FN(Mx_discordTNT) {
  GET_STORAGE(void *, sto);
  
  if(type == TICK) {
    // "clear" the discordant dyads
    ((UnsrtEL *)sto[0])->nedges = 0;
    ((UnsrtEL *)sto[1])->nedges = 0;
  }
}

MH_P_FN(MH_discordTNT) {
  GET_AUX_STORAGE(StoreTimeAndLasttoggle, dur_inf);
  
  GET_STORAGE(void *, sto);
  UnsrtEL *discordantEdges = sto[0];
  UnsrtEL *discordantNonEdges = sto[1];
  
  int in_discord;
  int in_network;
  
  int nedges = EDGECOUNT(nwp);
  int nddyads = kh_size(dur_inf->discord);
  
  if(nddyads == 0 || unif_rand() < 0.5) {
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
        UnsrtELGetRand(Mtail, Mhead, discordantEdges);
      } else {
        UnsrtELGetRand(Mtail, Mhead, discordantNonEdges);
      }
    }
  } else {
    // propose from discord
    if(unif_rand() < discordantEdges->nedges/((double) nddyads)) {
      UnsrtELGetRand(Mtail, Mhead, discordantEdges);
      in_network = TRUE;
    } else {
      UnsrtELGetRand(Mtail, Mhead, discordantNonEdges);
      in_network = FALSE;
    }
    in_discord = TRUE;
  }
  
  // compute logratio
  
  int ndyads = DYADCOUNT(nwp);
  
  // these ignore overall factor of 1/2 which cancels out in ratio
  double forward_discord = in_discord ? 1.0/nddyads : 0;
  double backward_discord = in_discord ? 0 : 1.0/(1 + nddyads);
  
  double forward_network = in_network ? (0.5/nedges + 0.5/ndyads) : (nedges == 0 ? 1.0/ndyads : 0.5/ndyads);
  double backward_network = in_network ? (nedges == 1 ? 1.0/ndyads : 0.5/ndyads) : (0.5/(nedges + 1) + 0.5/ndyads);
  
  if(nddyads == 0) forward_network *= 2;
  if(nddyads == 1 && in_discord) backward_network *= 2;
  
  double forward = forward_discord + forward_network;
  double backward = backward_discord + backward_network;

  MHp->logratio = log(backward/forward);

  ((int *)sto[2])[0] = in_discord;
}

MH_U_FN(Mu_discordTNT) {
  // add or remove the dyad from the appropriate discordance edgelist
  
  GET_STORAGE(void *, sto);
  
  if(((int *)sto[2])[0]) {
    // currently in discord; take it out
    if(edgeflag) {        
      UnsrtELDelete(tail, head, sto[0]);
    } else {
      UnsrtELDelete(tail, head, sto[1]);  
    }
  } else {
    // not currently in discord; add it in
    if(edgeflag) {
      UnsrtELInsert(tail, head, sto[1]);
    } else {
      UnsrtELInsert(tail, head, sto[0]);        
    }
  }
}

MH_F_FN(Mf_discordTNT) {
  GET_STORAGE(void *, sto);
  
  Free(sto[2]);
  UnsrtELDestroy(sto[1]);
  UnsrtELDestroy(sto[0]);
}
