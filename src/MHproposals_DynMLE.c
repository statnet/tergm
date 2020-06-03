#include "ergm_MHproposal.h"
#include "ergm_edgelist.h"
#include "ergm_changestat.h"
#include "ergm_rlebdm.h"
#include "ergm_MHstorage.h"
#include "ergm_unsorted_edgelist.h"
#include "changestats_lasttoggle.h"
#include "tergm_model.h"

typedef struct {
  UnsrtEL *discordantEdges;
  UnsrtEL *discordantNonEdges;

  int in_discord;

} staticDiscordTNTStorage;

MH_I_FN(Mi_staticDiscordTNT){
  MHp->ntoggles = 1;

  ALLOC_STORAGE(1, staticDiscordTNTStorage, sto);
  sto->discordantEdges = UnsrtELInitialize(0, NULL, NULL, FALSE);
  sto->discordantNonEdges = UnsrtELInitialize(0, NULL, NULL, FALSE);
  EXEC_THROUGH_NET_EDGES(t, h, e, {
      if(iEdgeListSearch(t, h, MH_IINPUTS)==0)
        UnsrtELInsert(t, h, sto->discordantEdges);
    });

  for(Edge i=0; i<EDGECOUNT(nwp); i++){
    Vertex t=MH_IINPUTS[1+i], h=MH_IINPUTS[1+EDGECOUNT(nwp)+i];
    if(IS_OUTEDGE(t, h)){
      UnsrtELInsert(t, h, sto->discordantNonEdges);      
    }
  }
}

MH_P_FN(MH_staticDiscordTNT) {
  GET_STORAGE(staticDiscordTNTStorage, sto);

  int in_discord;
  int in_network;

  int nedges = EDGECOUNT(nwp);
  int nddyads = sto->discordantEdges->nedges+sto->discordantNonEdges->nedges;

  if(nddyads == 0 || unif_rand() < 0.5) {
    // propose from network
    if(nedges == 0 || unif_rand() < 0.5) {
      // propose toggling a random dyad
      GetRandDyad(Mtail, Mhead, nwp);
      in_network = IS_OUTEDGE(Mtail[0], Mhead[0]);
      in_discord = XOR(in_network, iEdgeListSearch(Mtail[0], Mhead[0], MH_IINPUTS));
    } else {
      // propose toggling off an edge in network
      GetRandEdge(Mtail, Mhead, nwp);
      in_network = TRUE;
      in_discord = XOR(in_network, iEdgeListSearch(Mtail[0], Mhead[0], MH_IINPUTS));
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

  sto->in_discord = in_discord;
}

MH_U_FN(Mu_staticDiscordTNT) {
  // add or remove the dyad from the appropriate discordance edgelist

  GET_STORAGE(staticDiscordTNTStorage, sto);

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

MH_F_FN(Mf_staticDiscordTNT) {
  GET_STORAGE(staticDiscordTNTStorage, sto);

  UnsrtELDestroy(sto->discordantNonEdges);
  UnsrtELDestroy(sto->discordantEdges);
}
