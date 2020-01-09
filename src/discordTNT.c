#include "ergm_MHproposal.h"
#include "ergm_edgelist.h"
#include "ergm_changestat.h"
#include "ergm_rlebdm.h"
#include "ergm_MHstorage.h"
#include "ergm_unsorted_edgelist.h"
#include "changestats_lasttoggle.h"

MH_I_FN(Mi_discordTNT) {
  MHp->ntoggles = 1;
}

MH_P_FN(MH_discordTNT) {
  GET_AUX_STORAGE(StoreTimeAndLasttoggle, dur_inf);
  
  int nedges = EDGECOUNT(nwp);
  int ndedges = kh_size(dur_inf->discord);
  
  if(unif_rand() < 0.5) {
    // propose from discord
    if(ndedges == 0 || unif_rand() < 0.5) {
      // propose toggling a random dyad
      GetRandDyad(Mtail, Mhead, nwp);
    } else {
      // propose toggling off an edge in discord
      int index = unif_rand()*ndedges;
      
      khint_t iter = kh_begin(dur_inf->discord);
      while(true) {
        if(kh_exist(dur_inf->discord, iter)) {
          if(index == 0) {
            break;
          } else {
            index--;
          }
        }
        iter++;
      }
      
      TailHead edge = kh_key(dur_inf->discord, iter);
      Mtail[0] = edge.tail;
      Mhead[0] = edge.head;
    }
  } else {
    // propose from network
    if(nedges == 0 || unif_rand() < 0.5) {
      // propose toggling a random dyad
      GetRandDyad(Mtail, Mhead, nwp);
    } else {
      // propose toggling off an edge in network
      GetRandEdge(Mtail, Mhead, nwp);
    }
  }
  
  // compute logratio
  int in_network = IS_OUTEDGE(Mtail[0], Mhead[0]);
  int in_discord = kh_get(DyadMapInt, dur_inf->discord, THKey(dur_inf->discord, Mtail[0], Mhead[0])) != kh_none;
  
  int ndyads = DYADCOUNT(nwp);
  
  // these ignore overall factor of 1/2 which cancels out in ratio
  double forward_discord = in_discord ? (0.5/ndedges + 0.5/ndyads) : (ndedges == 0 ? 1.0/ndyads : 0.5/ndyads);
  double backward_discord = in_discord ? (ndedges == 1 ? 1.0/ndyads : 0.5/ndyads) : (0.5/(ndedges + 1) + 0.5/ndyads);
  
  double forward_network = in_network ? (0.5/nedges + 0.5/ndyads) : (nedges == 0 ? 1.0/ndyads : 0.5/ndyads);
  double backward_network = in_network ? (nedges == 1 ? 1.0/ndyads : 0.5/(ndyads)) : (0.5/(nedges + 1) + 0.5/ndyads);
  
  double forward = forward_discord + forward_network;
  double backward = backward_discord + backward_network;
  
/*  Rprintf("tail head: %d %d\n", Mtail[0], Mhead[0]);
  Rprintf("logratio: %f\n", log(backward/forward));
  
  Rprintf("kh_size: %d\n", ndedges);
  Rprintf("network edgecount: %d\n", nedges);
  Rprintf("network dyadcount: %d\n", ndyads);
  
  Rprintf("forward: %f\n", forward);
  Rprintf("backward: %f\n", backward);
*/  
  MHp->logratio = log(backward/forward);
}
