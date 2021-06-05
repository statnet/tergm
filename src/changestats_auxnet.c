/*  File src/changestats_auxnet.c in package tergm, part of the Statnet suite
 *  of packages for network analysis, https://statnet.org .
 *
 *  This software is distributed under the GPL-3 license.  It is free,
 *  open source, and has the attribution requirements (GPL Section 7) at
 *  https://statnet.org/attribution
 *
 *  Copyright 2008-2020 Statnet Commons
 */
#include "ergm_changestat_auxnet.h"
#include "tergm_changestats_auxnet.h"
#include "ergm_dyad_hashmap.h"
#include "ergm_dyad_hashmap_utils.h"
#include "tergm_model.h"

// Initialize aux network to an empty network. Then, loop through the changes
// from last time step and any changed edges. The storage
// auxnet->onwp should be initialized as y0 xor y1 at the end.
I_CHANGESTAT_FN(i__discord_lt_net_Network){
  I_AUXNET(NetworkInitialize(NULL, NULL, 0, N_NODES, DIRECTED, BIPARTITE, FALSE, 0, NULL));
  GET_AUX_STORAGE_NUM(StoreTimeAndLasttoggle, dur_inf, 1);
  TailHead dyad;
  kh_foreach_key(dur_inf->discord, dyad, {
      AddEdgeToTrees(dyad.tail,dyad.head, auxnet->onwp);
    });
}

U_CHANGESTAT_FN(u__discord_lt_net_Network){
  GET_AUX_STORAGE(StoreAuxnet, auxnet);
  GET_AUX_STORAGE_NUM(StoreTimeAndLasttoggle, dur_inf, 1);
  if(dur_inf->ticktock){
      ToggleKnownEdge(tail, head, auxnet->onwp, edgestate);
  } // In fiat mode, discord never changes.
}

X_CHANGESTAT_FN(x__discord_lt_net_Network){
  switch(type){
  case TICK:
    {
      GET_AUX_STORAGE(StoreAuxnet, auxnet);
      // There are two ways to do this. One is to drop the old network
      // and make a copy of the new one. The other is to make
      // incremental changes. TODO: Benchmark.

      /* NetworkDestroy(auxnet->onwp); */
      /* auxnet->onwp = NetworkCopy(nwp); */

      GET_AUX_STORAGE_NUM(StoreTimeAndLasttoggle, dur_inf, 1);
      TailHead dyad;
      // Here, we need to toggle all edges that were toggled in the current
      // (soon to be previous) time step.
      kh_foreach_key(dur_inf->discord, dyad, {
          ToggleEdge(dyad.tail, dyad.head, auxnet->onwp);
        });
    }
    break;
  default: break;
  }
}

F_CHANGESTAT_FN(f__discord_lt_net_Network){
  GET_AUX_STORAGE(StoreAuxnet, auxnet);
  NetworkDestroy(auxnet->onwp);
}


// Initialize aux network equal to y0. Then, loop through the changes
// from last time step and toggle any changed edges off. The storage
// auxnet->onwp should be initialized as y0&y1 at the end.
I_CHANGESTAT_FN(i__intersect_lt_net_Network){
  I_AUXNET(NetworkCopy(nwp));
  GET_AUX_STORAGE_NUM(StoreTimeAndLasttoggle, dur_inf, 1);
  TailHead dyad;
  kh_foreach_key(dur_inf->discord, dyad, {
      DeleteEdgeFromTrees(dyad.tail,dyad.head, auxnet->onwp);
    });
}

U_CHANGESTAT_FN(u__intersect_lt_net_Network){
  GET_AUX_STORAGE(StoreAuxnet, auxnet);
  GET_AUX_STORAGE_NUM(StoreTimeAndLasttoggle, dur_inf, 1);
  if(dur_inf->ticktock){
    // If the edge is not in y0, changing y1 won't matter. We infer that the
    // edge was in y0 if either it's in y1 *and* last toggle time is
    // not current time *or* it's not in y1 *and* last toggle time is
    // current time.  Note that if the key is not found in lasttoggle,
    // then it will return t+1 and therefore != t.
    if(edgestate!=JUST_CHANGED(dur_inf,tail,head))
      ToggleKnownEdge(tail, head, auxnet->onwp, edgestate);
  }else{
    // If we are in the "fiat" mode, then if the dyad just changed,
    // then it can only flip between (0,1) and (1,0) and so stays at
    // 0, but if not, then it flips between (0,0) and (1,1).
    if(!JUST_CHANGED(dur_inf,tail,head))
      ToggleKnownEdge(tail, head, auxnet->onwp, edgestate);
  }
}

X_CHANGESTAT_FN(x__intersect_lt_net_Network){
  switch(type){
  case TICK:
    {
      GET_AUX_STORAGE(StoreAuxnet, auxnet);
      // There are two ways to do this. One is to drop the old network
      // and make a copy of the new one. The other is to make
      // incremental changes. TODO: Benchmark.

      /* NetworkDestroy(auxnet->onwp); */
      /* auxnet->onwp = NetworkCopy(nwp); */

      GET_AUX_STORAGE_NUM(StoreTimeAndLasttoggle, dur_inf, 1);
      TailHead dyad;
      // Here, we need to add all edges that were added in the current
      // (soon to be previous) time step. One way to discover them is
      // to go through all the recently toggled edges and test if they
      // are present in the current network.
      kh_foreach_key(dur_inf->discord, dyad, {
        if(IS_OUTEDGE(dyad.tail, dyad.head))
          AddEdgeToTrees(dyad.tail, dyad.head, auxnet->onwp);
        });
    }
    break;
  default: break;
  }
}

F_CHANGESTAT_FN(f__intersect_lt_net_Network){
  GET_AUX_STORAGE(StoreAuxnet, auxnet);
  NetworkDestroy(auxnet->onwp);
}

// Initialize aux network to y1. Then loop through the changed edges
// and add any that are absent. The storage auxnet->onwp should be
// initialized as y0|y1 at the end.
I_CHANGESTAT_FN(i__union_lt_net_Network){
  I_AUXNET(NetworkCopy(nwp));
  GET_AUX_STORAGE_NUM(StoreTimeAndLasttoggle, dur_inf, 1);
  TailHead dyad;
  kh_foreach_key(dur_inf->discord, dyad, {
    if(EdgetreeSearch(dyad.tail, dyad.head, auxnet->onwp->outedges) == 0) {
      AddEdgeToTrees(dyad.tail,dyad.head, auxnet->onwp);
    }
  });
}

U_CHANGESTAT_FN(u__union_lt_net_Network){
  GET_AUX_STORAGE(StoreAuxnet, auxnet);
  GET_AUX_STORAGE_NUM(StoreTimeAndLasttoggle, dur_inf, 1);
  if(dur_inf->ticktock){
    // If the edge is in y0, changing y1 won't matter. We infer that the
    // edge was not in y0 if either it's in y1 *and* last toggle time is
    // current time *or* it's not in y1 *and* last toggle time is not
    // current time.  Note that if the key is not found in lasttoggle,
    // then it will return t+1 and therefore != t.
    if(edgestate == JUST_CHANGED(dur_inf,tail,head))
      ToggleKnownEdge(tail, head, auxnet->onwp, edgestate);
  }else{
    // If we are in the "fiat" mode, then if the dyad just changed,
    // then it can only flip between (0,1) and (1,0) and so stays at
    // 1, but if not, then it flips between (0,0) and (1,1).
    if(!JUST_CHANGED(dur_inf,tail,head))
      ToggleKnownEdge(tail, head, auxnet->onwp, edgestate);
  }
}

X_CHANGESTAT_FN(x__union_lt_net_Network){
  switch(type){
  case TICK:
    {
      GET_AUX_STORAGE(StoreAuxnet, auxnet);
      // There are two ways to do this. One is to drop the old network
      // and make a copy of the new one. The other is to make
      // incremental changes. TODO: Benchmark.

      /* NetworkDestroy(auxnet->onwp); */
      /* auxnet->onwp = NetworkCopy(nwp); */

      GET_AUX_STORAGE_NUM(StoreTimeAndLasttoggle, dur_inf, 1);
      TailHead dyad;
      // Here, we need to delete all edges that were toggled off in
      // the current (soon to be previous) time step. One way to
      // discover them is to go through all the recently toggled edges
      // and test if they are absent in the current network.
      kh_foreach_key(dur_inf->discord, dyad, {
          if(!IS_OUTEDGE(dyad.tail, dyad.head))
            DeleteEdgeFromTrees(dyad.tail, dyad.head, auxnet->onwp);
        });
    }
    break;
  default: break;
  }
}

F_CHANGESTAT_FN(f__union_lt_net_Network){
  GET_AUX_STORAGE(StoreAuxnet, auxnet);
  NetworkDestroy(auxnet->onwp);
}

// Initialize aux network to y1. Then loop through the last toggle
// structure. If it just toggled, then reverse it. Then storage should
// be initialized as y0 at the end.
I_CHANGESTAT_FN(i__previous_lt_net_Network){
  I_AUXNET(NetworkCopy(nwp));
  GET_AUX_STORAGE_NUM(StoreTimeAndLasttoggle, dur_inf, 1);
  TailHead dyad;
  kh_foreach_key(dur_inf->discord, dyad, {
      ToggleEdge(dyad.tail,dyad.head, auxnet->onwp);
    });
}

U_CHANGESTAT_FN(u__previous_lt_net_Network){
  GET_AUX_STORAGE(StoreAuxnet, auxnet);
  GET_AUX_STORAGE_NUM(StoreTimeAndLasttoggle, dur_inf, 1);

  // If we are within a time step (between a TICK and a TOCK), then
  // the previous network state is fixed. Otherwise, toggles of this
  // network tell us what the previous time step *must have been*.
  if(!dur_inf->ticktock) ToggleEdge(tail, head, auxnet->onwp);
}

X_CHANGESTAT_FN(x__previous_lt_net_Network){
  switch(type){
  case TICK: // We are about to increment the clock. This means that the current network is about to become the previous network.
    {
      GET_AUX_STORAGE(StoreAuxnet, auxnet);
      // There are two ways to do this. One is to drop the old network
      // and make a copy of the new one. The other is to make
      // incremental changes. TODO: Benchmark.

      /* NetworkDestroy(auxnet->onwp); */
      /* auxnet->onwp = NetworkCopy(nwp); */

      GET_AUX_STORAGE_NUM(StoreTimeAndLasttoggle, dur_inf, 1);
      TailHead dyad;
      // I.e., we toggle what was just toggled.
      kh_foreach_key(dur_inf->discord, dyad, {
            ToggleEdge(dyad.tail, dyad.head, auxnet->onwp);
        });
    }
    break;
  default: break;
  }
}

F_CHANGESTAT_FN(f__previous_lt_net_Network){
  GET_AUX_STORAGE(StoreAuxnet, auxnet);
  NetworkDestroy(auxnet->onwp);
}
