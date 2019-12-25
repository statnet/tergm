/*  File src/changestats_test.c in package ergm, part of the Statnet suite
 *  of packages for network analysis, http://statnet.org .
 *
 *  This software is distributed under the GPL-3 license.  It is free,
 *  open source, and has the attribution requirements (GPL Section 7) at
 *  http://statnet.org/attribution
 *
 *  Copyright 2003-2017 Statnet Commons
 */
#include "ergm_changestat_auxnet.h"
#include "tergm_changestats_auxnet.h"
#include "ergm_dyad_hashmap.h"
#include "ergm_dyad_hashmap_utils.h"
#include "tergm_model.h"

// Initialize empty aux network. Then loop through the edges of y0 (ref_el).
// If the edge also exists in y1, then keep it in auxnet->onwp.
// The storage auxnet->onwp should be initialized as y0&y1 at the end.
I_CHANGESTAT_FN(i__intersect_lt_net_Network){
  I_AUXNET(NetworkInitialize(NULL, NULL, 0, N_NODES, DIRECTED, BIPARTITE, FALSE, 0, NULL));
  GET_AUX_STORAGE_NUM(StoreTimeAndLasttoggle, dur_inf, 1);
  TailHead dyad;
  int lt;
  const int t = dur_inf->time;
  kh_foreach(dur_inf->lasttoggle, dyad, lt, {
      // Time difference of 0 means that the edge was absent either in this time step or in the previous time step.
      if(t-lt!=0 && IS_OUTEDGE(dyad.tail, dyad.head)){
        ToggleKnownEdge(dyad.tail,dyad.head, auxnet->onwp, FALSE);
      }
    });
}
U_CHANGESTAT_FN(u__intersect_lt_net_Network){
  GET_AUX_STORAGE(StoreAuxnet, auxnet);
  GET_AUX_STORAGE_NUM(StoreTimeAndLasttoggle, dur_inf, 1);
  // If the edge is not in y0, changing y1 won't matter. We infer that the
  // edge was in y0 if either it's in y1 *and* last toggle time is
  // not current time *or* it's not in y1 *and* last toggle time is
  // current time.  Note that if the key is not found in lasttoggle,
  // then it will return t+1 and therefore != t.
  const int t = dur_inf->time;
  if(edgeflag != (kh_getval(DyadMapInt, dur_inf->lasttoggle, THKey(dur_inf->lasttoggle,tail,head), t+1)==t))
    ToggleKnownEdge(tail, head, auxnet->onwp, edgeflag);
}

X_CHANGESTAT_FN(x__intersect_lt_net_Network){
  switch(type){
  case TOCK:
    {
      GET_AUX_STORAGE(StoreAuxnet, auxnet);
      GET_AUX_STORAGE_NUM(StoreTimeAndLasttoggle, dur_inf, 1);
      TailHead dyad;
      kh_foreach_key(dur_inf->discord, dyad, {
          if(IS_OUTEDGE(dyad.tail, dyad.head));
          ToggleKnownEdge(dyad.tail, dyad.head, auxnet->onwp, FALSE);
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

// Initialize aux network to y1. Then loop through the edges of y0 (ref_el).
// If the edge does not exists in y1, then toggle it on in aux network.
// The storage auxnet->onwp should be initialized as y0|y1 at the end.
I_CHANGESTAT_FN(i__union_lt_net_Network){
  I_AUXNET(NetworkCopy(nwp));
  GET_AUX_STORAGE_NUM(StoreTimeAndLasttoggle, dur_inf, 1);
  TailHead dyad;
  int lt;
  const int t = dur_inf->time;
  kh_foreach(dur_inf->lasttoggle, dyad, lt, {
      // Time difference of 0 means that edge was present either in this time step or in the previous time step.
      if(t-lt==0 || !IS_OUTEDGE(dyad.tail, dyad.head)){
        ToggleKnownEdge(dyad.tail,dyad.head, auxnet->onwp, FALSE);
      }
    });
}

U_CHANGESTAT_FN(u__union_lt_net_Network){
  GET_AUX_STORAGE(StoreAuxnet, auxnet);
  GET_AUX_STORAGE_NUM(StoreTimeAndLasttoggle, dur_inf, 1);
  // If the edge is in y0, changing y1 won't matter. We infer that the
  // edge was not in y0 if either it's in y1 *and* last toggle time is
  // current time *or* it's not in y1 *and* last toggle time is not
  // current time.  Note that if the key is not found in lasttoggle,
  // then it will return t+1 and therefore != t.
  const int t = dur_inf->time;
  if(edgeflag == (kh_getval(DyadMapInt, dur_inf->lasttoggle, THKey(dur_inf->lasttoggle,tail,head), t+1)==t))
    ToggleKnownEdge(tail, head, auxnet->onwp, edgeflag);
}

X_CHANGESTAT_FN(x__union_lt_net_Network){
  switch(type){
  case TOCK:
    {
      GET_AUX_STORAGE(StoreAuxnet, auxnet);
      GET_AUX_STORAGE_NUM(StoreTimeAndLasttoggle, dur_inf, 1);
      TailHead dyad;
      kh_foreach_key(dur_inf->discord, dyad, {
          if(!IS_OUTEDGE(dyad.tail, dyad.head));
          ToggleKnownEdge(dyad.tail, dyad.head, auxnet->onwp, TRUE);
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
