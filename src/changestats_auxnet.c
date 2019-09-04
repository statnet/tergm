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
#include "ergm_changestat_operator.h"

// Initialize empty aux network. Then loop through the edges of y0 (ref_el).
// If the edge also exists in y1, then keep it in auxnet->onwp.
// The storage auxnet->onwp should be initialized as y0&y1 at the end.
I_CHANGESTAT_FN(i__intersect_lt_net_Network){
  I_AUXNET(NetworkInitialize(NULL, NULL, 0, N_NODES, DIRECTED, BIPARTITE, FALSE, 0, NULL));
  TailHead dyad;
  kh_foreach_key(nwp->duration_info->lasttoggle, dyad, {
      if(IS_OUTEDGE(dyad.tail, dyad.head)){
        ToggleKnownEdge(dyad.tail,dyad.head, auxnet->onwp, FALSE);
      }
    });
  // Steal the duration_info data structure from input network.
  auxnet->onwp->duration_info = auxnet->inwp->duration_info;
}

U_CHANGESTAT_FN(u__intersect_lt_net_Network){
  GET_AUX_STORAGE(StoreAuxnet, auxnet);
  // only toggle if the edge is in y0. otherwise changing y1 won't matter.
  if(HASDMI(tail, head, auxnet->inwp->duration_info->lasttoggle))
    ToggleEdge(tail, head, auxnet->onwp);
}

F_CHANGESTAT_FN(f__intersect_lt_net_Network){
  GET_AUX_STORAGE(StoreAuxnet, auxnet);
  auxnet->onwp->duration_info = NULL; // So that we don't deallocate someone else's structure.
  NetworkDestroy(auxnet->onwp);
}

// Initialize aux network to y1. Then loop through the edges of y0 (ref_el).
// If the edge does not exists in y1, then toggle it on in aux network.
// The storage auxnet->onwp should be initialized as y0|y1 at the end.
I_CHANGESTAT_FN(i__union_lt_net_Network){
  I_AUXNET(NetworkCopy(nwp));
  TailHead dyad;
  kh_foreach_key(nwp->duration_info->lasttoggle, dyad, {
      if(!IS_OUTEDGE(dyad.tail, dyad.head)){
        ToggleKnownEdge(dyad.tail,dyad.head, auxnet->onwp, FALSE);
      }
    });
  // Steal the duration_info data structure from input network.
  auxnet->onwp->duration_info = auxnet->inwp->duration_info;
}

U_CHANGESTAT_FN(u__union_lt_net_Network){
  GET_AUX_STORAGE(StoreAuxnet, auxnet);
  // If the edge is in y0, changing y1 won't matter.
  if(!HASDMI(tail, head, auxnet->inwp->duration_info->lasttoggle))
    ToggleEdge(tail, head, auxnet->onwp);
}

F_CHANGESTAT_FN(f__union_lt_net_Network){
  GET_AUX_STORAGE(StoreAuxnet, auxnet);
  auxnet->onwp->duration_info = NULL; // So that we don't deallocate someone else's structure.
  NetworkDestroy(auxnet->onwp);
}
