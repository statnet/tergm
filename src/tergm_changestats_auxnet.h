/*  File src/tergm_changestats_auxnet.h in package tergm, part of the
 *  Statnet suite of packages for network analysis, https://statnet.org .
 *
 *  This software is distributed under the GPL-3 license.  It is free,
 *  open source, and has the attribution requirements (GPL Section 7) at
 *  https://statnet.org/attribution .
 *
 *  Copyright 2008-2021 Statnet Commons
 */
#ifndef _TERGM_CHANGESTATS_AUXNET_H_
#define _TERGM_CHANGESTATS_AUXNET_H_

#include "changestats_lasttoggle.h"
#include "ergm_changestat_auxnet.h"

#define map_toggle_maxtoggles__discord_lt_net_Network 1
MAP_TOGGLE_FN(map_toggle__discord_lt_net_Network){
  ModelTerm *mtp = auxnet->mtp;
  GET_AUX_STORAGE_NUM(StoreTimeAndLasttoggle, dur_inf, 1);
  // See comments for u__discord_lt_net_Network() for explanation.
  if(dur_inf->ticktock){
      MAP_TOGGLE_PROPAGATE;
  }else return 0;

}

#define map_toggle_maxtoggles__intersect_lt_net_Network 1
MAP_TOGGLE_FN(map_toggle__intersect_lt_net_Network){
  ModelTerm *mtp = auxnet->mtp;
  GET_AUX_STORAGE_NUM(StoreTimeAndLasttoggle, dur_inf, 1);
  // See comments for u__intersect_lt_net_Network() for explanation.
  if(dur_inf->ticktock){
      MAP_TOGGLE_PROPAGATE_IF(edgestate != JUST_CHANGED(dur_inf,tail,head));
  }else{
      MAP_TOGGLE_PROPAGATE_IF(!JUST_CHANGED(dur_inf,tail,head));
  }

}

#define map_toggle_maxtoggles__union_lt_net_Network 1
MAP_TOGGLE_FN(map_toggle__union_lt_net_Network){
  ModelTerm *mtp = auxnet->mtp;
  GET_AUX_STORAGE_NUM(StoreTimeAndLasttoggle, dur_inf, 1);
  // See comments for u__union_lt_net_Network() for explanation.
  if(dur_inf->ticktock){
      MAP_TOGGLE_PROPAGATE_IF(edgestate == JUST_CHANGED(dur_inf,tail,head));
  }else{
      MAP_TOGGLE_PROPAGATE_IF(!JUST_CHANGED(dur_inf,tail,head));
  }
}

#define map_toggle_maxtoggles__previous_lt_net_Network 1
MAP_TOGGLE_FN(map_toggle__previous_lt_net_Network){
  ModelTerm *mtp = auxnet->mtp;
  GET_AUX_STORAGE_NUM(StoreTimeAndLasttoggle, dur_inf, 1);
  // See comments for u__previous_lt_net_Network() for explanation.
  MAP_TOGGLE_PROPAGATE_IF(!dur_inf->ticktock);
}

#endif // _TERGM_CHANGESTATS_AUXNET_H_
