#ifndef _TERGM_CHANGESTATS_AUXNET_H_
#define _TERGM_CHANGESTATS_AUXNET_H_

#include "changestats_lasttoggle.h"
#include "ergm_changestat_auxnet.h"

#define map_toggle_maxtoggles__intersect_lt_net_Network 1
MAP_TOGGLE_FN(map_toggle__intersect_lt_net_Network){
  ModelTerm *mtp = auxnet->mtp;
  GET_AUX_STORAGE_NUM(StoreTimeAndLasttoggle, dur_inf, 1);
  // See comments for u__intersect_lt_net_Network() for explanation.
  const int t = dur_inf->time;
  MAP_TOGGLE_PROPAGATE_IF(edgeflag != (kh_getval(DyadMapInt, dur_inf->lasttoggle, THKey(dur_inf->lasttoggle,tail,head), t+1)==t));
}

#define map_toggle_maxtoggles__union_lt_net_Network 1
MAP_TOGGLE_FN(map_toggle__union_lt_net_Network){
  ModelTerm *mtp = auxnet->mtp;
  GET_AUX_STORAGE_NUM(StoreTimeAndLasttoggle, dur_inf, 1);
  // See comments for u__union_lt_net_Network() for explanation.
  const int t = dur_inf->time;
  MAP_TOGGLE_PROPAGATE_IF(edgeflag == (kh_getval(DyadMapInt, dur_inf->lasttoggle, THKey(dur_inf->lasttoggle,tail,head), t+1)==t));
}

#endif // _TERGM_CHANGESTATS_AUXNET_H_
