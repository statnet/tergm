#ifndef _TERGM_CHANGESTATS_AUXNET_H_
#define _TERGM_CHANGESTATS_AUXNET_H_

#include "changestats_lasttoggle.h"
#include "ergm_changestat_auxnet.h"

MAP_TOGGLE_MAXTOGGLES_FN(map_toggle_maxtoggles__intersect_lt_net_Network){
  return 1;
}

MAP_TOGGLE_FN(map_toggle__intersect_lt_net_Network){
  ModelTerm *mtp = auxnet->mtp;
  GET_AUX_STORAGE(StoreTimeAndLasttoggle, dur_inf);
  MAP_TOGGLE_PROPAGATE_IF(edgeflag != (kh_get(DyadMapInt, dur_inf->discord, THKey(dur_inf->discord,tail,head))!=kh_none));
}

MAP_TOGGLE_MAXTOGGLES_FN(map_toggle_maxtoggles__union_net_Network){
  return 1;
}

MAP_TOGGLE_FN(map_toggle__union_lt_net_Network){
  ModelTerm *mtp = auxnet->mtp;
  GET_AUX_STORAGE(StoreTimeAndLasttoggle, dur_inf);
  MAP_TOGGLE_PROPAGATE_IF(edgeflag == (kh_get(DyadMapInt, dur_inf->discord, THKey(dur_inf->discord,tail,head))!=kh_none));
}

#endif // _TERGM_CHANGESTATS_AUXNET_H_
