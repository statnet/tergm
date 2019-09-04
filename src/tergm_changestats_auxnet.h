#ifndef _TERGM_CHANGESTATS_AUXNET_H_
#define _TERGM_CHANGESTATS_AUXNET_H_

#include "ergm_changestat_auxnet.h"

MAP_TOGGLE_MAXTOGGLES_FN(map_toggle_maxtoggles__intersect_lt_net_Network){
  return 1;
}

MAP_TOGGLE_FN(map_toggle__intersect_lt_net_Network){
  MAP_TOGGLE_PROPAGATE_IF(HASDMI(tail, head, auxnet->inwp->duration_info->lasttoggle));
}

MAP_TOGGLE_MAXTOGGLES_FN(map_toggle_maxtoggles__union_net_Network){
  return 1;
}

MAP_TOGGLE_FN(map_toggle__union_lt_net_Network){
  MAP_TOGGLE_PROPAGATE_IF(!HASDMI(tail, head, auxnet->inwp->duration_info->lasttoggle));
}

#endif // _TERGM_CHANGESTATS_AUXNET_H_
