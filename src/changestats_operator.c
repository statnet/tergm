#include "tergm_changestats_auxnet.h"

I_CHANGESTAT_FN(i_on_union_lt_net_Network){
  GET_STORAGE(Model, m);
  GET_AUX_STORAGE(StoreAuxnet, auxnet);
  STORAGE = m = ModelInitialize(getListElement(mtp->R, "submodel"),  NULL, auxnet->onwp, FALSE);
}

C_CHANGESTAT_FN(c_on_union_lt_net_Network){
  GET_STORAGE(Model, m);
  GET_AUX_STORAGE(StoreAuxnet, auxnet);

  Vertex tails[map_toggle_maxtoggles__union_lt_net_Network], heads[map_toggle_maxtoggles__union_lt_net_Network];
  MAP_TOGGLE_THEN(_union_lt_net_Network, tail, head, edgeflag, auxnet, tails, heads){
    double *tmp = m->workspace;
    m->workspace = CHANGE_STAT;
    ChangeStats1(*tails, *heads, auxnet->onwp, m, IS_OUTEDGE(*tails, *heads, auxnet->onwp));
    m->workspace = tmp;
  }
}

U_CHANGESTAT_FN(u_on_union_lt_net_Network){
  GET_STORAGE(Model, m);
  GET_AUX_STORAGE(StoreAuxnet, auxnet);

  Vertex tails[map_toggle_maxtoggles__union_lt_net_Network], heads[map_toggle_maxtoggles__union_lt_net_Network];
  MAP_TOGGLE_THEN(_union_lt_net_Network, tail, head, edgeflag, auxnet, tails, heads){
    UPDATE_STORAGE(*tails, *heads, auxnet->onwp, m, NULL, IS_OUTEDGE(*tails, *heads, auxnet->onwp));
  }
}

Z_CHANGESTAT_FN(z_on_union_lt_net_Network){
  GET_STORAGE(Model, m);
  GET_AUX_STORAGE_NUM(StoreAuxnet, prevnet, 2);

  double *tmp = m->workspace;
  m->workspace = CHANGE_STAT;
  SummStats(0, NULL, NULL, prevnet->onwp, m);
  m->workspace = tmp;
}

F_CHANGESTAT_FN(f_on_union_lt_net_Network){
  GET_STORAGE(Model, m);
  GET_AUX_STORAGE(StoreAuxnet, auxnet);
  ModelDestroy(auxnet->onwp, m);
  STORAGE = NULL;
}

ON_AUXNET(_intersect_lt_net_Network)
