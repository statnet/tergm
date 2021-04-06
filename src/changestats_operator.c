/*  File src/changestats_operator.c in package tergm, part of the Statnet suite
 *  of packages for network analysis, https://statnet.org .
 *
 *  This software is distributed under the GPL-3 license.  It is free,
 *  open source, and has the attribution requirements (GPL Section 7) at
 *  https://statnet.org/attribution
 *
 *  Copyright 2008-2020 Statnet Commons
 */
#include "tergm_model.h"
#include "tergm_changestats_auxnet.h"

I_CHANGESTAT_FN(i_on_union_lt_net_Network){
  GET_STORAGE(Model, m);
  GET_AUX_STORAGE(StoreAuxnet, auxnet);
  STORAGE = m = ModelInitialize(getListElement(mtp->R, "submodel"),  NULL, auxnet->onwp, FALSE);
  DELETE_IF_UNUSED_IN_SUBMODEL(u_func, m);
}

C_CHANGESTAT_FN(c_on_union_lt_net_Network){
  GET_STORAGE(Model, m);
  GET_AUX_STORAGE(StoreAuxnet, auxnet);

  Vertex tails[map_toggle_maxtoggles__union_lt_net_Network], heads[map_toggle_maxtoggles__union_lt_net_Network];
  MAP_TOGGLE_THEN(_union_lt_net_Network, tail, head, edgestate, auxnet, tails, heads){
    double *tmp = m->workspace;
    m->workspace = CHANGE_STAT;
    ChangeStats1(*tails, *heads, auxnet->onwp, m, IS_OUTEDGE(*tails, *heads, auxnet->onwp));
    m->workspace = tmp;
  }
}

X_CHANGESTAT_FN(x_on_union_lt_net_Network){
  switch(type){
  case TICK:
    {
      GET_STORAGE(Model, m);
      GET_AUX_STORAGE(StoreAuxnet, auxnet);
      GET_AUX_STORAGE_NUM(StoreTimeAndLasttoggle, dur_inf, 1);
      TailHead dyad;
      // Here, we want (y0|y1) / y1: edges in y0 but not in y1.
      // TODO: Optimize.
      unsigned int nt=kh_size(dur_inf->discord), pos=0;
      Vertex *tails = Calloc(nt, Vertex);
      Vertex *heads = Calloc(nt, Vertex);
      kh_foreach_key(dur_inf->discord, dyad, {
          if(!IS_OUTEDGE(dyad.tail, dyad.head)){
            tails[pos] = dyad.tail;
            heads[pos] = dyad.head;
            pos++;
          }
        });

      ChangeStats(pos, tails, heads, auxnet->onwp, m);
      memcpy(CHANGE_STAT, m->workspace, m->n_stats*sizeof(double));
      Free(tails); Free(heads);
    }
    break;
  default: break;
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





I_CHANGESTAT_FN(i_on_intersect_lt_net_Network){
  GET_STORAGE(Model, m);
  GET_AUX_STORAGE(StoreAuxnet, auxnet);
  STORAGE = m = ModelInitialize(getListElement(mtp->R, "submodel"),  NULL, auxnet->onwp, FALSE);
  DELETE_IF_UNUSED_IN_SUBMODEL(u_func, m);
}

C_CHANGESTAT_FN(c_on_intersect_lt_net_Network){
  GET_STORAGE(Model, m);
  GET_AUX_STORAGE(StoreAuxnet, auxnet);

  Vertex tails[map_toggle_maxtoggles__intersect_lt_net_Network], heads[map_toggle_maxtoggles__intersect_lt_net_Network];
  MAP_TOGGLE_THEN(_intersect_lt_net_Network, tail, head, edgestate, auxnet, tails, heads){
    double *tmp = m->workspace;
    m->workspace = CHANGE_STAT;
    ChangeStats1(*tails, *heads, auxnet->onwp, m, IS_OUTEDGE(*tails, *heads, auxnet->onwp));
    m->workspace = tmp;
  }
}

X_CHANGESTAT_FN(x_on_intersect_lt_net_Network){
  switch(type){
  case TICK:
    {
      GET_STORAGE(Model, m);
      GET_AUX_STORAGE(StoreAuxnet, auxnet);
      GET_AUX_STORAGE_NUM(StoreTimeAndLasttoggle, dur_inf, 1);
      TailHead dyad;
      // Here, we want (y0&y1) / y0: edges in y1 but not in y0.
      // TODO: Optimize.
      unsigned int nt=kh_size(dur_inf->discord), pos=0;
      Vertex *tails = Calloc(nt, Vertex);
      Vertex *heads = Calloc(nt, Vertex);
      kh_foreach_key(dur_inf->discord, dyad, {
          if(IS_OUTEDGE(dyad.tail, dyad.head)){
            tails[pos] = dyad.tail;
            heads[pos] = dyad.head;
            pos++;
          }
        });

      ChangeStats(pos, tails, heads, auxnet->onwp, m);
      memcpy(CHANGE_STAT, m->workspace, m->n_stats*sizeof(double));
      Free(tails); Free(heads);
    }
    break;
  default: break;
  }
}

Z_CHANGESTAT_FN(z_on_intersect_lt_net_Network){
  GET_STORAGE(Model, m);
  GET_AUX_STORAGE(StoreAuxnet, auxnet);

  double *tmp = m->workspace;
  m->workspace = CHANGE_STAT;
  ZStats(auxnet->onwp, m, FALSE);
  m->workspace = tmp;
}

F_CHANGESTAT_FN(f_on_intersect_lt_net_Network){
  GET_STORAGE(Model, m);
  GET_AUX_STORAGE(StoreAuxnet, auxnet);
  ModelDestroy(auxnet->onwp, m);
  STORAGE = NULL;
}
