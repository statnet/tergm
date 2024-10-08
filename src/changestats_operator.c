/*  File src/changestats_operator.c in package tergm, part of the
 *  Statnet suite of packages for network analysis, https://statnet.org .
 *
 *  This software is distributed under the GPL-3 license.  It is free,
 *  open source, and has the attribution requirements (GPL Section 7) at
 *  https://statnet.org/attribution .
 *
 *  Copyright 2008-2024 Statnet Commons
 */
#include "ergm_util.h"
#include "tergm_model.h"
#include "tergm_changestats_auxnet.h"

I_CHANGESTAT_FN(i_on_union_lt_net_Network){
  GET_STORAGE(Model, m);
  GET_AUX_STORAGE(StoreAuxnet, auxnet);
  STORAGE = m = ModelInitialize(getListElement(mtp->R, "submodel"), isNULL(mtp->ext_state) ? NULL : mtp->ext_state, auxnet->onwp, FALSE);
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
  GET_STORAGE(Model, m);
  GET_AUX_STORAGE(StoreAuxnet, auxnet);

  switch(type){
  case TICK:
    {
      GET_AUX_STORAGE_NUM(StoreTimeAndLasttoggle, dur_inf, 1);
      TailHead dyad;
      // Here, we want (y0|y1) / y1: edges in y0 but not in y1.
      // TODO: Optimize.
      unsigned int nt=kh_size(dur_inf->discord), pos=0;
      Vertex *tails = R_Calloc(nt, Vertex);
      Vertex *heads = R_Calloc(nt, Vertex);
      kh_foreach_key(dur_inf->discord, dyad, {
          if(!IS_OUTEDGE(dyad.tail, dyad.head)){
            tails[pos] = dyad.tail;
            heads[pos] = dyad.head;
            pos++;
          }
        });

      ChangeStats(pos, tails, heads, auxnet->onwp, m);
      memcpy(CHANGE_STAT, m->workspace, m->n_stats*sizeof(double));
      R_Free(tails); R_Free(heads);
    }
    break;
  default: break;
  }

  PROPAGATE_X_SIGNAL_ADDONTO(auxnet->onwp, m, CHANGE_STAT);
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
  STORAGE = m = ModelInitialize(getListElement(mtp->R, "submodel"), isNULL(mtp->ext_state) ? NULL : mtp->ext_state, auxnet->onwp, FALSE);
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
  GET_STORAGE(Model, m);
  GET_AUX_STORAGE(StoreAuxnet, auxnet);

  switch(type){
  case TICK:
    {
      GET_AUX_STORAGE_NUM(StoreTimeAndLasttoggle, dur_inf, 1);
      TailHead dyad;
      // Here, we want (y0&y1) / y0: edges in y1 but not in y0.
      // TODO: Optimize.
      unsigned int nt=kh_size(dur_inf->discord), pos=0;
      Vertex *tails = R_Calloc(nt, Vertex);
      Vertex *heads = R_Calloc(nt, Vertex);
      kh_foreach_key(dur_inf->discord, dyad, {
          if(IS_OUTEDGE(dyad.tail, dyad.head)){
            tails[pos] = dyad.tail;
            heads[pos] = dyad.head;
            pos++;
          }
        });

      ChangeStats(pos, tails, heads, auxnet->onwp, m);
      memcpy(CHANGE_STAT, m->workspace, m->n_stats*sizeof(double));
      R_Free(tails); R_Free(heads);
    }
    break;
  default: break;
  }

  PROPAGATE_X_SIGNAL_ADDONTO(auxnet->onwp, m, CHANGE_STAT);
}

Z_CHANGESTAT_FN(z_on_intersect_lt_net_Network){
  GET_STORAGE(Model, m);
  GET_AUX_STORAGE(StoreAuxnet, auxnet);

  // Note: If we ever implement an s_ function for this operator, the
  // last arguments of EmptyNetworkStats() and ZStats() should
  // probably be changed to skip_s.
  double *tmp = m->workspace;
  m->workspace = CHANGE_STAT;
  EmptyNetworkStats(m, FALSE);
  m->workspace = tmp;
  ZStats(auxnet->onwp, m, FALSE);
  addonto(CHANGE_STAT, m->workspace, N_CHANGE_STATS);
}

F_CHANGESTAT_FN(f_on_intersect_lt_net_Network){
  GET_STORAGE(Model, m);
  GET_AUX_STORAGE(StoreAuxnet, auxnet);
  ModelDestroy(auxnet->onwp, m);
  STORAGE = NULL;
}




#define NEGATE_CHANGESTAT for(unsigned int i = 0; i < N_CHANGE_STATS; i++) CHANGE_STAT[i] = - CHANGE_STAT[i];

I_CHANGESTAT_FN(i_negate_on_intersect_lt_net_Network){
  i_on_intersect_lt_net_Network(mtp, nwp);
}

C_CHANGESTAT_FN(c_negate_on_intersect_lt_net_Network){
  c_on_intersect_lt_net_Network(tail, head, mtp, nwp, edgestate);
  NEGATE_CHANGESTAT;
}

X_CHANGESTAT_FN(x_negate_on_intersect_lt_net_Network){
  x_on_intersect_lt_net_Network(type, data, mtp, nwp);
  NEGATE_CHANGESTAT;
}

Z_CHANGESTAT_FN(z_negate_on_intersect_lt_net_Network){
  z_on_intersect_lt_net_Network(mtp, nwp, skip_s);
  NEGATE_CHANGESTAT;
}

F_CHANGESTAT_FN(f_negate_on_intersect_lt_net_Network){
  f_on_intersect_lt_net_Network(mtp, nwp);
}

#undef MINUS_WS_TO_CHANGESTAT



I_CHANGESTAT_FN(i_on_discord_lt_net_Network){
  GET_STORAGE(Model, m);
  GET_AUX_STORAGE(StoreAuxnet, auxnet);
  STORAGE = m = ModelInitialize(getListElement(mtp->R, "submodel"), isNULL(mtp->ext_state) ? NULL : mtp->ext_state, auxnet->onwp, FALSE);
}

C_CHANGESTAT_FN(c_on_discord_lt_net_Network){
  GET_STORAGE(Model, m);
  GET_AUX_STORAGE(StoreAuxnet, auxnet);

  Vertex tails[map_toggle_maxtoggles__discord_lt_net_Network], heads[map_toggle_maxtoggles__discord_lt_net_Network];
  MAP_TOGGLE_THEN(_discord_lt_net_Network, tail, head, edgestate, auxnet, tails, heads){
    double *tmp = m->workspace;
    m->workspace = CHANGE_STAT;
    ChangeStats1(*tails, *heads, auxnet->onwp, m, IS_OUTEDGE(*tails, *heads, auxnet->onwp));
    m->workspace = tmp;
  }
}

X_CHANGESTAT_FN(x_on_discord_lt_net_Network){
  GET_STORAGE(Model, m);
  GET_AUX_STORAGE(StoreAuxnet, auxnet);

  switch(type){
  case TICK:
    {
      GET_AUX_STORAGE_NUM(StoreTimeAndLasttoggle, dur_inf, 1);
      TailHead dyad;
      // TODO: Optimize.
      unsigned int nt=kh_size(dur_inf->discord), pos=0;
      Vertex *tails = R_Calloc(nt, Vertex);
      Vertex *heads = R_Calloc(nt, Vertex);
      kh_foreach_key(dur_inf->discord, dyad, {
          tails[pos] = dyad.tail;
          heads[pos] = dyad.head;
          pos++;
        });

      ChangeStats(pos, tails, heads, auxnet->onwp, m);
      memcpy(CHANGE_STAT, m->workspace, m->n_stats*sizeof(double));
      R_Free(tails); R_Free(heads);
    }
    break;
  default: break;
  }

  PROPAGATE_X_SIGNAL_ADDONTO(auxnet->onwp, m, CHANGE_STAT);
}

Z_CHANGESTAT_FN(z_on_discord_lt_net_Network){
  GET_STORAGE(Model, m);
  GET_AUX_STORAGE_NUM(StoreAuxnet, prevnet, 2);

  double *tmp = m->workspace;
  m->workspace = CHANGE_STAT;
  SummStats(0, NULL, NULL, prevnet->onwp, m);
  m->workspace = tmp;
}

F_CHANGESTAT_FN(f_on_discord_lt_net_Network){
  GET_STORAGE(Model, m);
  GET_AUX_STORAGE(StoreAuxnet, auxnet);
  ModelDestroy(auxnet->onwp, m);
  STORAGE = NULL;
}

/* Take a subset of the network statistics of the submodel and drop the rest. */

#define _REMAP_WS_ONTO_CS_                                              \
  for(unsigned int i = 0; i < N_CHANGE_STATS; i++)                      \
    CHANGE_STAT[i] = m->workspace[IINPUT_PARAM[i]]; // IINPUT_PARAM = the mapping.

I_CHANGESTAT_FN(i_subset_stats) {
  GET_STORAGE(Model, m);
  STORAGE = m = ModelInitialize(getListElement(mtp->R, "submodel"), mtp->ext_state, nwp, FALSE);

  SELECT_C_OR_D_BASED_ON_SUBMODEL(m);
  DELETE_IF_UNUSED_IN_SUBMODEL(x_func, m);
  DELETE_IF_UNUSED_IN_SUBMODEL(z_func, m);
}

D_CHANGESTAT_FN(d_subset_stats) {
  GET_STORAGE(Model, m);
  ChangeStats(ntoggles, tails, heads, nwp, m);
  _REMAP_WS_ONTO_CS_
}

C_CHANGESTAT_FN(c_subset_stats) {
  GET_STORAGE(Model, m);
  ChangeStats1(tail, head, nwp, m, edgestate);
  _REMAP_WS_ONTO_CS_
}

X_CHANGESTAT_FN(x_subset_stats) {
  GET_STORAGE(Model, m);
  PROPAGATE_X_SIGNAL_INTO(nwp, m, m->workspace);
  _REMAP_WS_ONTO_CS_
}

Z_CHANGESTAT_FN(z_subset_stats) {
  GET_STORAGE(Model, m);
  ZStats(nwp, m, FALSE);
  _REMAP_WS_ONTO_CS_
}

F_CHANGESTAT_FN(f_subset_stats) {
  GET_STORAGE(Model, m);
  ModelDestroy(nwp, m);
  STORAGE = NULL;
}

#undef _REMAP_WS_ONTO_CS_



/*****************
 EdgeAges

 Sum of (tie age) * (submodel on-toggle change stats) for all extant ties.

*****************/

typedef struct {
  Model *model;
  double *stats;
} EdgeAges_storage;

I_CHANGESTAT_FN(i_EdgeAges) {
  ALLOC_STORAGE(1, EdgeAges_storage, sto);
  sto->model = ModelInitialize(getListElement(mtp->R, "submodel"), mtp->ext_state, nwp, FALSE);
  sto->stats = R_Calloc(N_CHANGE_STATS, double);

  EXEC_THROUGH_NET_EDGES_PRE(tail, head, edge_var, {
    ChangeStats1(tail, head, nwp, sto->model, edge_var);

    for(int i = 0; i < N_CHANGE_STATS; i++) {
      sto->stats[i] -= sto->model->workspace[i];
    }
  });
}

X_CHANGESTAT_FN(x_EdgeAges) {
  GET_STORAGE(EdgeAges_storage, sto);

  if(type == TICK) {
    memcpy(CHANGE_STAT, sto->stats, N_CHANGE_STATS*sizeof(double));
  }

  // ignoring any change from this...
  PROPAGATE_X_SIGNAL(nwp, sto->model);
}

C_CHANGESTAT_FN(c_EdgeAges) {
  GET_STORAGE(EdgeAges_storage, sto);
  GET_AUX_STORAGE(StoreTimeAndLasttoggle, dur_inf);

  ChangeStats1(tail, head, nwp, sto->model, edgestate);
  int age = ElapsedTimeToggle(tail, head, dur_inf, edgestate) + 1;

  for(int i = 0; i < N_CHANGE_STATS; i++) {
    CHANGE_STAT[i] = age*sto->model->workspace[i];
  }
}

U_CHANGESTAT_FN(u_EdgeAges) {
  GET_STORAGE(EdgeAges_storage, sto);
  ChangeStats1(tail, head, nwp, sto->model, edgestate);
  for(int i = 0; i < N_CHANGE_STATS; i++) {
    sto->stats[i] += sto->model->workspace[i];
  }
}

S_CHANGESTAT_FN(s_EdgeAges) {
  Model *m = ModelInitialize(getListElement(mtp->R, "submodel"), mtp->ext_state, nwp, FALSE);
  GET_AUX_STORAGE(StoreTimeAndLasttoggle, dur_inf);

  EXEC_THROUGH_NET_EDGES_PRE(tail, head, edge_var, {
    ChangeStats1(tail, head, nwp, m, edge_var);
    int age = ElapsedTime(tail, head, dur_inf) + 1;
    for(int i = 0; i < N_CHANGE_STATS; i++) {
      CHANGE_STAT[i] -= age*m->workspace[i];
    }
  });
}

// no Z_FN as emptynwstats = 0 for EdgeAges, regardless of submodel

F_CHANGESTAT_FN(f_EdgeAges) {
  GET_STORAGE(EdgeAges_storage, sto);
  ModelDestroy(nwp, sto->model);
  R_Free(sto->stats);
}
