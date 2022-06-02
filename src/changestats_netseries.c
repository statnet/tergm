/*  File src/changestats_netseries.c in package tergm, part of the
 *  Statnet suite of packages for network analysis, https://statnet.org .
 *
 *  This software is distributed under the GPL-3 license.  It is free,
 *  open source, and has the attribution requirements (GPL Section 7) at
 *  https://statnet.org/attribution .
 *
 *  Copyright 2008-2022 Statnet Commons
 */
#include "changestat_netseries.h"
#include "ergm_changestat_operator.h"
#include "ergm_changestat.h"
#include "ergm_model.h"
#include "ergm_storage.h"
#include "ergm_util.h"


I_CHANGESTAT_FN(i__crossnets){
  int *iinputs = IINPUT_PARAM;
  ALLOC_AUX_STORAGE(1, StoreSubnets, sn);
  sn->ns = *(iinputs++);
  sn->inwp = nwp;
  sn->onwp = Calloc(sn->ns, Network *);
  sn->onwp--; // The -- is because Network IDs count from 1.

  /* Set up the layer information. */
  sn->sid = (Vertex *) iinputs - 1; // The -1 is because Vertex IDs count from 1.
  iinputs += N_NODES;
  sn->smap = (Vertex *) iinputs - 1;
  iinputs += N_NODES;

  for(unsigned int i=1; i<=sn->ns; i++){
    Vertex lnnodes, lbip;
    if(BIPARTITE){
      lbip = lnnodes = *(iinputs++);
      lnnodes += *(iinputs++);
    }else{
      lbip = 0;
      lnnodes = *(iinputs++);
    }

    sn->onwp[i] = NetworkInitialize(NULL, NULL, 0, lnnodes, DIRECTED, lbip, 0, 0, NULL);
  }

  EXEC_THROUGH_NET_EDGES_PRE(t, h, e, {
      ToggleKnownEdge(MN_IO_TAIL(sn, t), MN_IO_HEAD(sn, h), sn->onwp[MN_SID_TAIL(sn, t)], FALSE);
    });
}

U_CHANGESTAT_FN(u__crossnets){
  GET_AUX_STORAGE(StoreSubnets, sn);
  ToggleKnownEdge(MN_IO_TAIL(sn, tail), MN_IO_HEAD(sn, head),sn->onwp[MN_SID_TAIL(sn, tail)], edgestate);
}

F_CHANGESTAT_FN(f__crossnets){
  GET_AUX_STORAGE(StoreSubnets, sn);
  for(unsigned int i=1; i<=sn->ns; i++)
    NetworkDestroy(sn->onwp[i]);
  sn->onwp++;
  Free(sn->onwp);
}

// OnCrossNets: Take a networkwise sum of the networks' statistics.

I_CHANGESTAT_FN(i_OnCrossNets){
  GET_AUX_STORAGE(StoreSubnets, sn);
  unsigned int ns = sn->ns;

  ALLOC_STORAGE(ns, Model*, ms);

  SEXP submodels = getListElement(mtp->R, "submodels");
  for(unsigned int i=1; i<=sn->ns; i++){
    ms[i-1] = ModelInitialize(VECTOR_ELT(submodels, i-1), NULL, sn->onwp[i], FALSE);
  }
  DELETE_IF_UNUSED_IN_SUBMODELS(u_func, ms, sn->ns);
  DELETE_IF_UNUSED_IN_SUBMODELS(z_func, ms, sn->ns);
}

C_CHANGESTAT_FN(c_OnCrossNets){
  GET_AUX_STORAGE(StoreSubnets, sn);
  GET_STORAGE(Model*, ms);

  unsigned int i = MN_SID_TAIL(sn, tail);
  Model *m = ms[i-1];
  Vertex st = MN_IO_TAIL(sn, tail), sh = MN_IO_HEAD(sn, head);
  ChangeStats1(st, sh, sn->onwp[i], m, edgestate);
  memcpy(CHANGE_STAT, m->workspace, m->n_stats*sizeof(double));
}

F_CHANGESTAT_FN(f_OnCrossNets){
  GET_AUX_STORAGE(StoreSubnets, sn);
  GET_STORAGE(Model*, ms);

  for(unsigned int i=1; i<=sn->ns; i++){
    ModelDestroy(sn->onwp[i], ms[i-1]);
  }
}

Z_CHANGESTAT_FN(z_OnCrossNets){
  GET_AUX_STORAGE(StoreSubnets, sn);
  GET_STORAGE(Model*, ms);

  for(unsigned int i=1; i<=sn->ns; i++){
    Model *m = ms[i-1];
    ZStats(sn->onwp[i], m, FALSE);
    addonto(CHANGE_STAT, m->workspace, m->n_stats*sizeof(double));
  }
}
