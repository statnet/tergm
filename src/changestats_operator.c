#include "ergm_changestat_operator.h"
#include "ergm_changestats_auxnet.h"

/* formation */
I_CHANGESTAT_FN(i_formation){
  double *inputs = INPUT_ATTRIB;

  GET_AUX_STORAGE(StoreNetAndRefEL, fstorage);
  Model *m = STORAGE = unpack_Model_as_double(&inputs);
  InitStats(fstorage->nwp, m);
}

C_CHANGESTAT_FN(c_formation){
  GET_STORAGE(Model, m);
  GET_AUX_STORAGE(StoreNetAndRefEL, fstorage);

  Rboolean prevflag = dEdgeListSearch(tail, head, fstorage->ref_el);
  if(!prevflag){
    ChangeStats(1, &tail, &head, fstorage->nwp, m);
    memcpy(CHANGE_STAT, m->workspace, N_CHANGE_STATS*sizeof(double));
  }
}

U_CHANGESTAT_FN(u_formation){
  GET_STORAGE(Model, m);
  GET_AUX_STORAGE(StoreNetAndRefEL, fstorage);

  Rboolean prevflag = dEdgeListSearch(tail, head, fstorage->ref_el);
  if(!prevflag) GET_EDGE_UPDATE_STORAGE(tail, head, fstorage->nwp, m, NULL);
}

F_CHANGESTAT_FN(f_formation){
  GET_STORAGE(Model, m);
  GET_AUX_STORAGE(StoreNetAndRefEL, fstorage);
  ModelDestroy(fstorage->nwp, m);
  STORAGE = NULL;
}

/* dissolution */
I_CHANGESTAT_FN(i_dissolution){
  double *inputs = INPUT_ATTRIB;

  GET_AUX_STORAGE(StoreNetAndRefEL, dstorage);
  Model *m = STORAGE = unpack_Model_as_double(&inputs);
  InitStats(dstorage->nwp, m);
}

C_CHANGESTAT_FN(c_dissolution){
  GET_STORAGE(Model, m);
  GET_AUX_STORAGE(StoreNetAndRefEL, dstorage);

  Rboolean prevflag = dEdgeListSearch(tail, head, dstorage->ref_el);
  if(prevflag){
    ChangeStats(1, &tail, &head, dstorage->nwp, m);
    memcpy(CHANGE_STAT, m->workspace, N_CHANGE_STATS*sizeof(double));
  }
}

U_CHANGESTAT_FN(u_dissolution){
  GET_STORAGE(Model, m);
  GET_AUX_STORAGE(StoreNetAndRefEL, dstorage);

  Rboolean prevflag = dEdgeListSearch(tail, head, dstorage->ref_el);
  if(prevflag) GET_EDGE_UPDATE_STORAGE(tail, head, dstorage->nwp, m, NULL);
}

F_CHANGESTAT_FN(f_dissolution){
  GET_STORAGE(Model, m);
  GET_AUX_STORAGE(StoreNetAndRefEL, dstorage);
  ModelDestroy(dstorage->nwp, m);
  STORAGE = NULL;
}
