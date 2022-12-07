
/********************
    MH_maxTNT
********************/
typedef struct {
  double ne;
  Dyad nd;
} maxTNTstorage;

MH_I_FN(Mi_maxTNT) {
  MHp->ntoggles = 1;
  ALLOC_STORAGE(1, maxTNTstorage, sto);
  sto->ne = asReal(getListElement(MHp->R, "n"));
  sto->nd = DYADCOUNT(nwp);
}

MH_P_FN(MH_maxTNT) {
  GET_STORAGE(maxTNTstorage, sto);

  if(EDGECOUNT(nwp) >= sto->ne) {
    error("edgecount in maxTNT cannot reach or exceed maximum number of edges argument");
  }

  double u = unif_rand() - EDGECOUNT(nwp)/(2*sto->ne);
  if(u < 0) {
    GetRandEdge(Mtail, Mhead, nwp);
  } else if (u < 0.5) {
    GetRandDyad(Mtail, Mhead, nwp);
  } else {
    // do nothing
    MHp->toggletail[0] = MH_FAILED;
    MHp->togglehead[0] = MH_CONSTRAINT;
  }

  // need to set logratio
  if(IS_OUTEDGE(*Mtail, *Mhead)) {
    MHp->logratio = -log(1 + sto->nd/sto->ne);
  } else {
    MHp->logratio = log(1 + sto->nd/sto->ne);
  }
}
