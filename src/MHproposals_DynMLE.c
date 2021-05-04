/*  File src/MHproposals_DynMLE.c in package tergm, part of the Statnet suite
 *  of packages for network analysis, https://statnet.org .
 *
 *  This software is distributed under the GPL-3 license.  It is free,
 *  open source, and has the attribution requirements (GPL Section 7) at
 *  https://statnet.org/attribution
 *
 *  Copyright 2008-2020 Statnet Commons
 */
#include "ergm_MHproposal.h"
#include "ergm_MHproposal_bd.h"
#include "ergm_MHstorage.h"
#include "ergm_changestat.h"
#include "ergm_Rutil.h"
#include "ergm_dyadgen.h"

typedef struct{
  DyadGen *gen[2];
  DegreeBound *bd;
} Store2DyadGenAndDegreeBound;

MH_I_FN(Mi_staticDiscordTNT){
  MHp->ntoggles = 1;

  ALLOC_STORAGE(1, Store2DyadGenAndDegreeBound, storage);
  storage->gen[0] = DyadGenInitializeR(getListElement(MHp->R, "formable"), nwp, TRUE);
  storage->gen[1] = DyadGenInitializeR(getListElement(MHp->R, "dissolvable"), nwp, TRUE);
  storage->bd = DegreeBoundInitializeR(MHp->R, nwp);

  if(storage->gen[0]->ndyads==0 && storage->gen[1]->ndyads==0) error("At least one of the dyad sets has to have toggleable dyads.");
}


MH_P_FN(MH_staticDiscordTNT) {
  GET_STORAGE(Store2DyadGenAndDegreeBound, storage);
  double logratio = 0;

  DyadGen *activegen, *inactivegen;

  if(storage->gen[1]->ndyads==0 || (storage->gen[0]->ndyads!=0 && unif_rand()<MH_INPUTS[0])){ // Propose from the formables.
    activegen = storage->gen[0];
    inactivegen = storage->gen[1];
  }else{ // Propose from the dissolvables.
    activegen = storage->gen[1];
    inactivegen = storage->gen[0];
  }

  DyadGenWake(activegen);
  DyadGenSleep(inactivegen);

  const double P=0.5, Q=1-P;
  double DP = P*activegen->ndyads, DO = DP/Q;
  Edge nedges = DyadGenEdgecount(activegen);

  BD_LOOP(storage->bd, {
      if(unif_rand() < P && nedges > 0){ /* Select a tie at random from the network of eligibles */
        DyadGenRandEdge(Mtail, Mhead, activegen);
        logratio = TNT_LR_E(nedges, Q, DP, DO);
      }else{ /* Select a dyad at random from the list */
        DyadGenRandDyad(Mtail, Mhead, activegen);

        if(IS_OUTEDGE(Mtail[0],Mhead[0])){
          logratio = TNT_LR_DE(nedges, Q, DP, DO);
        }else{
          logratio = TNT_LR_DN(nedges, Q, DP, DO);
        }
      }
    });

  MHp->logratio += logratio;
}


MH_F_FN(Mf_staticDiscordTNT) {
  GET_STORAGE(Store2DyadGenAndDegreeBound, storage);

  DyadGenDestroy(storage->gen[0]);
  DyadGenDestroy(storage->gen[1]);
  DegreeBoundDestroy(storage->bd);
}
