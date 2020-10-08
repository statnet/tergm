#include "ergm_MHproposal.h"
#include "ergm_MHstorage.h"
#include "ergm_changestat.h"
#include "ergm_Rutil.h"
#include "ergm_dyadgen.h"

MH_I_FN(Mi_staticDiscordTNT){
  MHp->ntoggles = 1;

  ALLOC_STORAGE(2, DyadGen *, storage);
  storage[0] = DyadGenInitializeR(getListElement(MHp->R, "formable"), nwp, TRUE);
  storage[1] = DyadGenInitializeR(getListElement(MHp->R, "dissolvable"), nwp, TRUE);
}


MH_P_FN(MH_staticDiscordTNT) {
  GET_STORAGE(DyadGen *, storage);
  double logratio = 0;

  DyadGen *activegen, *inactivegen;

  if(unif_rand() < MH_INPUTS[0]){ // Propose from the formables.
    activegen = storage[0];
    inactivegen = storage[1];
  }else{ // Propose from the dissolvables.
    activegen = storage[1];
    inactivegen = storage[0];
  }

  DyadGenWake(activegen);
  DyadGenSleep(inactivegen);

  const double P=0.5, Q=1-P;
  double DP = P*activegen->ndyads, DO = DP/Q;
  Edge nedges = DyadGenEdgecount(activegen);

  BD_LOOP({
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
  GET_STORAGE(DyadGen *, storage);

  DyadGenDestroy(storage[0]);
  DyadGenDestroy(storage[1]);
}
