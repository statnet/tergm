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
  const double PT = 0.5;
  double logratio = 0;

  DyadGen *formgen = storage[0], *dissgen = storage[1];

  if(unif_rand() < PT){ // Propose from the formables.
    DyadGenWake(formgen);
    DyadGenSleep(dissgen);
    const double P=0.5, Q=1-P;
    double DP = P*formgen->ndyads, DO = DP/Q;

    Edge nedges = DyadGenEdgecount(formgen);
    BD_LOOP({
        if (unif_rand() < P && nedges > 0) { /* Select a tie at random from the network of eligibles */
          DyadGenRandEdge(Mtail, Mhead, formgen);
          logratio = TNT_LR_E(nedges, Q, DP, DO);
        }else{ /* Select a dyad at random from the list */
          DyadGenRandDyad(Mtail, Mhead, formgen);
          
          if(IS_OUTEDGE(Mtail[0],Mhead[0])){
            logratio = TNT_LR_DE(nedges, Q, DP, DO);
          }else{
            logratio = TNT_LR_DN(nedges, Q, DP, DO);
          }
        }
      });

  }else{ // Propose from the dissolvables.
    DyadGenWake(dissgen);
    DyadGenSleep(formgen);
    const double P=0.5, Q=1-P;
    double DP = P*dissgen->ndyads, DO = DP/Q;

    Edge nedges = DyadGenEdgecount(dissgen);
    BD_LOOP({
        if (unif_rand() < P && nedges > 0) { /* Select a tie at random from the network of eligibles */
          DyadGenRandEdge(Mtail, Mhead, dissgen);
          logratio = TNT_LR_E(nedges, Q, DP, DO);
        }else{ /* Select a dyad at random from the list */
          DyadGenRandDyad(Mtail, Mhead, dissgen);
          
          if(IS_OUTEDGE(Mtail[0],Mhead[0])){
            logratio = TNT_LR_DE(nedges, Q, DP, DO);
          }else{
            logratio = TNT_LR_DN(nedges, Q, DP, DO);
          }
        }
      });
  }
  MHp->logratio += logratio;
}


MH_F_FN(Mf_staticDiscordTNT) {
  GET_STORAGE(DyadGen *, storage);

  DyadGenDestroy(storage[0]);
  DyadGenDestroy(storage[1]);
}
