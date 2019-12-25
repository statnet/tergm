#include "tergm_model.h"
#include "ergm_storage.h"
#include "changestats_lasttoggle.h"
#include "ergm_Rutil.h"

/* _lasttoggle */

I_CHANGESTAT_FN(i__lasttoggle){
  ALLOC_AUX_STORAGE(1, StoreTimeAndLasttoggle, dur_inf);
  dur_inf->time = asInteger(getListElement(mtp->ext_state, "time"));
  dur_inf->lasttoggle = kh_init(DyadMapInt);
  dur_inf->discord = kh_init(DyadMapInt);
  dur_inf->discord->directed = dur_inf->lasttoggle->directed=DIRECTED;
  dur_inf->ticktock = FALSE;
  SEXP ltR = getListElement(mtp->ext_state, "lasttoggle");
  Edge nlt = length(ltR)/3, *lt = (Edge *) INTEGER(ltR);
  for(Edge i = 0; i < nlt; i++){
    Vertex tail=lt[i], head=lt[i+nlt];
    /* Note: we can't use helper macros here, since those treat 0 as deletion. */
    kh_set(DyadMapInt,dur_inf->lasttoggle,THKey(dur_inf->lasttoggle,tail,head), lt[i+nlt+nlt]);
  }
}

X_CHANGESTAT_FN(x__lasttoggle){
  switch(type){
  case TICK: // Start time step
    {
      GET_AUX_STORAGE(StoreTimeAndLasttoggle, dur_inf);
      if(dur_inf->ticktock) error("_lasttoggle was sent two TICK signals without a TOCK signal.");
      dur_inf->time++; // Advance the clock.
      if(dur_inf->time%TIMESTAMP_HORIZON_FREQ == 0) ExpireTimestamps(dur_inf, TIMESTAMP_HORIZON_EDGE, TIMESTAMP_HORIZON_NONEDGE, nwp);
      dur_inf->ticktock = TRUE; // We are in the middle of a time step.
    }
    break;
  case TOCK: // Finish time step
    {
      GET_AUX_STORAGE(StoreTimeAndLasttoggle, dur_inf);
      if(!dur_inf->ticktock) error("_lasttoggle was sent a TOCK signal without a previous TICK signal.");
      kh_clear(DyadMapInt, dur_inf->discord);
      dur_inf->ticktock = FALSE; // The time step is completed.
    }
    break;
  default:
    break;
  }
}

U_CHANGESTAT_FN(u__lasttoggle){
  GET_AUX_STORAGE(StoreTimeAndLasttoggle, dur_inf);
  if(!dur_inf->ticktock) return; // If the network is being modified outside of time, do not update lasttoggle or discord.

  TailHead dyad = THKey(dur_inf->discord,tail, head);
  khint_t i = kh_get(DyadMapInt,dur_inf->discord,dyad);
  if(i == kh_none){
    // If the dyad has *not* been toggled in this time step, then save its last toggle info (if any) and make a provisional change in lasttoggle.
    // Here, current time is used as a placeholder for no last toggle info.
    kh_set(DyadMapInt, dur_inf->discord, dyad, kh_getval(DyadMapInt, dur_inf->lasttoggle, dyad, dur_inf->time));
    kh_set(DyadMapInt, dur_inf->lasttoggle, dyad, dur_inf->time);
  }else{
    // If the dyad *has* been toggled in this timestep, then untoggle it by restoring its change in lasttoggle.
    if(kh_value(dur_inf->discord, i) != dur_inf->time){
      kh_set(DyadMapInt, dur_inf->lasttoggle, dyad, kh_value(dur_inf->discord, i));
    }else{
      kh_unset(DyadMapInt, dur_inf->lasttoggle, dyad);
    }

    kh_del(DyadMapInt, dur_inf->discord, i);
  }
}

W_CHANGESTAT_FN(w__lasttoggle){
  GET_AUX_STORAGE(StoreTimeAndLasttoggle, dur_inf);
  if(dur_inf->ticktock) error("_lasttoggle auxiliary cannot write out its extended state in the middle of a time step.");

  SEXP es = PROTECT(mkNamed(VECSXP, (const char *[]){"time","lasttoggle",""}));
  SET_VECTOR_ELT(es, 0, ScalarInteger(dur_inf->time));

  {
    Edge nlt = kh_size(dur_inf->lasttoggle);
    SEXP ltR = PROTECT(allocMatrix(INTSXP,nlt,3));
    int *lt = INTEGER(ltR);
    TailHead dyad;
    int ts;
    unsigned int i=0;
    kh_foreach(dur_inf->lasttoggle, dyad, ts, {
        lt[i] = dyad.tail;
        lt[i+nlt] = dyad.head;
        lt[i+nlt+nlt] = ts;
        i++;
      });
    SET_VECTOR_ELT(es, 1, ltR);
    UNPROTECT(1);
  }

  UNPROTECT(1);
  return es;
}


F_CHANGESTAT_FN(f__lasttoggle){
  GET_AUX_STORAGE(StoreTimeAndLasttoggle, dur_inf);
  kh_destroy(DyadMapInt, dur_inf->lasttoggle);
  kh_destroy(DyadMapInt, dur_inf->discord);
}

/*****************
 void ExpireTimestamps

 Walks through the lasttoggle structure, expiring edges and non-edges
 from the lasttoggle structure that had been toggled more than the
 specified number of steps ago.
 *****************/
void ExpireTimestamps(StoreTimeAndLasttoggle *dur_inf, unsigned int edges, unsigned int nonedges, Network *nwp){
  if(dur_inf->ticktock) error("_lasttoggle auxiliary cannot expire time stamps in the middle of a time step.");
  if(edges==nonedges){ // Same horizon for edges and non-edges means that we don't need to check if an edge exists.
    int lt;
    kh_foreach_value(dur_inf->lasttoggle, lt, {
        /* Note: This bit is implementation-dependent, relying on the
           fact that __i is the iterator variable. If we ever change to
           a backend different from khash, we would need to implement
           this differently. */
        if(dur_inf->time - lt > edges)
          kh_del(DyadMapInt, dur_inf->lasttoggle, __i);
      });
  }else{
    TailHead dyad;
    int lt;
    kh_foreach(dur_inf->lasttoggle, dyad, lt, {
        /* Note: This bit is implementation-dependent, relying on the
           fact that __i is the iterator variable. If we ever change to
           a backend different from khash, we would need to implement
           this differently. */
        if(dur_inf->time - lt > (EdgetreeSearch(dyad.tail,dyad.head,nwp->outedges) ? edges : nonedges))
          kh_del(DyadMapInt, dur_inf->lasttoggle, __i);
      });
  }
}

