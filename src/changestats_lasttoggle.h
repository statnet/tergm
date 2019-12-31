#ifndef _CHANGESTATS_LASTTOGGLE_H_
#define _CHANGESTATS_LASTTOGGLE_H_

#include "ergm_changestat.h"
#include "ergm_dyad_hashmap.h"

/* Dur_Inf is a structure containing information about durations of
edges in a network structure.
*/ 
typedef struct StoreTimeAndLasttoggle_struct {
  int time;
  StoreDyadMapInt *lasttoggle;
  StoreDyadMapInt *discord;
  Rboolean ticktock;
} StoreTimeAndLasttoggle;

void ExpireTimestamps(StoreTimeAndLasttoggle *dur_inf, unsigned int edges, unsigned int nonedges, Network *nwp);

/*****************
 long int ElapsedTime

 Return time since given (tail,head) was last toggled using
 ToggleEdgeWithTimestamp function
*****************/

/* *** don't forget tail->head, so this function now accepts tail before head */

static inline int ElapsedTime(Vertex tail, Vertex head, StoreTimeAndLasttoggle *dur_inf){
  khint_t i = kh_get(DyadMapInt, dur_inf->lasttoggle, THKey(dur_inf->lasttoggle,tail,head));
  if(i==kh_none) return INT_MAX;
  return dur_inf->time - kh_value(dur_inf->lasttoggle, i); // Possible int overflow here.
}

#define JUST_CHANGED(dur_inf, tail, head) (kh_getval(DyadMapInt, (dur_inf)->lasttoggle, THKey((dur_inf)->lasttoggle,(tail),(head)), (dur_inf)->time+1)==(dur_inf)->time)

static inline int ElapsedTimeToggle(Vertex tail, Vertex head, StoreTimeAndLasttoggle *dur_inf, Vertex toggletail, Vertex togglehead, int edgeflag){
  if(tail != toggletail || head != togglehead) {
    // if this *isn't* the dyad we're toggling then we can safely use ElapsedTime
    return ElapsedTime(tail, head, dur_inf);
  }
  
  TailHead dyad = THKey(dur_inf->discord,tail, head);
  khint_t i = kh_get(DyadMapInt,dur_inf->discord,dyad);
  if(i == kh_none){
    // not in discord -> even number of toggles this timestep (including the current toggle)
    if(edgeflag) return 0;
    else return ElapsedTime(tail, head, dur_inf);
  }else{
    // in discord -> odd number of toggles this timestep (including the current toggle)
    if(edgeflag) {
      int ltt = kh_value(dur_inf->discord, i);
      if(ltt == dur_inf->time) return INT_MAX;
      else return dur_inf->time - ltt;
    }
    else return 0;
  }
}

#endif // _CHANGESTATS_LASTTOGGLE_H_
