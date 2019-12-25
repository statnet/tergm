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

#endif // _CHANGESTATS_LASTTOGGLE_H_
