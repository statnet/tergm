#ifndef _CHANGESTATS_LASTTOGGLE_H_
#define _CHANGESTATS_LASTTOGGLE_H_

#include "ergm_changestat.h"
#include "ergm_dyad_hashmap.h"

/* Dur_Inf is a structure containing information about durations of
   edges in a network structure. It is managed by the _lasttoggle
   auxiliary.

   This data structure can be in one of two states, controlled by TICK
   and TOCK signals sent to the term and stored in element .ticktock,
   which is always initialized to FALSE. When a TICK signal is
   received, .ticktock is set to TRUE, .time is incremented, and
   .discord is cleared; and when TOCK signal is received, it is set to
   FALSE. The invariants of the two states are described in turn.

   .ticktock=FALSE:

   Changes to the network are assumed to be occurring "outside of
   time", such as when constructing the initial network. Thus, the
   .lasttoggle and .discord mappings are not updated, and client
   functions should treat a toggle in the network as toggling the
   current network state *and* its state in all prior networks and
   adjust the statistics accordingly.

   In this state, .discord contains the changes in the network from
   the immediately prior time step, and is equivalent to the subset of
   .lasttoggle for which the last toggle time equals to the current
   time.

   .ticktock=TRUE:
   
   Changes to the network are assumed to be occurring as a part of the
   TERGM process. This means that a toggle that changes a dyad
   relative to the start of the time step will result in .lasttoggle
   being updated with the new last toggle time for that dyad, and
   .discord will be used to "back up" its previous last toggle
   time. If the toggle is reverted within the same time step, the dyad
   will be deleted from .discord and .lasttoggle time will be
   reverted.
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

#define JUST_CHANGED(dur_inf, tail, head) (kh_get(DyadMapInt, (dur_inf)->discord, THKey((dur_inf)->discord,(tail),(head)))!=kh_none)

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
