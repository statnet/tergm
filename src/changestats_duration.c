/*  File src/changestats_duration.c in package tergm, part of the
 *  Statnet suite of packages for network analysis, https://statnet.org .
 *
 *  This software is distributed under the GPL-3 license.  It is free,
 *  open source, and has the attribution requirements (GPL Section 7) at
 *  https://statnet.org/attribution .
 *
 *  Copyright 2008-2024 Statnet Commons
 */
#include "changestats_duration.h"

#define CSD_TRANSFORM_ET(et)                    \
  double ett=0;                                \
  double ett1=1;                        \
  switch(transform) {                        \
    case 0: ett = et; ett1 = et+1; break;                \
    case 1: ett = log(et); ett1 = log(et+1); break;        \
    default: error("Unrecognized dyad age transformation code."); \
  }                                \
  (void) ett; (void) ett1; // Get rid of unused variable warnings, since either ett or ett1 may be unused.

/*****************
 edges_ageinterval

 This is essentially the edges statistic, which only counts dyads with "age"
 (time steps spent in the current state) in the interval [inputparams0,inputparams1).
*****************/

X_CHANGESTAT_FN(x_edges_ageinterval){
  ZERO_ALL_CHANGESTATS();

  if(type == TICK) {
    GET_AUX_STORAGE(StoreTimeAndLasttoggle, dur_inf);
  
    EXEC_THROUGH_NET_EDGES_PRE(tail, head, edge_var, {
      int age = ElapsedTime(tail,head,dur_inf) + 1; // lasttoggle auxiliary has *not* yet been updated via TICK
      for(unsigned int j=0; j<N_CHANGE_STATS; j++) {
        unsigned int from = INPUT_PARAM[j*2];
        unsigned int to = INPUT_PARAM[j*2+1];
        if(age+1 == from) CHANGE_STAT[j]++; // The tie "ages" into the interval.
        if(to!=0 && age+1 == to) CHANGE_STAT[j]--; // The tie "ages" out of the interval.
      }        
    });    
  }
}

C_CHANGESTAT_FN(c_edges_ageinterval){
  GET_AUX_STORAGE(StoreTimeAndLasttoggle, dur_inf);

  int age = ElapsedTimeToggle(tail,head,dur_inf,edgestate);
  // Only count if the age is in [from,to). ( to=0 ==> to=Inf )

  for(unsigned int j=0; j<N_CHANGE_STATS; j++){
    unsigned int from = INPUT_PARAM[j*2], to = INPUT_PARAM[j*2+1];
    if(edgestate){ // If already an edge, we are dissolving.
      if(from<=age+1 && (to==0 || age+1<to)) CHANGE_STAT[j]--; // Statistic only changes if it's in the interval.
    }else{ // If not already an edge, we are forming.
      if(from<=age+1 && (to==0 || age+1<to)) CHANGE_STAT[j]++; // Statistic only changes if it's in the interval.
    }
  }
}

S_CHANGESTAT_FN(s_edges_ageinterval){
  GET_AUX_STORAGE(StoreTimeAndLasttoggle, dur_inf);
  
  ZERO_ALL_CHANGESTATS(i);
  EXEC_THROUGH_NET_EDGES_PRE(tail, head, edge_var, {
    int age = ElapsedTime(tail,head,dur_inf) + 1;
    for(unsigned int j=0; j<N_CHANGE_STATS; j++){
      unsigned int from = INPUT_PARAM[j*2];
      unsigned int to = INPUT_PARAM[j*2+1];
      if(from<=age && (to==0 || age<to)) CHANGE_STAT[j]++;
    }
  });
}

/*****************
 edge_ages

 Sum of ages of all extant ties.

*****************/

X_CHANGESTAT_FN(x_edge_ages){ 
  ZERO_ALL_CHANGESTATS();
  if(type == TICK) {
    CHANGE_STAT[0] = N_EDGES;
  }
}


C_CHANGESTAT_FN(c_edge_ages){
  GET_AUX_STORAGE(StoreTimeAndLasttoggle, dur_inf);
  
  int age = ElapsedTimeToggle(tail,head,dur_inf,edgestate);

  CHANGE_STAT[0] += edgestate ? - age - 1 : age + 1;
}

S_CHANGESTAT_FN(s_edge_ages){
  GET_AUX_STORAGE(StoreTimeAndLasttoggle, dur_inf);
  
  CHANGE_STAT[0] = 0;
  EXEC_THROUGH_NET_EDGES_PRE(tail, head, edge_var, {
    int age = ElapsedTime(tail,head,dur_inf) + 1;
    CHANGE_STAT[0] += age;
  });
}

/*****************
 edgecov_ages

 Weighted sum of ages of all extant ties.

*****************/

X_CHANGESTAT_FN(x_edgecov_ages){ 
  ZERO_ALL_CHANGESTATS();
  if(type == TICK) {
    int noffset = BIPARTITE, nrow;
    if(noffset > 0){
      nrow = noffset;
    }else{
      nrow = INPUT_PARAM[0];
    }
    
    // Sum of weights.
    EXEC_THROUGH_NET_EDGES_PRE(tail, head, edge_var, {
      double val = INPUT_ATTRIB[(head-1-noffset)*nrow+(tail-1)];
      CHANGE_STAT[0] += val;
    });
  }
}

C_CHANGESTAT_FN(c_edgecov_ages){
  GET_AUX_STORAGE(StoreTimeAndLasttoggle, dur_inf);
  
  int noffset = BIPARTITE, nrow;
  if(noffset > 0){
    nrow = noffset;
  }else{
    nrow = INPUT_PARAM[0];
  }

  double val = INPUT_ATTRIB[(head - 1 - noffset) * nrow + (tail - 1)];
  if(val!=0){
    int age = ElapsedTimeToggle(tail, head, dur_inf,edgestate);

    CHANGE_STAT[0] += edgestate ? - age*val - val : age*val + val;
  }
}

S_CHANGESTAT_FN(s_edgecov_ages){
  GET_AUX_STORAGE(StoreTimeAndLasttoggle, dur_inf);
  
  int noffset = BIPARTITE, nrow;
  if(noffset > 0){
    nrow = noffset;
  }else{
    nrow = INPUT_PARAM[0];
  }

  CHANGE_STAT[0] = 0;
  EXEC_THROUGH_NET_EDGES_PRE(tail, head, edge_var, {
    double val = INPUT_ATTRIB[(head-1-noffset)*nrow+(tail-1)];
    int age = ElapsedTime(tail,head,dur_inf) + 1;
    CHANGE_STAT[0] += age*val;
  });
}

/*****************
 mean_age

 Mean of (optionally log-) ages of all extant ties.

 The mean_ages of an empty network is defined to be emptyval.

 *****************/

typedef struct {
  double age; // sum of edge ages in current network
  double prop_age; // sum of edge ages in proposed network
} mean_age_storage;

I_CHANGESTAT_FN(i_mean_age){
  ALLOC_STORAGE(1, mean_age_storage, sto);
  GET_AUX_STORAGE(StoreTimeAndLasttoggle, dur_inf);  
  int transform = INPUT_PARAM[1];
  
  EXEC_THROUGH_NET_EDGES_PRE(tail, head, edge_var, {
    int et = ElapsedTime(tail,head,dur_inf);
    CSD_TRANSFORM_ET(et);
    sto->age += ett1;
  });
}

X_CHANGESTAT_FN(x_mean_age){
  ZERO_ALL_CHANGESTATS();
  GET_AUX_STORAGE(StoreTimeAndLasttoggle, dur_inf);  
  
  int transform = INPUT_PARAM[1];
  
  if(type == TICK) {
    GET_STORAGE(mean_age_storage, sto);
    
    if(transform == 0) {
      sto->age += N_EDGES;
      CHANGE_STAT[0] = N_EDGES ? 1 : 0;
    } else {
      double oldval = sto->age;
      sto->age = 0;
      EXEC_THROUGH_NET_EDGES_PRE(tail, head, edge_var, {
        int et = ElapsedTime(tail,head,dur_inf) + 1;
        CSD_TRANSFORM_ET(et);
        sto->age += ett1;
      });
      CHANGE_STAT[0] = N_EDGES ? (sto->age - oldval)/N_EDGES : 0;
    }
  }
}

void process_toggle_mean_age(Vertex tail, Vertex head, ModelTerm *mtp, Network *nwp, Rboolean edgestate, Rboolean write_changestats) {
  GET_STORAGE(mean_age_storage, sto);
  GET_AUX_STORAGE(StoreTimeAndLasttoggle, dur_inf);
    
  double s0 = sto->age, s1 = sto->age; // Sum of age values of initial and final network.
  double zeroval = INPUT_PARAM[0]; // Empty network value.
  int transform = INPUT_PARAM[1]; // Transformation code.
  Edge e0, e1; // Number of edges in initial and final network.
  
  e0 = e1 = N_EDGES;
  
  int change = edgestate ? -1 : 1;
  int et = ElapsedTimeToggle(tail,head,dur_inf,edgestate);
  CSD_TRANSFORM_ET(et);
  s1 += change*ett1;
  e1 += change;
  sto->prop_age = s1;

  if(write_changestats) {
    CHANGE_STAT[0]=(e1==0?zeroval:s1/e1)-(e0==0?zeroval:s0/e0);
  }
}

C_CHANGESTAT_FN(c_mean_age){
  process_toggle_mean_age(tail, head, mtp, nwp, edgestate, TRUE);
}

U_CHANGESTAT_FN(u_mean_age){
  process_toggle_mean_age(tail, head, mtp, nwp, edgestate, FALSE);

  GET_STORAGE(mean_age_storage, sto);
  sto->age = sto->prop_age;
}

S_CHANGESTAT_FN(s_mean_age){
  GET_AUX_STORAGE(StoreTimeAndLasttoggle, dur_inf);
  
  CHANGE_STAT[0] = 0;
  double zeroval = INPUT_PARAM[0];
  int transform = INPUT_PARAM[1]; // Transformation code.

  if(N_EDGES>0){
    EXEC_THROUGH_NET_EDGES_PRE(tail, head, edge_var, {
      int et = ElapsedTime(tail,head,dur_inf);
      CSD_TRANSFORM_ET(et);
      CHANGE_STAT[0] += ett1;
    });
    
    CHANGE_STAT[0] /= N_EDGES;
  }else{
    CHANGE_STAT[0] = zeroval;
  }
}

// mean age of half-ties by nodal attribute
typedef struct {
  int *nodecov;
  int *edges;
  double *ages;
  double *newages;
  double *emptyvals;
  int log;
} nodefactor_mean_age_storage;

I_CHANGESTAT_FN(i_nodefactor_mean_age) {
  ALLOC_STORAGE(1, nodefactor_mean_age_storage, sto);
  
  sto->nodecov = INTEGER(getListElement(mtp->R, "nodecov"));
  sto->log = asInteger(getListElement(mtp->R, "log"));
  sto->emptyvals = REAL(getListElement(mtp->R, "emptynwstats"));
  sto->edges = R_Calloc(N_CHANGE_STATS, int);
  sto->ages = R_Calloc(N_CHANGE_STATS, double);
  sto->newages = R_Calloc(N_CHANGE_STATS, double);

  // populate fields with initial edge set
  GET_AUX_STORAGE(StoreTimeAndLasttoggle, dur_inf);  
  
  int transform = sto->log;
  
  EXEC_THROUGH_NET_EDGES_PRE(tail, head, edge_var, {
    int et = ElapsedTime(tail,head,dur_inf);
    CSD_TRANSFORM_ET(et);
    int index = sto->nodecov[tail];
    if(index >= 0) {
      sto->ages[index] += ett1;
      sto->edges[index]++;
    }
    index = sto->nodecov[head];
    if(index >= 0) {
      sto->ages[index] += ett1;
      sto->edges[index]++;
    }
  });  
}

X_CHANGESTAT_FN(x_nodefactor_mean_age) {
  ZERO_ALL_CHANGESTATS();
  GET_AUX_STORAGE(StoreTimeAndLasttoggle, dur_inf);  
    
  if(type == TICK) {
    GET_STORAGE(nodefactor_mean_age_storage, sto);

    if(sto->log == 0) {
      for(int i = 0; i < N_CHANGE_STATS; i++) {
        sto->ages[i] += sto->edges[i];
        CHANGE_STAT[i] = sto->edges[i] ? 1 : 0;
      }
    } else {
      int transform = sto->log;

      double *oldages = R_Calloc(N_CHANGE_STATS, double);
      memcpy(oldages, sto->ages, N_CHANGE_STATS*sizeof(double));      
      memset(sto->ages, 0, N_CHANGE_STATS*sizeof(double));
      
      EXEC_THROUGH_NET_EDGES_PRE(tail, head, edge_var, {
        int et = ElapsedTime(tail, head, dur_inf) + 1;
        CSD_TRANSFORM_ET(et);
        
        int tailindex = sto->nodecov[tail];
        int headindex = sto->nodecov[head];
        if(tailindex >= 0) {
          sto->ages[tailindex] += ett1;
        }
        if(headindex >= 0) {
          sto->ages[headindex] += ett1;
        }
      });

      for(int i = 0; i < N_CHANGE_STATS; i++) {
        CHANGE_STAT[i] = sto->edges[i] ? (sto->ages[i] - oldages[i])/sto->edges[i] : 0;          
      }
      
      R_Free(oldages);
    }
  }  
}
    
void process_toggle_nodefactor_mean_age(Vertex tail, Vertex head, ModelTerm *mtp, Network *nwp, Rboolean edgestate, Rboolean write_changestats) {
  GET_STORAGE(nodefactor_mean_age_storage, sto);
  GET_AUX_STORAGE(StoreTimeAndLasttoggle, dur_inf);
  
  // NB: need to consider case tailtype == headtype
  int tailindex = sto->nodecov[tail];
  int headindex = sto->nodecov[head];
  
  if(tailindex == headindex && tailindex >= 0) {
    double s0 = sto->ages[tailindex], s1 = sto->ages[tailindex]; // Sum of age values of initial and final network.
    Edge e0 = sto->edges[tailindex], e1 = sto->edges[tailindex]; // Number of edges in initial and final network.
          
    int change = edgestate ? -2 : 2;
    int et = ElapsedTimeToggle(tail, head, dur_inf, edgestate);
    
    int transform = sto->log; // Transformation code.      
    CSD_TRANSFORM_ET(et);
    
    s1 += change*ett1;
    e1 += change;
    
    sto->newages[tailindex] = s1;
    
    if(write_changestats) {
      CHANGE_STAT[tailindex] = (e1 == 0 ? sto->emptyvals[tailindex] : s1/e1) - (e0 == 0 ? sto->emptyvals[tailindex] : s0/e0);      
    }
  } else {
    if(tailindex >= 0) {
      double s0 = sto->ages[tailindex], s1 = sto->ages[tailindex]; // Sum of age values of initial and final network.
      Edge e0 = sto->edges[tailindex], e1 = sto->edges[tailindex]; // Number of edges in initial and final network.
            
      int change = edgestate ? -1 : 1;
      int et = ElapsedTimeToggle(tail, head, dur_inf, edgestate);
      
      int transform = sto->log; // Transformation code.      
      CSD_TRANSFORM_ET(et);
      
      s1 += change*ett1;
      e1 += change;
      
      sto->newages[tailindex] = s1;
    
      if(write_changestats) {
        CHANGE_STAT[tailindex] = (e1 == 0 ? sto->emptyvals[tailindex] : s1/e1) - (e0 == 0 ? sto->emptyvals[tailindex] : s0/e0);
      }
    }
    if(headindex >= 0) {
      double s0 = sto->ages[headindex], s1 = sto->ages[headindex]; // Sum of age values of initial and final network.
      Edge e0 = sto->edges[headindex], e1 = sto->edges[headindex]; // Number of edges in initial and final network.
            
      int change = edgestate ? -1 : 1;
      int et = ElapsedTimeToggle(tail, head, dur_inf, edgestate);
      
      int transform = sto->log; // Transformation code.      
      CSD_TRANSFORM_ET(et);
      
      s1 += change*ett1;
      e1 += change;
      
      sto->newages[headindex] = s1;
    
      if(write_changestats) {
        CHANGE_STAT[headindex] = (e1 == 0 ? sto->emptyvals[headindex] : s1/e1) - (e0 == 0 ? sto->emptyvals[headindex] : s0/e0);      
      }
    }
  }
}

C_CHANGESTAT_FN(c_nodefactor_mean_age) {
  process_toggle_nodefactor_mean_age(tail, head, mtp, nwp, edgestate, TRUE);
}


U_CHANGESTAT_FN(u_nodefactor_mean_age) {
  process_toggle_nodefactor_mean_age(tail, head, mtp, nwp, edgestate, FALSE);

  GET_STORAGE(nodefactor_mean_age_storage, sto);

  // note the following is fine even if tailindex == headindex,
  // although it may be confusing to assign sto->ages[index] twice 
  // in that case
  int index = sto->nodecov[tail];
  if(index >= 0) {
    sto->ages[index] = sto->newages[index];
    sto->edges[index] += edgestate ? -1 : +1;
  }
  index = sto->nodecov[head];
  if(index >= 0) {
    sto->ages[index] = sto->newages[index];
    sto->edges[index] += edgestate ? -1 : +1;
  }  
}

F_CHANGESTAT_FN(f_nodefactor_mean_age) {
  GET_STORAGE(nodefactor_mean_age_storage, sto);

  R_Free(sto->edges);
  R_Free(sto->ages);
  R_Free(sto->newages);
}

S_CHANGESTAT_FN(s_nodefactor_mean_age) {  
  GET_AUX_STORAGE(StoreTimeAndLasttoggle, dur_inf);

  int *nodecov = INTEGER(getListElement(mtp->R, "nodecov"));
  double *emptyvals = REAL(getListElement(mtp->R, "emptynwstats"));
  int transform = asInteger(getListElement(mtp->R, "log"));

  int *edges = R_Calloc(N_CHANGE_STATS, int);
  double *ages = R_Calloc(N_CHANGE_STATS, double);
  
  EXEC_THROUGH_NET_EDGES_PRE(tail, head, edge_var, {
    int et = ElapsedTime(tail,head,dur_inf);
    CSD_TRANSFORM_ET(et);
    
    int tailindex = nodecov[tail];
    int headindex = nodecov[head];

    if(tailindex >= 0) {
      ages[tailindex] += ett1;
      edges[tailindex]++;
    }
    if(headindex >= 0) {
      ages[headindex] += ett1;
      edges[headindex]++;        
    }
  });
  
  for(int i = 0; i < N_CHANGE_STATS; i++) {
    if(edges[i] > 0) {
      CHANGE_STAT[i] = ages[i]/edges[i];
    } else {
      CHANGE_STAT[i] = emptyvals[i];
    }
  }

  R_Free(edges);
  R_Free(ages);
}

// mean age of ties by mixing type

// mean age of half-ties by nodal attribute
typedef struct {
  int *nodecov;
  int *edges;
  int **indmat;
  double *ages;
  double *newages;
  double *emptyvals;
  int log;
} nodemix_mean_age_storage;

I_CHANGESTAT_FN(i_nodemix_mean_age) {
  ALLOC_STORAGE(1, nodemix_mean_age_storage, sto);
  
  sto->nodecov = INTEGER(getListElement(mtp->R, "nodecov"));
  sto->log = asInteger(getListElement(mtp->R, "log"));
  sto->emptyvals = REAL(getListElement(mtp->R, "emptynwstats"));
  sto->edges = R_Calloc(N_CHANGE_STATS, int);
  sto->ages = R_Calloc(N_CHANGE_STATS, double);
  sto->newages = R_Calloc(N_CHANGE_STATS, double);

  int nr = asInteger(getListElement(mtp->R, "nr"));
  int nc = asInteger(getListElement(mtp->R, "nc"));
  
  sto->indmat = R_Calloc(nr, int *);
  sto->indmat[0] = INTEGER(getListElement(mtp->R, "indmat"));
  for(int i = 1; i < nr; i++) {
    sto->indmat[i] = sto->indmat[i - 1] + nc;
  }
  
  // populate fields with initial edge set
  GET_AUX_STORAGE(StoreTimeAndLasttoggle, dur_inf);  
  
  int transform = sto->log;
  
  EXEC_THROUGH_NET_EDGES_PRE(tail, head, edge_var, {
    int et = ElapsedTime(tail,head,dur_inf);
    CSD_TRANSFORM_ET(et);
    int index = sto->indmat[sto->nodecov[tail]][sto->nodecov[head]];
    if(index >= 0) {
      sto->ages[index] += ett1;
      sto->edges[index]++;
    }
  });  
}

X_CHANGESTAT_FN(x_nodemix_mean_age) {
  ZERO_ALL_CHANGESTATS();
  GET_AUX_STORAGE(StoreTimeAndLasttoggle, dur_inf);  
    
  if(type == TICK) {
    GET_STORAGE(nodemix_mean_age_storage, sto);

    if(sto->log == 0) {
      for(int i = 0; i < N_CHANGE_STATS; i++) {
        sto->ages[i] += sto->edges[i];
        CHANGE_STAT[i] = sto->edges[i] ? 1 : 0;
      }
    } else {
      int transform = sto->log;

      double *oldages = R_Calloc(N_CHANGE_STATS, double);
      memcpy(oldages, sto->ages, N_CHANGE_STATS*sizeof(double));      
      memset(sto->ages, 0, N_CHANGE_STATS*sizeof(double));
      
      EXEC_THROUGH_NET_EDGES_PRE(tail, head, edge_var, {
        int et = ElapsedTime(tail, head, dur_inf) + 1;
        CSD_TRANSFORM_ET(et);
        
        int index = sto->indmat[sto->nodecov[tail]][sto->nodecov[head]];
        if(index >= 0) {
          sto->ages[index] += ett1;
        }
      });

      for(int i = 0; i < N_CHANGE_STATS; i++) {
        CHANGE_STAT[i] = sto->edges[i] ? (sto->ages[i] - oldages[i])/sto->edges[i] : 0;          
      }
      
      R_Free(oldages);
    }
  }  
}
    
void process_toggle_nodemix_mean_age(Vertex tail, Vertex head, ModelTerm *mtp, Network *nwp, Rboolean edgestate, Rboolean write_changestats) {
  GET_STORAGE(nodemix_mean_age_storage, sto);
  GET_AUX_STORAGE(StoreTimeAndLasttoggle, dur_inf);
  
  int index = sto->indmat[sto->nodecov[tail]][sto->nodecov[head]];
  if(index >= 0) {
    double s0 = sto->ages[index], s1 = sto->ages[index]; // Sum of age values of initial and final network.
    Edge e0 = sto->edges[index], e1 = sto->edges[index]; // Number of edges in initial and final network.
          
    int change = edgestate ? -1 : 1;
    int et = ElapsedTimeToggle(tail, head, dur_inf, edgestate);
    
    int transform = sto->log; // Transformation code.      
    CSD_TRANSFORM_ET(et);
    
    s1 += change*ett1;
    e1 += change;
    
    sto->newages[index] = s1;
  
    if(write_changestats) {
      CHANGE_STAT[index] = (e1 == 0 ? sto->emptyvals[index] : s1/e1) - (e0 == 0 ? sto->emptyvals[index] : s0/e0);
    }
  }  
}

C_CHANGESTAT_FN(c_nodemix_mean_age) {
  process_toggle_nodemix_mean_age(tail, head, mtp, nwp, edgestate, TRUE);
}

U_CHANGESTAT_FN(u_nodemix_mean_age) {
  process_toggle_nodemix_mean_age(tail, head, mtp, nwp, edgestate, FALSE);

  GET_STORAGE(nodemix_mean_age_storage, sto);

  int index = sto->indmat[sto->nodecov[tail]][sto->nodecov[head]];
  if(index >= 0) {
    sto->ages[index] = sto->newages[index];
    sto->edges[index] += edgestate ? -1 : +1;
  }
}

F_CHANGESTAT_FN(f_nodemix_mean_age) {
  GET_STORAGE(nodemix_mean_age_storage, sto);

  R_Free(sto->edges);
  R_Free(sto->ages);
  R_Free(sto->newages);
  R_Free(sto->indmat);  
}

S_CHANGESTAT_FN(s_nodemix_mean_age) {
  GET_AUX_STORAGE(StoreTimeAndLasttoggle, dur_inf);

  int *nodecov = INTEGER(getListElement(mtp->R, "nodecov"));
  double *emptyvals = REAL(getListElement(mtp->R, "emptynwstats"));
  int transform = asInteger(getListElement(mtp->R, "log"));

  int nr = asInteger(getListElement(mtp->R, "nr"));
  int nc = asInteger(getListElement(mtp->R, "nc"));
  
  int **indmat = R_Calloc(nr, int *);
  indmat[0] = INTEGER(getListElement(mtp->R, "indmat"));
  for(int i = 1; i < nr; i++) {
    indmat[i] = indmat[i - 1] + nc;
  }

  int *edges = R_Calloc(N_CHANGE_STATS, int);
  double *ages = R_Calloc(N_CHANGE_STATS, double);
  
  EXEC_THROUGH_NET_EDGES_PRE(tail, head, edge_var, {
    int et = ElapsedTime(tail,head,dur_inf);
    CSD_TRANSFORM_ET(et);
    
    int index = indmat[nodecov[tail]][nodecov[head]];

    if(index >= 0) {
      ages[index] += ett1;
      edges[index]++;
    }
  });
  
  for(int i = 0; i < N_CHANGE_STATS; i++) {
    if(edges[i] > 0) {
      CHANGE_STAT[i] = ages[i]/edges[i];
    } else {
      CHANGE_STAT[i] = emptyvals[i];
    }
  }

  R_Free(indmat);
  R_Free(edges);
  R_Free(ages);
}

/*****************
 edgecov_mean_age

 Weighted mean of ages of all extant ties.

 The edgecov_mean_ages of an empty network is defined to be emptyval.

 *****************/
 
typedef struct {
  double agewts; // sum of age*wt over edges in current network
  double wts; // sum of wt over edges in current network
  double prop_agewts; // sum of age*wt over edges in proposed network
  double prop_wts; // sum of wt over edges in proposed network
} edgecov_mean_age_storage;

I_CHANGESTAT_FN(i_edgecov_mean_age) {
  ALLOC_STORAGE(1, edgecov_mean_age_storage, sto);
  GET_AUX_STORAGE(StoreTimeAndLasttoggle, dur_inf);  
  int transform = INPUT_PARAM[1];
  
  int noffset = BIPARTITE, nrow;
  if(noffset > 0){
    nrow = noffset;
  }else{
    nrow = INPUT_PARAM[2];
  }
  
  EXEC_THROUGH_NET_EDGES_PRE(tail, head, edge_var, {
    double val = INPUT_ATTRIB[(head - 1 - noffset) * nrow + (tail - 1)];   
    if(val!=0){   
      int et = ElapsedTime(tail,head,dur_inf);
      CSD_TRANSFORM_ET(et);
      sto->agewts += ett1*val;
      sto->wts += val;
    }
  });
}

X_CHANGESTAT_FN(x_edgecov_mean_age) {
  ZERO_ALL_CHANGESTATS();
  if(type == TICK) {
    GET_AUX_STORAGE(StoreTimeAndLasttoggle, dur_inf);  

    int noffset = BIPARTITE, nrow;
    if(noffset > 0){
      nrow = noffset;
    }else{
      nrow = INPUT_PARAM[2];
    }
  
    int transform = INPUT_PARAM[1];
    
    GET_STORAGE(edgecov_mean_age_storage, sto);
    
    if(sto->wts != 0) {
      if(transform == 0) {
        sto->agewts += sto->wts;
        CHANGE_STAT[0] = 1;
      } else {
        double oldval = sto->agewts;
        sto->agewts = 0;
        EXEC_THROUGH_NET_EDGES_PRE(tail, head, edge_var, {
          double val = INPUT_ATTRIB[(head - 1 - noffset) * nrow + (tail - 1)];   
          if(val!=0) {
            int et = ElapsedTime(tail,head,dur_inf) + 1;
            CSD_TRANSFORM_ET(et);
            sto->agewts += ett1*val;
          }
        });
        CHANGE_STAT[0] = (sto->agewts - oldval)/sto->wts;
      }
    }
  }
}

void process_toggle_edgecov_mean_age(Vertex tail, Vertex head, ModelTerm *mtp, Network *nwp, Rboolean edgestate, Rboolean write_changestats) {
  GET_AUX_STORAGE(StoreTimeAndLasttoggle, dur_inf);
  
  int noffset = BIPARTITE, nrow;
  if(noffset > 0){
    nrow = noffset;
  }else{
    nrow = INPUT_PARAM[2];
  }

  GET_STORAGE(edgecov_mean_age_storage, sto);

  double s0 = sto->agewts, s1 = sto->agewts; // Sum of age values times weights in initial and final network.
  double zeroval = INPUT_PARAM[0]; // Empty network value.
  int transform = INPUT_PARAM[1]; // Transformation code.
  double e0 = sto->wts, e1 = sto->wts; // Sum of edge weights in initial and final network.
  
  double val = INPUT_ATTRIB[(head - 1 - noffset) * nrow + (tail - 1)];   
  if(val!=0){
    if(edgestate){
      int et = ElapsedTimeToggle(tail,head,dur_inf,edgestate);
      CSD_TRANSFORM_ET(et);
      s1 -= ett1*val;
      e1 -= val;
    }else{
      int et = ElapsedTimeToggle(tail,head,dur_inf,edgestate);
      CSD_TRANSFORM_ET(et);
      s1 += ett1*val;
      e1 += val;
    }
  }
  
  if(write_changestats) {
    CHANGE_STAT[0]=(e1==0?zeroval:s1/e1)-(e0==0?zeroval:s0/e0);
  }
  
  sto->prop_agewts = s1;
  sto->prop_wts = e1;
}

C_CHANGESTAT_FN(c_edgecov_mean_age) {
  process_toggle_edgecov_mean_age(tail, head, mtp, nwp, edgestate, TRUE);
}

U_CHANGESTAT_FN(u_edgecov_mean_age){
  process_toggle_edgecov_mean_age(tail, head, mtp, nwp, edgestate, FALSE);

  GET_STORAGE(edgecov_mean_age_storage, sto);
  sto->agewts = sto->prop_agewts;
  sto->wts = sto->prop_wts;
}

S_CHANGESTAT_FN(s_edgecov_mean_age){
  GET_AUX_STORAGE(StoreTimeAndLasttoggle, dur_inf);
  
  CHANGE_STAT[0] = 0;
  double zeroval = INPUT_PARAM[0], s=0, e=0;
  int transform = INPUT_PARAM[1]; // Transformation code.
  int noffset = BIPARTITE, nrow;
  if(noffset > 0){
    nrow = noffset;
  }else{
    nrow = INPUT_PARAM[2];
  }

  EXEC_THROUGH_NET_EDGES_PRE(tail, head, edge_var, {
    double val = INPUT_ATTRIB[(head - 1 - noffset) * nrow + (tail - 1)];
    if(val!=0){
      int et = ElapsedTime(tail,head,dur_inf);    
      CSD_TRANSFORM_ET(et);
      s += ett1 * val;
      e += val;
    }
  });
   
  if(e!=0){
    CHANGE_STAT[0] = s/e;
  }else{
    CHANGE_STAT[0] = zeroval;
  }
}

/*****************
 degree_mean_age

 Mean of ages of all extant ties with a particular degree.

 The degree_mean_age of a network with no actors with degree of interest is defined to be emptyval.

 *****************/

typedef struct {
  double *ages;
  int *counts;
  double *prop_ages;
  int *prop_counts;
} degree_mean_age_storage;

I_CHANGESTAT_FN(i_degree_mean_age){
  ALLOC_STORAGE(1, degree_mean_age_storage, sto);
  
  sto->ages = R_Calloc(N_CHANGE_STATS, double);
  sto->counts = R_Calloc(N_CHANGE_STATS, int);
  
  sto->prop_ages = R_Calloc(N_CHANGE_STATS, double);
  sto->prop_counts = R_Calloc(N_CHANGE_STATS, int);
  
  GET_AUX_STORAGE(StoreTimeAndLasttoggle, dur_inf);
    
  Vertex *id=IN_DEG, *od=OUT_DEG;
  int transform = INPUT_PARAM[1]; // Transformation code.
    
  for(unsigned int j = 0; j < N_CHANGE_STATS; j++){
    double s0 = 0;
    Edge e0 = 0;
  
    Vertex deg = INPUT_PARAM[j+2];
      
    EXEC_THROUGH_NET_EDGES_PRE(tail, head, edge_var, {
      unsigned int w = (od[tail]+id[tail]==deg) + (od[head]+id[head]==deg);
        
      if(w){
        int et = ElapsedTime(tail,head,dur_inf);
        CSD_TRANSFORM_ET(et);
        s0 += ett1*w;
        e0+=w;
      }
    });
    
    sto->ages[j] = s0;
    sto->counts[j] = e0;
  }
}
 
X_CHANGESTAT_FN(x_degree_mean_age){
  ZERO_ALL_CHANGESTATS();
  if(type == TICK) {
    GET_AUX_STORAGE(StoreTimeAndLasttoggle, dur_inf);
    GET_STORAGE(degree_mean_age_storage, sto);
       
    Vertex *id=IN_DEG, *od=OUT_DEG;
    double zeroval = INPUT_PARAM[0];
    int transform = INPUT_PARAM[1]; // Transformation code.

    for(unsigned int j = 0; j < N_CHANGE_STATS; j++) {
      double s0, s1;
      int e0;
      if(transform == 0) { // do it the fast way
        s0 = sto->ages[j];
        e0 = sto->counts[j];
        
        s1 = s0 + e0;
      } else { // transform == 1 and we need to do it the old way
        s0 = 0;
        s1 = 0;
        e0 = 0;
      
        Vertex deg = INPUT_PARAM[j+2];
        
        EXEC_THROUGH_NET_EDGES_PRE(tail, head, edge_var, {
          unsigned int w = (od[tail]+id[tail]==deg) + (od[head]+id[head]==deg);
          
          if(w){
            int et = ElapsedTime(tail,head,dur_inf) + 1;
            CSD_TRANSFORM_ET(et);
            s0 += ett*w;
            s1 += ett1*w;
            e0+=w;
          }
        });      
      }
      
      CHANGE_STAT[j]=(e0==0?zeroval:s1/e0)-(e0==0?zeroval:s0/e0);
      
      sto->ages[j] = s1;      
    }
  }
}

void process_toggle_degree_mean_age(Vertex tail, Vertex head, ModelTerm *mtp, Network *nwp, Rboolean edgestate, Rboolean write_changestats) {
  GET_AUX_STORAGE(StoreTimeAndLasttoggle, dur_inf);
  
  Vertex *id=IN_DEG, *od=OUT_DEG;
  double zeroval = INPUT_PARAM[0];
  int transform = INPUT_PARAM[1]; // Transformation code.
  
  GET_STORAGE(degree_mean_age_storage, sto);
    
  for(unsigned int j = 0; j < N_CHANGE_STATS; j++){
    double s0 = sto->ages[j], s1 = sto->ages[j];
    Edge e0 = sto->counts[j], e1 = sto->counts[j];

    Vertex deg = INPUT_PARAM[j+2];
    
    int change = edgestate ? -1 : +1;
    int taildiff = (od[tail]+id[tail] + change == deg)-(od[tail]+id[tail] == deg);
    int headdiff = (od[head]+id[head] + change == deg)-(od[head]+id[head] == deg);

    Edge e;
    Vertex head1, tail1;
    
    switch(taildiff){
      case -1: // tail was previously counted, but is no longer
        STEP_THROUGH_OUTEDGES(tail, e, head1){
          int et = ElapsedTime(tail,head1,dur_inf);
          CSD_TRANSFORM_ET(et);
          s1 -= ett1;
          e1--;
        }
        STEP_THROUGH_INEDGES(tail, e, head1){
          int et = ElapsedTime(head1,tail,dur_inf);
          CSD_TRANSFORM_ET(et);
          s1 -= ett1;
          e1--;
        }
        break;
      // We don't need to do anything special for the focus dyad here:
      // if it's formed, then it wasn't counted in the first place;
      // if it's dissolved, then it will have been subtracted off by the previous two loops.
      
      case +1: // tail was previously not counted, but is now
        STEP_THROUGH_OUTEDGES(tail, e, head1){
          int et = ElapsedTime(tail,head1,dur_inf);
          CSD_TRANSFORM_ET(et);
          s1 += ett1;
          e1++;
        }
        STEP_THROUGH_INEDGES(tail, e, head1){
          int et = ElapsedTime(head1,tail,dur_inf);
          CSD_TRANSFORM_ET(et);
          s1 += ett1;
          e1++;
        }
        // Here, we need to handle the focus dyad:
        if(change==+1){// if it's formed, add to s1
          int et = ElapsedTimeToggle(tail,head,dur_inf,edgestate);
          CSD_TRANSFORM_ET(et);
          s1 += ett1;
          e1++;
        }else{// if it's dissolved, it had been counted in the previous two loops, and it should be subtracted
          int et = ElapsedTimeToggle(tail,head,dur_inf,edgestate);
          CSD_TRANSFORM_ET(et);
          s1 -= ett1;
          e1--;
        }
        break;
    }

    switch(headdiff){
      case -1: // head was previously counted, but is no longer
        STEP_THROUGH_OUTEDGES(head, e, tail1){
          int et = ElapsedTime(head,tail1,dur_inf);
          CSD_TRANSFORM_ET(et);
          s1 -= ett1;
          e1--;
        }
        STEP_THROUGH_INEDGES(head, e, tail1){
          int et = ElapsedTime(tail1,head,dur_inf);
          CSD_TRANSFORM_ET(et);
          s1 -= ett1;
          e1--;
        }
        break;
      // We don't need to do anything special for the focus dyad here:
      // if it's formed, then it wasn't counted in the first place;
      // if it's dissolved, then it will have been subtracted off by the previous two loops.

      case +1: // head was previously not counted, but is now
        STEP_THROUGH_OUTEDGES(head, e, tail1){
          int et = ElapsedTime(head,tail1,dur_inf);
          CSD_TRANSFORM_ET(et);
          s1 += ett1;
          e1++;
        }
        STEP_THROUGH_INEDGES(head, e, tail1){
          int et = ElapsedTime(tail1,head,dur_inf);
          CSD_TRANSFORM_ET(et);
          s1 += ett1;
          e1++;
        }
        // Here, we need to handle the focus dyad:
        if(change==+1){// if it's formed, add to s1
          int et = ElapsedTimeToggle(tail,head,dur_inf,edgestate);
          CSD_TRANSFORM_ET(et);
          s1 += ett1;
          e1++;
        }else{// if it's dissolved, it had been counted in the previous two loops, and it should be subtracted
          int et = ElapsedTimeToggle(tail,head,dur_inf,edgestate);
          CSD_TRANSFORM_ET(et);
          s1 -= ett1;
          e1--;
        }
        break;
    }
  
    sto->prop_ages[j] = s1;
    sto->prop_counts[j] = e1;
  
    if(write_changestats) {
      CHANGE_STAT[j]=(e1==0?zeroval:s1/e1)-(e0==0?zeroval:s0/e0);
    }
  }
}

C_CHANGESTAT_FN(c_degree_mean_age) {
  process_toggle_degree_mean_age(tail, head, mtp, nwp, edgestate, TRUE);
}

U_CHANGESTAT_FN(u_degree_mean_age){
  process_toggle_degree_mean_age(tail, head, mtp, nwp, edgestate, FALSE);

  GET_STORAGE(degree_mean_age_storage, sto);
  
  memcpy(sto->ages, sto->prop_ages, N_CHANGE_STATS*sizeof(double));
  memcpy(sto->counts, sto->prop_counts, N_CHANGE_STATS*sizeof(int));
}

F_CHANGESTAT_FN(f_degree_mean_age){
  GET_STORAGE(degree_mean_age_storage, sto);

  R_Free(sto->ages);
  R_Free(sto->counts);
  R_Free(sto->prop_ages);
  R_Free(sto->prop_counts);  
}

S_CHANGESTAT_FN(s_degree_mean_age){
  GET_AUX_STORAGE(StoreTimeAndLasttoggle, dur_inf);
  
  Vertex *id=IN_DEG, *od=OUT_DEG;
  double zeroval = INPUT_PARAM[0];
  int transform = INPUT_PARAM[1]; // Transformation code.

  ZERO_ALL_CHANGESTATS();

  for(unsigned int j = 0; j < N_CHANGE_STATS; j++){
    Vertex deg = INPUT_PARAM[j+2];
    Edge e=0;

    EXEC_THROUGH_NET_EDGES_PRE(tail, head, edge_var, {
      unsigned int w = (od[tail]+id[tail]==deg ? 1:0) + (od[head]+id[head]==deg ? 1:0);
      
      if(w){
        int et = ElapsedTime(tail,head,dur_inf);
        CSD_TRANSFORM_ET(et);
        CHANGE_STAT[j] += ett1*w;
        e+=w;
      }
    });
    
    if(e) CHANGE_STAT[j] /= e;
    else CHANGE_STAT[j] = zeroval;
  }
}

/*****************
 degree_by_attr_mean_age

 Mean of ages of all extant ties with a particular degree, by actor attribute.

 The degree_by_attr_mean_age of a network with no actors with degree of interest is defined to be emptyval.

 *****************/

typedef degree_mean_age_storage degree_by_attr_mean_age_storage;


I_CHANGESTAT_FN(i_degree_by_attr_mean_age){
  ALLOC_STORAGE(1, degree_by_attr_mean_age_storage, sto);
  
  sto->ages = R_Calloc(N_CHANGE_STATS, double);
  sto->counts = R_Calloc(N_CHANGE_STATS, int);
  
  sto->prop_ages = R_Calloc(N_CHANGE_STATS, double);
  sto->prop_counts = R_Calloc(N_CHANGE_STATS, int);
  
  GET_AUX_STORAGE(StoreTimeAndLasttoggle, dur_inf);
    
  Vertex *id=IN_DEG, *od=OUT_DEG;
  int transform = INPUT_PARAM[1]; // Transformation code.
    
  for(unsigned int j = 0; j < N_CHANGE_STATS; j++){
    double s0 = 0;
    Edge e0 = 0;
  
    Vertex deg = INPUT_PARAM[2*j+2];
    int testattr = INPUT_PARAM[2*j+3];
      
    EXEC_THROUGH_NET_EDGES_PRE(tail, head, edge_var, {
      Vertex taildeg = od[tail]+id[tail];
      Vertex headdeg = od[head]+id[head];
      int tailattr = INPUT_PARAM[2*N_CHANGE_STATS + tail + 1]; 
      int headattr = INPUT_PARAM[2*N_CHANGE_STATS + head + 1]; 
  
      unsigned int w = ((taildeg==deg && tailattr==testattr) ? 1 : 0) +
        ((headdeg==deg && headattr==testattr) ? 1 : 0);
        
      if(w){
        int et = ElapsedTime(tail,head,dur_inf);
        CSD_TRANSFORM_ET(et);
        s0 += ett1*w;
        e0+=w;
      }
    });
    
    sto->ages[j] = s0;
    sto->counts[j] = e0;
  }
}

X_CHANGESTAT_FN(x_degree_by_attr_mean_age){
  ZERO_ALL_CHANGESTATS();
  if(type == TICK) {
    GET_AUX_STORAGE(StoreTimeAndLasttoggle, dur_inf);
    GET_STORAGE(degree_by_attr_mean_age_storage, sto);
       
    Vertex *id=IN_DEG, *od=OUT_DEG;
    double zeroval = INPUT_PARAM[0];
    int transform = INPUT_PARAM[1]; // Transformation code.

    for(unsigned int j = 0; j < N_CHANGE_STATS; j++) {
      double s0, s1;
      int e0;
      if(transform == 0) { // do it the fast way
        s0 = sto->ages[j];
        e0 = sto->counts[j];
        
        s1 = s0 + e0;
      } else { // transform == 1 and we need to do it the old way
        s0 = 0;
        s1 = 0;
        e0 = 0;
      
        Vertex deg = INPUT_PARAM[2*j+2];
        int testattr = INPUT_PARAM[2*j+3];
        
        EXEC_THROUGH_NET_EDGES_PRE(tail, head, edge_var, {
          Vertex taildeg = od[tail]+id[tail];
          Vertex headdeg = od[head]+id[head];
          int tailattr = INPUT_PARAM[2*N_CHANGE_STATS + tail + 1]; 
          int headattr = INPUT_PARAM[2*N_CHANGE_STATS + head + 1]; 
  
          unsigned int w = ((taildeg==deg && tailattr==testattr) ? 1 : 0) +
            ((headdeg==deg && headattr==testattr) ? 1 : 0);
          
          if(w){
            int et = ElapsedTime(tail,head,dur_inf) + 1;
            CSD_TRANSFORM_ET(et);
            s0 += ett*w;
            s1 += ett1*w;
            e0+=w;
          }
        });      
      }
      
      CHANGE_STAT[j]=(e0==0?zeroval:s1/e0)-(e0==0?zeroval:s0/e0);
      
      sto->ages[j] = s1;      
    }
  }
}

void process_toggle_degree_by_attr_mean_age(Vertex tail, Vertex head, ModelTerm *mtp, Network *nwp, Rboolean edgestate, Rboolean write_changestats) {
  GET_AUX_STORAGE(StoreTimeAndLasttoggle, dur_inf);

  GET_STORAGE(degree_by_attr_mean_age_storage, sto);
      
  Vertex *id=IN_DEG, *od=OUT_DEG;
  double zeroval = INPUT_PARAM[0];
  int transform = INPUT_PARAM[1]; // Transformation code.
  
  for(unsigned int j = 0; j < N_CHANGE_STATS; j++){
    double s0 = sto->ages[j], s1 = sto->ages[j];
    Edge e0 = sto->counts[j], e1 = sto->counts[j];

    Vertex deg = INPUT_PARAM[2*j+2];
    int testattr = INPUT_PARAM[2*j+3];

    int tailattr = INPUT_PARAM[2*N_CHANGE_STATS + tail + 1]; 
    int headattr = INPUT_PARAM[2*N_CHANGE_STATS + head + 1]; 

    // If neither attribute matches, this toggle has no effect on the statistic.
    if(tailattr!=testattr && headattr!=testattr){
      sto->prop_ages[j] = sto->ages[j];
      sto->prop_counts[j] = sto->counts[j];
      continue; 
    }

    int change = edgestate ? -1 : +1;
    int taildiff = (od[tail]+id[tail] + change == deg)-(od[tail]+id[tail] == deg);
    int headdiff = (od[head]+id[head] + change == deg)-(od[head]+id[head] == deg);

    Edge e;
    Vertex head1, tail1;
    
    switch(taildiff * (tailattr==testattr)){ // If tailattr!=testattr, it'll look for case 0, i.e., do nothing.
      case -1: // tail was previously counted, but is no longer
        STEP_THROUGH_OUTEDGES(tail, e, head1){
          int et = ElapsedTime(tail,head1,dur_inf);
          CSD_TRANSFORM_ET(et);
          s1 -= ett1;
          e1--;
        }
        STEP_THROUGH_INEDGES(tail, e, head1){
          int et = ElapsedTime(head1,tail,dur_inf);
          CSD_TRANSFORM_ET(et);
          s1 -= ett1;
          e1--;
        }
        break;
      // We don't need to do anything special for the focus dyad here:
      // if it's formed, then it wasn't counted in the first place;
      // if it's dissolved, then it will have been subtracted off by the previous two loops.

      case +1: // tail was previously not counted, but is now
        STEP_THROUGH_OUTEDGES(tail, e, head1){
          int et = ElapsedTime(tail,head1,dur_inf);
          CSD_TRANSFORM_ET(et);
          s1 += ett1;
          e1++;
        }
        STEP_THROUGH_INEDGES(tail, e, head1){
          int et = ElapsedTime(head1,tail,dur_inf);
          CSD_TRANSFORM_ET(et);
          s1 += ett1;
          e1++;
        }
        // Here, we need to handle the focus dyad:
        if(change==+1){// if it's formed, add to s1
          int et = ElapsedTimeToggle(tail,head,dur_inf,edgestate);
          CSD_TRANSFORM_ET(et);
          s1 += ett1;
          e1++;
        }else{// if it's dissolved, it had been counted in the previous two loops, and it should be subtracted
          int et = ElapsedTimeToggle(tail,head,dur_inf,edgestate);
          CSD_TRANSFORM_ET(et);
          s1 -= ett1;
          e1--;
        }
        break;
    }

    switch(headdiff * (headattr==testattr)){ // If headattr!=testattr, it'll look for case 0, i.e., do nothing.
      case -1: // head was previously counted, but is no longer
        STEP_THROUGH_OUTEDGES(head, e, tail1){
          int et = ElapsedTime(head,tail1,dur_inf);
          CSD_TRANSFORM_ET(et);
          s1 -= ett1;
          e1--;
        }
        STEP_THROUGH_INEDGES(head, e, tail1){
          int et = ElapsedTime(tail1,head,dur_inf);
          CSD_TRANSFORM_ET(et);
          s1 -= ett1;
          e1--;
        }
        break;
      // We don't need to do anything special for the focus dyad here:
      // if it's formed, then it wasn't counted in the first place;
      // if it's dissolved, then it will have been subtracted off by the previous two loops.

      case +1: // head was previously not counted, but is now
        STEP_THROUGH_OUTEDGES(head, e, tail1){
          int et = ElapsedTime(head,tail1,dur_inf);
          CSD_TRANSFORM_ET(et);
          s1 += ett1;
          e1++;
        }
        STEP_THROUGH_INEDGES(head, e, tail1){
          int et = ElapsedTime(tail1,head,dur_inf);
          CSD_TRANSFORM_ET(et);
          s1 += ett1;
          e1++;
        }
        // Here, we need to handle the focus dyad:
        if(change==+1){// if it's formed, add to s1
          int et = ElapsedTimeToggle(tail,head,dur_inf,edgestate);
          CSD_TRANSFORM_ET(et);
          s1 += ett1;
          e1++;
        }else{// if it's dissolved, it had been counted in the previous two loops, and it should be subtracted
          int et = ElapsedTimeToggle(tail,head,dur_inf,edgestate);
          CSD_TRANSFORM_ET(et);
          s1 -= ett1;
          e1--;
        }
        break;
    }
  
    if(write_changestats) {
      CHANGE_STAT[j]=(e1==0?zeroval:s1/e1)-(e0==0?zeroval:s0/e0);
    }
    
    sto->prop_counts[j] = e1;
    sto->prop_ages[j] = s1;
  }
}

C_CHANGESTAT_FN(c_degree_by_attr_mean_age){
  process_toggle_degree_by_attr_mean_age(tail, head, mtp, nwp, edgestate, TRUE);
}

U_CHANGESTAT_FN(u_degree_by_attr_mean_age){
  process_toggle_degree_by_attr_mean_age(tail, head, mtp, nwp, edgestate, FALSE);

  GET_STORAGE(degree_by_attr_mean_age_storage, sto);
  
  memcpy(sto->ages, sto->prop_ages, N_CHANGE_STATS*sizeof(double));
  memcpy(sto->counts, sto->prop_counts, N_CHANGE_STATS*sizeof(int));
}

F_CHANGESTAT_FN(f_degree_by_attr_mean_age){
  GET_STORAGE(degree_by_attr_mean_age_storage, sto);

  R_Free(sto->ages);
  R_Free(sto->counts);
  R_Free(sto->prop_ages);
  R_Free(sto->prop_counts);  
}


S_CHANGESTAT_FN(s_degree_by_attr_mean_age){
  GET_AUX_STORAGE(StoreTimeAndLasttoggle, dur_inf);
  
  Vertex *id=IN_DEG, *od=OUT_DEG;
  double zeroval = INPUT_PARAM[0];
  int transform = INPUT_PARAM[1]; // Transformation code.
  
  ZERO_ALL_CHANGESTATS(i);

  for(unsigned int j = 0; j < N_CHANGE_STATS; j++){
    Edge e=0;

    EXEC_THROUGH_NET_EDGES_PRE(tail, head, edge_var, {
      Vertex taildeg = od[tail]+id[tail];
      Vertex headdeg = od[head]+id[head];
      int tailattr = INPUT_PARAM[2*N_CHANGE_STATS + tail + 1]; 
      int headattr = INPUT_PARAM[2*N_CHANGE_STATS + head + 1]; 
    
      Vertex deg = INPUT_PARAM[2*j+2];
      int testattr = INPUT_PARAM[2*j+3];

      unsigned int w = ((taildeg==deg && tailattr==testattr) ? 1 : 0) +
        ((headdeg==deg && headattr==testattr) ? 1 : 0);
      
      if(w){
        int et = ElapsedTime(tail,head,dur_inf);
        CSD_TRANSFORM_ET(et);
        CHANGE_STAT[j] += ett1*w;
        e+=w;
      }
    });
    
    if(e) CHANGE_STAT[j] /= e;
    else CHANGE_STAT[j] = zeroval;
  }
}


/*****************
 degrange_mean_age

 Mean of ages of all extant ties with a particular degree range.

 The degrange_mean_age of a network with no actors with degree of interest is defined to be emptyval.

 *****************/

// A macro indicating whether x is in [from,to)
#define FROM_TO(x, from, to) ((x)>=(from) && (x)<(to))

typedef degree_mean_age_storage degrange_mean_age_storage;

I_CHANGESTAT_FN(i_degrange_mean_age){
  ALLOC_STORAGE(1, degrange_mean_age_storage, sto);
  
  sto->ages = R_Calloc(N_CHANGE_STATS, double);
  sto->counts = R_Calloc(N_CHANGE_STATS, int);
  
  sto->prop_ages = R_Calloc(N_CHANGE_STATS, double);
  sto->prop_counts = R_Calloc(N_CHANGE_STATS, int);
  
  GET_AUX_STORAGE(StoreTimeAndLasttoggle, dur_inf);
    
  Vertex *id=IN_DEG, *od=OUT_DEG;
  int transform = INPUT_PARAM[1]; // Transformation code.
    
  for(unsigned int j = 0; j < N_CHANGE_STATS; j++){
    double s0 = 0;
    Edge e0 = 0;
  
    Vertex from = INPUT_PARAM[j*2+2], to = INPUT_PARAM[j*2+3];
      
    EXEC_THROUGH_NET_EDGES_PRE(tail, head, edge_var, {
      unsigned int w = FROM_TO(od[tail]+id[tail],from,to) + FROM_TO(od[head]+id[head],from,to);
        
      if(w){
        int et = ElapsedTime(tail,head,dur_inf);
        CSD_TRANSFORM_ET(et);
        s0 += ett1*w;
        e0+=w;
      }
    });
    
    sto->ages[j] = s0;
    sto->counts[j] = e0;
  }
}

X_CHANGESTAT_FN(x_degrange_mean_age){
  ZERO_ALL_CHANGESTATS();
  if(type == TICK) {
    GET_AUX_STORAGE(StoreTimeAndLasttoggle, dur_inf);
    GET_STORAGE(degrange_mean_age_storage, sto);
       
    Vertex *id=IN_DEG, *od=OUT_DEG;
    double zeroval = INPUT_PARAM[0];
    int transform = INPUT_PARAM[1]; // Transformation code.

    for(unsigned int j = 0; j < N_CHANGE_STATS; j++) {
      double s0, s1;
      int e0;
      if(transform == 0) { // do it the fast way
        s0 = sto->ages[j];
        e0 = sto->counts[j];
        
        s1 = s0 + e0;
      } else { // transform == 1 and we need to do it the old way
        s0 = 0;
        s1 = 0;
        e0 = 0;
      
        Vertex from = INPUT_PARAM[j*2+2], to = INPUT_PARAM[j*2+3];        
        
        EXEC_THROUGH_NET_EDGES_PRE(tail, head, edge_var, {
          unsigned int w = FROM_TO(od[tail]+id[tail],from,to) + FROM_TO(od[head]+id[head],from,to);
          
          if(w){
            int et = ElapsedTime(tail,head,dur_inf) + 1;
            CSD_TRANSFORM_ET(et);
            s0 += ett*w;
            s1 += ett1*w;
            e0+=w;
          }
        });      
      }
      
      CHANGE_STAT[j]=(e0==0?zeroval:s1/e0)-(e0==0?zeroval:s0/e0);
      
      sto->ages[j] = s1;      
    }
  }
}

void process_toggle_degrange_mean_age(Vertex tail, Vertex head, ModelTerm *mtp, Network *nwp, Rboolean edgestate, Rboolean write_changestats) {
  GET_AUX_STORAGE(StoreTimeAndLasttoggle, dur_inf);
  
  Vertex *id=IN_DEG, *od=OUT_DEG;
  double zeroval = INPUT_PARAM[0];
  int transform = INPUT_PARAM[1]; // Transformation code.
  
  GET_STORAGE(degrange_mean_age_storage, sto);  
  
  for(unsigned int j = 0; j < N_CHANGE_STATS; j++){
    double s0 = sto->ages[j], s1 = sto->ages[j];
    Edge e0 = sto->counts[j], e1 = sto->counts[j];

    Vertex from = INPUT_PARAM[j*2+2], to = INPUT_PARAM[j*2+3];

    int change = edgestate ? -1 : +1;
    // In the degree range case, it's possible to gain or lose a tie without entering or exiting a given degree range.
    unsigned int tailin1 = FROM_TO(od[tail]+id[tail] + change, from, to),
      tailin0 = FROM_TO(od[tail]+id[tail], from, to),
      headin1 = FROM_TO(od[head]+id[head] + change, from, to),
      headin0 = FROM_TO(od[head]+id[head], from, to);
    
    Edge e;
    Vertex head1, tail1;

    if(tailin0 && !tailin1){ // tail was previously counted, but is no longer
      STEP_THROUGH_OUTEDGES(tail, e, head1){
        int et = ElapsedTime(tail,head1,dur_inf);
        CSD_TRANSFORM_ET(et);
        s1 -= ett1;
        e1--;
      }
      STEP_THROUGH_INEDGES(tail, e, head1){
        int et = ElapsedTime(head1,tail,dur_inf);
        CSD_TRANSFORM_ET(et);
        s1 -= ett1;
        e1--;
      }
      // We don't need to do anything special for the focus dyad here:
      // if it's formed, then it wasn't counted in the first place;
      // if it's dissolved, then it will have been subtracted off by the previous two loops.
    }else if(!tailin0 && tailin1){ // tail was previously not counted, but is now
      STEP_THROUGH_OUTEDGES(tail, e, head1){
        int et = ElapsedTime(tail,head1,dur_inf);
        CSD_TRANSFORM_ET(et);
        s1 += ett1;
        e1++;
      }
      STEP_THROUGH_INEDGES(tail, e, head1){
        int et = ElapsedTime(head1,tail,dur_inf);
        CSD_TRANSFORM_ET(et);
        s1 += ett1;
        e1++;
      }
      // Here, we need to handle the focus dyad:
      if(change==+1){// if it's formed, add to s1
        int et = ElapsedTimeToggle(tail,head,dur_inf,edgestate);
        CSD_TRANSFORM_ET(et);
        s1 += ett1;
        e1++;
      }else{// if it's dissolved, it had been counted in the previous two loops, and it should be subtracted
        int et = ElapsedTimeToggle(tail,head,dur_inf,edgestate);
        CSD_TRANSFORM_ET(et);
        s1 -= ett1;
        e1--;
      }
    }else if(tailin0 && tailin1){ // tail was counted both times, but we need to handle the focus dyad
      if(change==+1){// if it's formed, add to s1
        int et = ElapsedTimeToggle(tail,head,dur_inf,edgestate);
        CSD_TRANSFORM_ET(et);
        s1 += ett1;
        e1++;
      }else{// if it's dissolved, it must be subtracted from s1
        int et = ElapsedTimeToggle(tail,head,dur_inf,edgestate);
        CSD_TRANSFORM_ET(et);
        s1 -= ett1;
        e1--;
      }
    }
    // If !tailin0 && !tailin1, then it made no difference.
    
    if(headin0 && !headin1){ // head was previously counted, but is no longer
      STEP_THROUGH_OUTEDGES(head, e, tail1){
        int et = ElapsedTime(head,tail1,dur_inf);
        CSD_TRANSFORM_ET(et);
        s1 -= ett1;
        e1--;
      }
      STEP_THROUGH_INEDGES(head, e, tail1){
        int et = ElapsedTime(tail1,head,dur_inf);
        CSD_TRANSFORM_ET(et);
        s1 -= ett1;
        e1--;
      }
      // We don't need to do anything special for the focus dyad here:
      // if it's formed, then it wasn't counted in the first place;
      // if it's dissolved, then it will have been subtracted off by the previous two loops.
    }else if(!headin0 && headin1){ // head was previously not counted, but is now
      STEP_THROUGH_OUTEDGES(head, e, tail1){
        int et = ElapsedTime(head,tail1,dur_inf);
        CSD_TRANSFORM_ET(et);
        s1 += ett1;
        e1++;
      }
      STEP_THROUGH_INEDGES(head, e, tail1){
        int et = ElapsedTime(tail1,head,dur_inf);
        CSD_TRANSFORM_ET(et);
        s1 += ett1;
        e1++;
      }
      // Here, we need to handle the focus dyad:
      if(change==+1){// if it's formed, add to s1
        int et = ElapsedTimeToggle(tail,head,dur_inf,edgestate);
        CSD_TRANSFORM_ET(et);
        s1 += ett1;
        e1++;
      }else{// if it's dissolved, it had been counted in the previous two loops, and it should be subtracted
        int et = ElapsedTimeToggle(tail,head,dur_inf,edgestate);
        CSD_TRANSFORM_ET(et);
        s1 -= ett1;
        e1--;
      }
    }else if(headin0 && headin1){ // tail was counted both times, but we need to handle the focus dyad
      if(change==+1){// if it's formed, add to s1
        int et = ElapsedTimeToggle(tail,head,dur_inf,edgestate);
        CSD_TRANSFORM_ET(et);
        s1 += ett1;
        e1++;
      }else{// if it's dissolved, it must be subtracted from s1
        int et = ElapsedTimeToggle(tail,head,dur_inf,edgestate);
        CSD_TRANSFORM_ET(et);
        s1 -= ett1;
        e1--;
      }
    }
    // If !headin0 && !headin1, then it made no difference.
    
    if(write_changestats) {
      CHANGE_STAT[j]=(e1==0?zeroval:s1/e1)-(e0==0?zeroval:s0/e0);
    }
    
    sto->prop_ages[j] = s1;
    sto->prop_counts[j] = e1;
  }
}

C_CHANGESTAT_FN(c_degrange_mean_age) {
  process_toggle_degrange_mean_age(tail, head, mtp, nwp, edgestate, TRUE);
}

U_CHANGESTAT_FN(u_degrange_mean_age){
  process_toggle_degrange_mean_age(tail, head, mtp, nwp, edgestate, FALSE);

  GET_STORAGE(degrange_mean_age_storage, sto);

  memcpy(sto->ages, sto->prop_ages, N_CHANGE_STATS*sizeof(double));
  memcpy(sto->counts, sto->prop_counts, N_CHANGE_STATS*sizeof(int));
}

F_CHANGESTAT_FN(f_degrange_mean_age){
  GET_STORAGE(degrange_mean_age_storage, sto);

  R_Free(sto->ages);
  R_Free(sto->counts);
  R_Free(sto->prop_ages);
  R_Free(sto->prop_counts);  
}


S_CHANGESTAT_FN(s_degrange_mean_age){
  GET_AUX_STORAGE(StoreTimeAndLasttoggle, dur_inf);
  
  Vertex *id=IN_DEG, *od=OUT_DEG;
  double zeroval = INPUT_PARAM[0];
  int transform = INPUT_PARAM[1]; // Transformation code.
  
  ZERO_ALL_CHANGESTATS(i);

  for(unsigned int j = 0; j < N_CHANGE_STATS; j++){
    Vertex from = INPUT_PARAM[j*2+2], to = INPUT_PARAM[j*2+3];
    Edge e=0;

    EXEC_THROUGH_NET_EDGES_PRE(tail, head, edge_var, {
      unsigned int w = FROM_TO(od[tail]+id[tail],from,to) + FROM_TO(od[head]+id[head],from,to);
      
      if(w){
        int et = ElapsedTime(tail,head,dur_inf);    
        CSD_TRANSFORM_ET(et);
        CHANGE_STAT[j] += ett1*w;
        e+=w;
      }
    });
    
    if(e) CHANGE_STAT[j] /= e;
    else CHANGE_STAT[j] = zeroval;
  }
}

/*****************
 degrange_by_attr_mean_age

 Mean of ages of all extant ties with a particular degree, by actor attribute.

 The degrange_by_attr_mean_age of a network with no actors with degree of interest is defined to be emptyval.

 *****************/

typedef degree_mean_age_storage degrange_by_attr_mean_age_storage;

I_CHANGESTAT_FN(i_degrange_by_attr_mean_age){
  ALLOC_STORAGE(1, degrange_by_attr_mean_age_storage, sto);
  
  sto->ages = R_Calloc(N_CHANGE_STATS, double);
  sto->counts = R_Calloc(N_CHANGE_STATS, int);
  
  sto->prop_ages = R_Calloc(N_CHANGE_STATS, double);
  sto->prop_counts = R_Calloc(N_CHANGE_STATS, int);
  
  GET_AUX_STORAGE(StoreTimeAndLasttoggle, dur_inf);
    
  Vertex *id=IN_DEG, *od=OUT_DEG;
  int transform = INPUT_PARAM[1]; // Transformation code.
    
  for(unsigned int j = 0; j < N_CHANGE_STATS; j++){
    double s0 = 0;
    Edge e0 = 0;
  
    Vertex from = INPUT_PARAM[3*j+2], to = INPUT_PARAM[3*j+3];
    int testattr = INPUT_PARAM[3*j+4];

      
    EXEC_THROUGH_NET_EDGES_PRE(tail, head, edge_var, {
      Vertex taildeg = od[tail]+id[tail];
      Vertex headdeg = od[head]+id[head];
      int tailattr = INPUT_PARAM[3*N_CHANGE_STATS + tail + 1]; 
      int headattr = INPUT_PARAM[3*N_CHANGE_STATS + head + 1]; 
   
      unsigned int w = (FROM_TO(taildeg, from, to) && tailattr==testattr) +
        (FROM_TO(headdeg, from, to) && headattr==testattr);
        
      if(w){
        int et = ElapsedTime(tail,head,dur_inf);
        CSD_TRANSFORM_ET(et);
        s0 += ett1*w;
        e0+=w;
      }
    });
    
    sto->ages[j] = s0;
    sto->counts[j] = e0;
  }
}

X_CHANGESTAT_FN(x_degrange_by_attr_mean_age){
  ZERO_ALL_CHANGESTATS();
  if(type == TICK) {
    GET_AUX_STORAGE(StoreTimeAndLasttoggle, dur_inf);
    GET_STORAGE(degrange_by_attr_mean_age_storage, sto);
       
    Vertex *id=IN_DEG, *od=OUT_DEG;
    double zeroval = INPUT_PARAM[0];
    int transform = INPUT_PARAM[1]; // Transformation code.

    for(unsigned int j = 0; j < N_CHANGE_STATS; j++) {
      double s0, s1;
      int e0;
      if(transform == 0) { // do it the fast way
        s0 = sto->ages[j];
        e0 = sto->counts[j];
        
        s1 = s0 + e0;
      } else { // transform == 1 and we need to do it the old way
        s0 = 0;
        s1 = 0;
        e0 = 0;
      
        Vertex from = INPUT_PARAM[3*j+2], to = INPUT_PARAM[3*j+3];
        int testattr = INPUT_PARAM[3*j+4];
        
        EXEC_THROUGH_NET_EDGES_PRE(tail, head, edge_var, {
          Vertex taildeg = od[tail]+id[tail];
          Vertex headdeg = od[head]+id[head];
          int tailattr = INPUT_PARAM[3*N_CHANGE_STATS + tail + 1]; 
          int headattr = INPUT_PARAM[3*N_CHANGE_STATS + head + 1]; 
   
          unsigned int w = (FROM_TO(taildeg, from, to) && tailattr==testattr) +
            (FROM_TO(headdeg, from, to) && headattr==testattr);
          
          if(w){
            int et = ElapsedTime(tail,head,dur_inf) + 1;
            CSD_TRANSFORM_ET(et);
            s0 += ett*w;
            s1 += ett1*w;
            e0+=w;
          }
        });      
      }
      
      CHANGE_STAT[j]=(e0==0?zeroval:s1/e0)-(e0==0?zeroval:s0/e0);
      
      sto->ages[j] = s1;      
    }
  }
}

void process_toggle_degrange_by_attr_mean_age(Vertex tail, Vertex head, ModelTerm *mtp, Network *nwp, Rboolean edgestate, Rboolean write_changestats) {
  GET_AUX_STORAGE(StoreTimeAndLasttoggle, dur_inf);
  
  Vertex *id=IN_DEG, *od=OUT_DEG;
  double zeroval = INPUT_PARAM[0];
  int transform = INPUT_PARAM[1]; // Transformation code.
  
  GET_STORAGE(degrange_by_attr_mean_age_storage, sto);
    
  for(unsigned int j = 0; j < N_CHANGE_STATS; j++){
    double s0 = sto->ages[j], s1 = sto->ages[j];
    Edge e0 = sto->counts[j], e1 = sto->counts[j];

    Vertex from = INPUT_PARAM[3*j+2], to = INPUT_PARAM[3*j+3];
    int testattr = INPUT_PARAM[3*j+4];
    
    int tailattr = INPUT_PARAM[3*N_CHANGE_STATS + tail + 1]; 
    int headattr = INPUT_PARAM[3*N_CHANGE_STATS + head + 1]; 

    // If neither attribute matches, this toggle has no effect on the statistic.
    if(tailattr!=testattr && headattr!=testattr){
      sto->prop_ages[j] = sto->ages[j];
      sto->prop_counts[j] = sto->counts[j];
      continue; 
    }

    int change = edgestate ? -1 : +1;
    // In the degree range case, it's possible to gain or lose a tie without entering or exiting a given degree range.
    unsigned int tailin1 = FROM_TO(od[tail]+id[tail] + change, from, to),
      tailin0 = FROM_TO(od[tail]+id[tail], from, to),
      headin1 = FROM_TO(od[head]+id[head] + change, from, to),
      headin0 = FROM_TO(od[head]+id[head], from, to);

    Edge e;
    Vertex head1, tail1;
    
    if(tailattr==testattr){
      if(tailin0 && !tailin1){ // tail was previously counted, but is no longer
        STEP_THROUGH_OUTEDGES(tail, e, head1){
          int et = ElapsedTime(tail,head1,dur_inf);
          CSD_TRANSFORM_ET(et);
          s1 -= ett1;
          e1--;
        }
        STEP_THROUGH_INEDGES(tail, e, head1){
          int et = ElapsedTime(head1,tail,dur_inf);
          CSD_TRANSFORM_ET(et);
          s1 -= ett1;
          e1--;
        }
        // We don't need to do anything special for the focus dyad here:
        // if it's formed, then it wasn't counted in the first place;
        // if it's dissolved, then it will have been subtracted off by the previous two loops.
      }else if(!tailin0 && tailin1){ // tail was previously not counted, but is now
        STEP_THROUGH_OUTEDGES(tail, e, head1){
          int et = ElapsedTime(tail,head1,dur_inf);
          CSD_TRANSFORM_ET(et);
          s1 += ett1;
          e1++;
        }
        STEP_THROUGH_INEDGES(tail, e, head1){
          int et = ElapsedTime(head1,tail,dur_inf);
          CSD_TRANSFORM_ET(et);
          s1 += ett1;
          e1++;
        }
        // Here, we need to handle the focus dyad:
        if(change==+1){// if it's formed, add to s1
          int et = ElapsedTimeToggle(tail,head,dur_inf,edgestate);
          CSD_TRANSFORM_ET(et);
          s1 += ett1;
          e1++;
        }else{// if it's dissolved, it had been counted in the previous two loops, and it should be subtracted
          int et = ElapsedTimeToggle(tail,head,dur_inf,edgestate);
          CSD_TRANSFORM_ET(et);
          s1 -= ett1;
          e1--;
        }
      }else if(tailin0 && tailin1){ // tail was counted both times, but we need to handle the focus dyad
        if(change==+1){// if it's formed, add to s1
          int et = ElapsedTimeToggle(tail,head,dur_inf,edgestate);
          CSD_TRANSFORM_ET(et);
          s1 += ett1;
          e1++;
        }else{// if it's dissolved, it must be subtracted from s1
          int et = ElapsedTimeToggle(tail,head,dur_inf,edgestate);
          CSD_TRANSFORM_ET(et);
          s1 -= ett1;
          e1--;
        }
      }
      // If !tailin0 && !tailin1, then it made no difference.
    }
    
    if(headattr==testattr){
      if(headin0 && !headin1){ // head was previously counted, but is no longer
        STEP_THROUGH_OUTEDGES(head, e, tail1){
          int et = ElapsedTime(head,tail1,dur_inf);
          CSD_TRANSFORM_ET(et);
          s1 -= ett1;
          e1--;
        }
        STEP_THROUGH_INEDGES(head, e, tail1){
          int et = ElapsedTime(tail1,head,dur_inf);
          CSD_TRANSFORM_ET(et);
          s1 -= ett1;
          e1--;
        }
        // We don't need to do anything special for the focus dyad here:
        // if it's formed, then it wasn't counted in the first place;
        // if it's dissolved, then it will have been subtracted off by the previous two loops.
      }else if(!headin0 && headin1){ // head was previously not counted, but is now
        STEP_THROUGH_OUTEDGES(head, e, tail1){
          int et = ElapsedTime(head,tail1,dur_inf);
          CSD_TRANSFORM_ET(et);
          s1 += ett1;
          e1++;
        }
        STEP_THROUGH_INEDGES(head, e, tail1){
          int et = ElapsedTime(tail1,head,dur_inf);
          CSD_TRANSFORM_ET(et);
          s1 += ett1;
          e1++;
        }
        // Here, we need to handle the focus dyad:
        if(change==+1){// if it's formed, add to s1
          int et = ElapsedTimeToggle(tail,head,dur_inf,edgestate);
          CSD_TRANSFORM_ET(et);
          s1 += ett1;
          e1++;
        }else{// if it's dissolved, it had been counted in the previous two loops, and it should be subtracted
          int et = ElapsedTimeToggle(tail,head,dur_inf,edgestate);
          CSD_TRANSFORM_ET(et);
          s1 -= ett1;
          e1--;
        }
      }else if(headin0 && headin1){ // tail was counted both times, but we need to handle the focus dyad
        if(change==+1){// if it's formed, add to s1
          int et = ElapsedTimeToggle(tail,head,dur_inf,edgestate);
          CSD_TRANSFORM_ET(et);
          s1 += ett1;
          e1++;
        }else{// if it's dissolved, it must be subtracted from s1
          int et = ElapsedTimeToggle(tail,head,dur_inf,edgestate);
          CSD_TRANSFORM_ET(et);
          s1 -= ett1;
          e1--;
        }
      }
      // If !headin0 && !headin1, then it made no difference.
    }
  
    if(write_changestats) {
      CHANGE_STAT[j]=(e1==0?zeroval:s1/e1)-(e0==0?zeroval:s0/e0);
    }
    
    sto->prop_ages[j] = s1;
    sto->prop_counts[j] = e1;
  }
}

C_CHANGESTAT_FN(c_degrange_by_attr_mean_age){
  process_toggle_degrange_by_attr_mean_age(tail, head, mtp, nwp, edgestate, TRUE);
}


U_CHANGESTAT_FN(u_degrange_by_attr_mean_age){
  process_toggle_degrange_by_attr_mean_age(tail, head, mtp, nwp, edgestate, FALSE);

  GET_STORAGE(degrange_by_attr_mean_age_storage, sto);

  memcpy(sto->ages, sto->prop_ages, N_CHANGE_STATS*sizeof(double));
  memcpy(sto->counts, sto->prop_counts, N_CHANGE_STATS*sizeof(int));
}

F_CHANGESTAT_FN(f_degrange_by_attr_mean_age){
  GET_STORAGE(degrange_by_attr_mean_age_storage, sto);

  R_Free(sto->ages);
  R_Free(sto->counts);
  R_Free(sto->prop_ages);
  R_Free(sto->prop_counts);  
}


S_CHANGESTAT_FN(s_degrange_by_attr_mean_age){
  GET_AUX_STORAGE(StoreTimeAndLasttoggle, dur_inf);
  
  Vertex *id=IN_DEG, *od=OUT_DEG;
  double zeroval = INPUT_PARAM[0];
  int transform = INPUT_PARAM[1]; // Transformation code.
  ZERO_ALL_CHANGESTATS(i);

  for(unsigned int j = 0; j < N_CHANGE_STATS; j++){
    Edge e=0;

    EXEC_THROUGH_NET_EDGES_PRE(tail, head, edge_var, {
      Vertex taildeg = od[tail]+id[tail];
      Vertex headdeg = od[head]+id[head];
      int tailattr = INPUT_PARAM[3*N_CHANGE_STATS + tail + 1]; 
      int headattr = INPUT_PARAM[3*N_CHANGE_STATS + head + 1]; 
    
      Vertex from = INPUT_PARAM[3*j+2];
      Vertex to = INPUT_PARAM[3*j+3];
      int testattr = INPUT_PARAM[3*j+4];

      unsigned int w = (FROM_TO(taildeg, from, to) && tailattr==testattr) +
        (FROM_TO(headdeg, from, to) && headattr==testattr);
      
      if(w){
        int et = ElapsedTime(tail,head,dur_inf);
        CSD_TRANSFORM_ET(et);
        CHANGE_STAT[j] += ett1*w;
        e+=w;
      }
    });
    
    if(e) CHANGE_STAT[j] /= e;
    else CHANGE_STAT[j] = zeroval;
  }
}

#undef FROM_TO
#undef CSD_TRANSFORM_ET
