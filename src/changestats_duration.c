#include "changestats_duration.h"

/*****************
 void d_edges_ageinterval

 This is essentially the edges statistic, which only counts dyads with "age"
 (time steps spent in the current state) in the interval [inputparams0,inputparams1).
*****************/
D_CHANGESTAT_FN(d_edges_ageinterval){
  int edgeflag, i;
  Vertex tail, head;
  int from = INPUT_PARAM[0], to = INPUT_PARAM[1];
  
  ZERO_ALL_CHANGESTATS(i);
  FOR_EACH_TOGGLE(i){
    int age = ElapsedTime(tail=TAIL(i),head=HEAD(i),nwp);
    // Only count if the age is in [from,to). ( to=0 ==> to=Inf )
    if(from<=age && (to==0 || age<to)){
      edgeflag = IS_OUTEDGE(tail, head);
      CHANGE_STAT[0] += edgeflag ? - 1 : 1;
    }
    TOGGLE_IF_MORE_TO_COME(i);
  }
  UNDO_PREVIOUS_TOGGLES(i);
}

D_CHANGESTAT_FN(d_edges_ageinterval_mon){
  int edgeflag, i;
  Vertex tail, head;
  int from = INPUT_PARAM[0], to = INPUT_PARAM[1];
    
  ZERO_ALL_CHANGESTATS(i);

  for(Edge k=1; k <= N_EDGES; k++){
    Vertex tail, head;
    FindithEdge(&tail, &head, k, nwp);
    int age = ElapsedTime(tail,head,nwp); // Every tie "starts out" at age 1.
    if(age+1 == from) CHANGE_STAT[0]++; // The tie "ages" into the interval.
    if(to!=0 && age+1 == to) CHANGE_STAT[0]--; // The tie "ages" out of the interval.
  }

  FOR_EACH_TOGGLE(i){
    int age = ElapsedTime(tail=TAIL(i),head=HEAD(i),nwp);
    // Only count if the age is in [from,to). ( to=0 ==> to=Inf )
    if(from<=age+1 && (to==0 || age+1<to)){
      edgeflag = IS_OUTEDGE(tail, head);
      CHANGE_STAT[0] += edgeflag ? - 1 : +1;
    }
  }
}

S_CHANGESTAT_FN(s_edges_ageinterval_mon){
  int from=INPUT_PARAM[0], to=INPUT_PARAM[1];
  CHANGE_STAT[0] = 0;
  for (Edge k=1; k <= N_EDGES; k++){
    Vertex tail, head;
    FindithEdge(&tail, &head, k, nwp);
    int age = ElapsedTime(tail,head,nwp);
    if(from<=age && (to==0 || age<to)) CHANGE_STAT[0]++;
  }
}

/*****************
 void d_edge_ages

 Sum of ages of all extant ties.

*****************/

D_CHANGESTAT_FN(d_edge_ages_mon){
  int edgeflag, i;
  Vertex tail, head;
  
  // Changes
  // 0->0 ->  0
  // 0->1 -> +1
  // 1->0 -> -e
  // 1->1 -> +1
  
  // Basic idea:
  // 1. Start with a +1 for every edge in y0.
  // 2. For every edge being added, add +1.
  // 3. For every edge being removed, add -e, and, also -1, to cancel Step 1. 
  CHANGE_STAT[0] = N_EDGES;
  FOR_EACH_TOGGLE(i){
    int age = ElapsedTime(tail=TAIL(i),head=HEAD(i),nwp);
    edgeflag = IS_OUTEDGE(tail, head);
    CHANGE_STAT[0] += edgeflag ? - age - 1 : +1;
  }
}

S_CHANGESTAT_FN(s_edge_ages_mon){
  CHANGE_STAT[0] = 0;
  for (Edge k=1; k <= N_EDGES; k++){
    Vertex tail, head;
    FindithEdge(&tail, &head, k, nwp);
    int age = ElapsedTime(tail,head,nwp);
    CHANGE_STAT[0] += age;
  }
}

/*****************
 void d_edgecov_ages

 Weighted sum of ages of all extant ties. This quantity changes when the clock
 advances, so it needs a t_??? function.

*****************/

D_CHANGESTAT_FN(d_edgecov_ages_mon){
  int noffset = BIPARTITE, nrow;
  if(noffset > 0){
    nrow = noffset;
  }else{
    nrow = INPUT_PARAM[0];
  }

  int i;

  CHANGE_STAT[0]=0;
  
  // Sum of weights.
  for (Edge k=1; k <= N_EDGES; k++){
    Vertex tail, head;
    FindithEdge(&tail, &head, k, nwp);
    double val = INPUT_ATTRIB[(head-1-noffset)*nrow+(tail-1)];
    CHANGE_STAT[0] += val;
  }

  FOR_EACH_TOGGLE(i){
    Vertex tail=TAIL(i),head=HEAD(i);
    double val = INPUT_ATTRIB[(head - 1 - noffset) * nrow + (tail - 1)];
    if(val!=0){
      int age = ElapsedTime(tail, head, nwp);
      int edgeflag = IS_OUTEDGE(tail, head);
      CHANGE_STAT[0] += edgeflag ? - age*val - val : + val;
    }
  }
}

S_CHANGESTAT_FN(s_edgecov_ages_mon){
  int noffset = BIPARTITE, nrow;
  if(noffset > 0){
    nrow = noffset;
  }else{
    nrow = INPUT_PARAM[0];
  }

  CHANGE_STAT[0] = 0;
  for (Edge k=1; k <= N_EDGES; k++){
    Vertex tail, head;
    FindithEdge(&tail, &head, k, nwp);
    double val = INPUT_ATTRIB[(head-1-noffset)*nrow+(tail-1)];
    int age = ElapsedTime(tail,head,nwp);
    CHANGE_STAT[0] += age*val;
  }
}

/*****************
 void d_mean_age

 Mean of ages of all extant ties. This quantity changes when the clock
 advances, so it needs a t_??? function.

 The mean_ages of an empty network is defined to be emptyval.

 *****************/

D_CHANGESTAT_FN(d_mean_age_mon){
  int i;
  
  double s0 = 0, s1 = 0; // Sum of age values of initial and final network.
  double zeroval = INPUT_PARAM[0]; // Empty network value.
  Edge e0, e1; // Number of edges in initial and final network.
  
  e0 = e1 = N_EDGES;
  
  for(Edge k=1; k <= e0; k++){
    Vertex tail, head;
    FindithEdge(&tail, &head, k, nwp);
    s0 += ElapsedTime(tail,head,nwp);
    s1 += ElapsedTime(tail,head,nwp) + 1;
  }
  
  FOR_EACH_TOGGLE(i){
    Vertex tail = tails[i], head = heads[i];
    int et = ElapsedTime(tail,head,nwp); 
    if(IS_OUTEDGE(tail, head)){
      s1 -= et + 1;
      e1--;
    }else{
      s1 += 1;
      e1++;
    }
  }
  
  CHANGE_STAT[0]=(e1==0?zeroval:s1/e1)-(e0==0?zeroval:s0/e0);
}

S_CHANGESTAT_FN(s_mean_age_mon){
  CHANGE_STAT[0] = 0;
  double zeroval = INPUT_PARAM[0];

  if(N_EDGES>0){
    for (Edge k=1; k <= N_EDGES; k++){
      Vertex tail, head;
      FindithEdge(&tail, &head, k, nwp);
      int age = ElapsedTime(tail,head,nwp);
      CHANGE_STAT[0] += age;
    }
    
    CHANGE_STAT[0] /= N_EDGES;
  }else{
    CHANGE_STAT[0] = zeroval;
  }
}

/*****************
 void d_edgecov_mean_age

 Weigted mean of ages of all extant ties. This quantity changes when the clock
 advances, so it needs a t_??? function.

 The edgecov_mean_ages of an empty network is defined to be emptyval.

 *****************/

D_CHANGESTAT_FN(d_edgecov_age_mon){
  int noffset = BIPARTITE, nrow;
  if(noffset > 0){
    nrow = noffset;
  }else{
    nrow = INPUT_PARAM[1];
  }

  int i;

  double s0 = 0, s1 = 0; // Sum of age values of initial and final network.
  double zeroval = INPUT_PARAM[0]; // Empty network value.
  double e0 = 0, e1 = 0; // Sum of edge weights in initial and final network.

  for(Edge k=1; k <= N_EDGES; k++){
    Vertex tail, head;
    FindithEdge(&tail, &head, k, nwp);
    double val = INPUT_ATTRIB[(head - 1 - noffset) * nrow + (tail - 1)];   
    if(val!=0){
      s0 += ElapsedTime(tail,head,nwp)*val;
      s1 += (ElapsedTime(tail,head,nwp) + 1)*val;
      e0 += val;
    }
  }
  
  e1 = e0;
  
  FOR_EACH_TOGGLE(i){
    Vertex tail = TAIL(i), head = HEAD(i);
    double val = INPUT_ATTRIB[(head - 1 - noffset) * nrow + (tail - 1)];   
    if(val!=0){
      int et = ElapsedTime(tail,head,nwp); 
      if(IS_OUTEDGE(tail, head)){
	s1 -= (et+1)*val;
	e1 -= val;
      }else{
	s1 += 1;
	e1 += val;
      }
    }
  }
  
  CHANGE_STAT[0]=(e1==0?zeroval:s1/e1)-(e0==0?zeroval:s0/e0);
}

S_CHANGESTAT_FN(s_edgecov_mean_age_mon){
  CHANGE_STAT[0] = 0;
  double zeroval = INPUT_PARAM[0], s=0, e=0;
  int noffset = BIPARTITE, nrow;
  if(noffset > 0){
    nrow = noffset;
  }else{
    nrow = INPUT_PARAM[1];
  }

  for (Edge k=1; k <= N_EDGES; k++){
    Vertex tail, head;
    FindithEdge(&tail, &head, k, nwp);
    double val = INPUT_ATTRIB[(head - 1 - noffset) * nrow + (tail - 1)];
    if(val!=0){
      int age = ElapsedTime(tail,head,nwp);
      s += age * val;
      e += val;
    }
  }
   
  if(e!=0){
    CHANGE_STAT[0] = s/e;
  }else{
    CHANGE_STAT[0] = zeroval;
  }
}
