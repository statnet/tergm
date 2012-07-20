#include "changestats_duration.h"

/*****************
 void d_edges_ageinterval

 This is essentially the edges statistic, which only counts dyads with "age"
 (time steps spent in the current state) in the interval [inputparams0,inputparams1).

 It changes whenever the clock advances, so it needs a t_??? function.
*****************/
D_CHANGESTAT_FN(d_edges_ageinterval){
  int edgeflag, i;
  Vertex tail, head;
  int from = mtp->inputparams[0], to = mtp->inputparams[1];
  
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
  int from = mtp->inputparams[0], to = mtp->inputparams[1];
  
  ZERO_ALL_CHANGESTATS(i);
  FOR_EACH_TOGGLE(i){
    int age = ElapsedTime(tail=TAIL(i),head=HEAD(i),nwp);
    // Only count if the age is in [from,to). ( to=0 ==> to=Inf )
    if(from<=age && (to==0 || age<to)){
      edgeflag = IS_OUTEDGE(tail, head);
      CHANGE_STAT[0] += edgeflag ? - 1 : 0;
    }
    TOGGLE_IF_MORE_TO_COME(i);
  }
  UNDO_PREVIOUS_TOGGLES(i);
}

T_CHANGESTAT_FN(t_edges_ageinterval_mon){
  int from=mtp->inputparams[0], to=mtp->inputparams[1];
  CHANGE_STAT[0] = 0;
  for (Edge k=1; k <= N_EDGES; k++){
    Vertex tail, head;
    FindithEdge(&tail, &head, k, nwp);
    int age = ElapsedTime(tail,head,nwp) + 1; // Every tie "starts out" at age 1.
    if(age == from) CHANGE_STAT[0]++; // The tie "ages" into the interval.
    if(to!=0 && age == to) CHANGE_STAT[0]--; // The tie "ages" out of the interval.
  }
}

S_CHANGESTAT_FN(s_edges_ageinterval_mon){
  int from=mtp->inputparams[0], to=mtp->inputparams[1];
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

 Sum of ages of all extant ties. This quantity changes when the clock
 advances, so it needs a t_??? function.

*****************/

D_CHANGESTAT_FN(d_edge_ages_mon){
  int edgeflag, i;
  Vertex tail, head;
  
  ZERO_ALL_CHANGESTATS(i);
  FOR_EACH_TOGGLE(i){
    int age = ElapsedTime(tail=TAIL(i),head=HEAD(i),nwp)+1;
    edgeflag = IS_OUTEDGE(tail, head);
    CHANGE_STAT[0] += edgeflag ? - age : 0;
    TOGGLE_IF_MORE_TO_COME(i);
  }
  UNDO_PREVIOUS_TOGGLES(i);
}

T_CHANGESTAT_FN(t_edge_ages_mon){
  CHANGE_STAT[0] = +N_EDGES; // Each extant edge's age increases by 1, so their sum increases by their number.
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
 void d_mean_age

 Mean of ages of all extant ties. This quantity changes when the clock
 advances, so it needs a t_??? function.

 The mean_ages of an empty network is defined to be 0.

 It's slower than edge_ages and is likely to have severe problems on very sparse models, but it's more numerically stable.



 *****************/

D_CHANGESTAT_FN(d_mean_age_mon){
  int i;
  /* double s0 = 0, // Sum of age values of dissolved edges and */
  /* 		s01 = 0; // Sum of age values of preserved edges. */
  /* Edge e0, e1; // Number of edges in initial and final network. */
  
  /* e0 = e1 = N_EDGES; */
  
  /* for(Edge k=1; k <= e0; k++){ */
  /* 	Vertex tail, head; */
  /* 	FindithEdge(&tail, &head, k, nwp); */
  /* 	s01 += ElapsedTime(tail,head,nwp)+1; */
  /* } */
  
  /* FOR_EACH_TOGGLE(i){ */
  /* 	Vertex tail = tails[i], head = heads[i]; */
  /* 	if(IS_OUTEDGE(tail, head)){ */
  /* 		int et = ElapsedTime(tail,head,nwp)+1;  */
  /* 		s0 += et; */
  /* 		s01 -= et; */
  /* 		e1--; */
  /* 	}else{ */
  /* 		e1++; */
  /* 	} */
  /* 	// Note that adding ties doesn't increase mean ages until the clock advances. */
  /* 	// We know that toggles are unique, since it's a monitoring function, so we don't have to toggle and toggle back. */
  /* 	//TOGGLE(tail,head); */
  /* } */
  
  //FOR_EACH_TOGGLE(i) TOGGLE(tails[(i)],heads[(i)]);
  
  //CHANGE_STAT[0]=(1.0/e1-1.0/e0)*s01-1.0/e0*s0;
  
  double s0 = 0, s1 = 0; // Sum of age values of initial and final network.
  Edge e0, e1; // Number of edges in initial and final network.
  
  e0 = e1 = N_EDGES;
  
  for(Edge k=1; k <= e0; k++){
    Vertex tail, head;
    FindithEdge(&tail, &head, k, nwp);
    s0 += ElapsedTime(tail,head,nwp);
  }
  
  s1=s0;
  
  FOR_EACH_TOGGLE(i){
    Vertex tail = tails[i], head = heads[i];
    if(IS_OUTEDGE(tail, head)){
      int et = ElapsedTime(tail,head,nwp); 
      s1 -= et;
      e1--;
    }else{
      e1++;
    }
  }
  
  CHANGE_STAT[0]=s1/e1-s0/e0;
}

T_CHANGESTAT_FN(t_mean_age_mon){
	CHANGE_STAT[0] = +1; // Each extant edge's age increases by 1, so their mean increases by 1.
}

S_CHANGESTAT_FN(s_mean_age_mon){
  CHANGE_STAT[0] = 0;

  if(N_EDGES>0){
    for (Edge k=1; k <= N_EDGES; k++){
      Vertex tail, head;
      FindithEdge(&tail, &head, k, nwp);
      int age = ElapsedTime(tail,head,nwp);
      CHANGE_STAT[0] += age;
    }

    CHANGE_STAT[0]/=N_EDGES;
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
  Vertex tail, head;

  ZERO_ALL_CHANGESTATS(i);
  FOR_EACH_TOGGLE(i){
    int age = ElapsedTime(tail=TAIL(i),head=HEAD(i),nwp)+1;
    int edgeflag = IS_OUTEDGE(tail, head);
    double val = INPUT_ATTRIB[(head - 1 - noffset) * nrow + (tail - 1)];
    CHANGE_STAT[0] += edgeflag ? - age*val : 0;
    TOGGLE_IF_MORE_TO_COME(i);
  }
  UNDO_PREVIOUS_TOGGLES(i);
}

T_CHANGESTAT_FN(t_edgecov_ages_mon){
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
    CHANGE_STAT[0] += val;
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

