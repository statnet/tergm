#include "changestats_duration.h"

/*****************
 void d_edges_ageinterval

 This is essentially the edges statistic, which only counts dyads with "age"
 (time steps spent in the current state) in the interval [inputparams0,inputparams1).
*****************/
D_CHANGESTAT_FN(d_edges_ageinterval){
  int i;
  
  ZERO_ALL_CHANGESTATS(i);

  FOR_EACH_TOGGLE(i){
    Vertex tail, head;
    int age = ElapsedTime(tail=TAIL(i),head=HEAD(i),nwp);
    // Only count if the age is in [from,to). ( to=0 ==> to=Inf )
    for(unsigned int j=0; j<N_CHANGE_STATS; j++){
      unsigned int from = INPUT_PARAM[j*2], to = INPUT_PARAM[j*2+1];
      if(from<=age && (to==0 || age<to)){
	unsigned int edgeflag = IS_OUTEDGE(tail, head);
	CHANGE_STAT[j] += edgeflag ? - 1 : 1;
      }
    }
    TOGGLE_IF_MORE_TO_COME(i);
  }
  UNDO_PREVIOUS_TOGGLES(i);
}

D_CHANGESTAT_FN(d_edges_ageinterval_mon){
  int i;

  ZERO_ALL_CHANGESTATS(i);

  for(Edge k=1; k <= N_EDGES; k++){
    Vertex tail, head;
    FindithEdge(&tail, &head, k, nwp);
    int age = ElapsedTime(tail,head,nwp); // Every tie "starts out" at age 1.
    for(unsigned int j=0; j<N_CHANGE_STATS; j++){
      unsigned int from = INPUT_PARAM[j*2], to = INPUT_PARAM[j*2+1];
      if(age+1 == from) CHANGE_STAT[j]++; // The tie "ages" into the interval.
      if(to!=0 && age+1 == to) CHANGE_STAT[j]--; // The tie "ages" out of the interval.
    }
  }

  FOR_EACH_TOGGLE(i){
    Vertex tail, head;
    int age = ElapsedTime(tail=TAIL(i),head=HEAD(i),nwp);
    // Only count if the age is in [from,to). ( to=0 ==> to=Inf )

    for(unsigned int j=0; j<N_CHANGE_STATS; j++){
      unsigned int from = INPUT_PARAM[j*2], to = INPUT_PARAM[j*2+1];
      if(IS_OUTEDGE(tail, head)){ // If already an edge, we are dissolving.
	if(from<=age+1 && (to==0 || age+1<to)) CHANGE_STAT[j]--; // Statistic only changes if it's in the interval.
      }else{ // If not already an edge, we are forming.
	if(from==1) CHANGE_STAT[j]++; // Statistic only changes if it starts at just-formed edges.
      }
    }
  }
}

S_CHANGESTAT_FN(s_edges_ageinterval_mon){
  int i;
  ZERO_ALL_CHANGESTATS(i);
  for (Edge k=1; k <= N_EDGES; k++){
    Vertex tail, head;
    FindithEdge(&tail, &head, k, nwp);
    int age = ElapsedTime(tail,head,nwp);
    for(unsigned int j=0; j<N_CHANGE_STATS; j++){
      unsigned int from = INPUT_PARAM[j*2], to = INPUT_PARAM[j*2+1];
      if(from<=age && (to==0 || age<to)) CHANGE_STAT[0]++;
    }
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

 Weighted sum of ages of all extant ties.

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

 Mean of ages of all extant ties.

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
    int et = ElapsedTime(tail,head,nwp);
    s0 += et;
    s1 += et + 1;
  }
  
  FOR_EACH_TOGGLE(i){
    Vertex tail = tails[i], head = heads[i];
    if(IS_OUTEDGE(tail, head)){
      s1 -= ElapsedTime(tail,head,nwp) + 1;
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

 Weigted mean of ages of all extant ties.

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
      int et = ElapsedTime(tail,head,nwp);
      s0 += et*val;
      s1 += (et + 1)*val;
      e0 += val;
    }
  }
  
  e1 = e0;
  
  FOR_EACH_TOGGLE(i){
    Vertex tail = TAIL(i), head = HEAD(i);
    double val = INPUT_ATTRIB[(head - 1 - noffset) * nrow + (tail - 1)];   
    if(val!=0){
      if(IS_OUTEDGE(tail, head)){
	s1 -= (ElapsedTime(tail,head,nwp)+1)*val;
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

/*****************
 void d_degree_mean_age

 Mean of ages of all extant ties with a particular degree.

 The degree_mean_age of a network with no actors with degree of interest is defined to be emptyval.

 *****************/

D_CHANGESTAT_FN(d_degree_mean_age_mon){
  int i;
  Vertex *id=IN_DEG, *od=OUT_DEG;
  double zeroval = INPUT_PARAM[0];
  
  for(unsigned int j = 0; j < N_CHANGE_STATS; j++){
    double s0 = 0, s1 = 0;
    Edge e0 = 0, e1 = 0;

    Vertex deg = INPUT_PARAM[j+1];
    
    for (Edge k=1; k <= N_EDGES; k++){
      Vertex tail, head;
      FindithEdge(&tail, &head, k, nwp);
      
      unsigned int w = (od[tail]+id[tail]==deg) + (od[head]+id[head]==deg);
      
      if(w){
	int et = ElapsedTime(tail,head,nwp);
	s0 += et*w;
	s1 += (et+1)*w;
	e0+=w;
	e1+=w;
      }
    }

    FOR_EACH_TOGGLE(i){
      Vertex tail, head;
      int change = IS_OUTEDGE(tail=TAIL(i), head=HEAD(i)) ? -1 : +1;
      int taildiff = (od[tail]+id[tail] + change == deg)-(od[tail]+id[tail] == deg);
      int headdiff = (od[head]+id[head] + change == deg)-(od[head]+id[head] == deg);

      Edge e;
      Vertex head1, tail1;

      // This kludge is necessary because ElapsedTime will return the wrong result for a timestamp that hasn't been updated, which it hasn't been yet.
#define TMP_SET_ET(t,h)							\
      unsigned int just_formed = FALSE;					\
      for(unsigned int l=0; l<i; l++){					\
	if(t==TAIL(l) && h==HEAD(l)){					\
	  just_formed = TRUE;						\
	  break;							\
	}								\
      }									\
      int et = just_formed ? 0 : ElapsedTime(t,h,nwp);
      
      switch(taildiff){
      case -1: // tail was previously counted, but is no longer
	STEP_THROUGH_OUTEDGES(tail, e, head1){
	  TMP_SET_ET(tail,head1);
	  s1 -= et+1;
	  e1--;
	}
	STEP_THROUGH_INEDGES(tail, e, head1){
	  TMP_SET_ET(head1,tail);
	  s1 -= et+1;
	  e1--;
	}
	break;
	// We don't need to do anything special for the focus dyad here:
	// if it's formed, then it wasn't counted in the first place;
	// if it's dissolved, then it will have been subtracted off by the previous two loops.

      case +1: // tail was previously not counted, but is now
	STEP_THROUGH_OUTEDGES(tail, e, head1){
	  TMP_SET_ET(tail,head1);
	  s1 += et+1;
	  e1++;
	}
	STEP_THROUGH_INEDGES(tail, e, head1){
	  TMP_SET_ET(head1,tail);
	  s1 += et+1;
	  e1++;
	}
	// Here, we need to handle the focus dyad:
	if(change==+1){// if it's formed, add 1 to s1
	  s1 += 1;
	  e1++;
	}else{// if it's dissolved, it had been counted in the previous two loops, and it should be subtracted
	  int et = ElapsedTime(tail,head,nwp);
	  s1 -= et+1;
	  e1--;
	}
	break;
      }

      switch(headdiff){
      case -1: // head was previously counted, but is no longer
	STEP_THROUGH_OUTEDGES(head, e, tail1){
	  TMP_SET_ET(head,tail1);
	  s1 -= et+1;
	  e1--;
	}
	STEP_THROUGH_INEDGES(head, e, tail1){
	  TMP_SET_ET(tail1,head);
	  s1 -= et+1;
	  e1--;
	}
	break;
	// We don't need to do anything special for the focus dyad here:
	// if it's formed, then it wasn't counted in the first place;
	// if it's dissolved, then it will have been subtracted off by the previous two loops.

      case +1: // head was previously not counted, but is now
	STEP_THROUGH_OUTEDGES(head, e, tail1){
	  TMP_SET_ET(head,tail1);
	  s1 += et+1;
	  e1++;
	}
	STEP_THROUGH_INEDGES(head, e, tail1){
	  TMP_SET_ET(tail1,head);
	  s1 += et+1;
	  e1++;
	}
	// Here, we need to handle the focus dyad:
	if(change==+1){// if it's formed, add 1 to s1
	  s1 += 1;
	  e1++;
	}else{// if it's dissolved, it had been counted in the previous two loops, and it should be subtracted
	  int et = ElapsedTime(tail,head,nwp);
	  s1 -= et+1;
	  e1--;
	}
	break;
      }
#undef TMP_SET_ET
      TOGGLE_IF_MORE_TO_COME(i);
    }

    UNDO_PREVIOUS_TOGGLES(i);
  
    CHANGE_STAT[j]=(e1==0?zeroval:s1/e1)-(e0==0?zeroval:s0/e0);
  }
}


S_CHANGESTAT_FN(s_degree_mean_age_mon){
  int i;
  Vertex *id=IN_DEG, *od=OUT_DEG;
  double zeroval = INPUT_PARAM[0];

  ZERO_ALL_CHANGESTATS(i);

  for(unsigned int j = 0; j < N_CHANGE_STATS; j++){
    Vertex deg = INPUT_PARAM[j+1];
    Edge e=0;

    for (Edge k=1; k <= N_EDGES; k++){
      Vertex tail, head;
      FindithEdge(&tail, &head, k, nwp);
      
      unsigned int w = (od[tail]+id[tail]==deg ? 1:0) + (od[head]+id[head]==deg ? 1:0);
      
      if(w){
	CHANGE_STAT[j] += ElapsedTime(tail,head,nwp)*w;
	e+=w;
      }
    }
    
    if(e) CHANGE_STAT[j] /= e;
    else CHANGE_STAT[j] = zeroval;
  }
}

/*****************
 void d_degree_by_attr_mean_age

 Mean of ages of all extant ties with a particular degree, by actor attribute.

 The degree_by_attr_mean_age of a network with no actors with degree of interest is defined to be emptyval.

 *****************/

D_CHANGESTAT_FN(d_degree_by_attr_mean_age_mon){
  int i;
  Vertex *id=IN_DEG, *od=OUT_DEG;
  double zeroval = INPUT_PARAM[0];
  
  for(unsigned int j = 0; j < N_CHANGE_STATS; j++){
    double s0 = 0, s1 = 0;
    Edge e0 = 0, e1 = 0;

    Vertex deg = INPUT_PARAM[2*j+1];
    int testattr = INPUT_PARAM[2*j+2];
    
    for (Edge k=1; k <= N_EDGES; k++){
      Vertex tail, head;
      FindithEdge(&tail, &head, k, nwp);
      
      Vertex taildeg = od[tail]+id[tail], headdeg = od[head]+id[head];
      int tailattr = INPUT_PARAM[2*N_CHANGE_STATS + tail]; 
      int headattr = INPUT_PARAM[2*N_CHANGE_STATS + head]; 

      unsigned int w = ((taildeg==deg && tailattr==testattr) ? 1 : 0) +
	((headdeg==deg && headattr==testattr) ? 1 : 0);
      
      if(w){
	int et = ElapsedTime(tail,head,nwp);
	s0 += et*w;
	s1 += (et+1)*w;
	e0+=w;
	e1+=w;
      }
    }

    FOR_EACH_TOGGLE(i){
      Vertex tail=TAIL(i), head=HEAD(i);
      int tailattr = INPUT_PARAM[2*N_CHANGE_STATS + tail]; 
      int headattr = INPUT_PARAM[2*N_CHANGE_STATS + head]; 

      // If neither attribute matches, this toggle has no effect on the statistic.
      if(tailattr!=testattr && headattr!=testattr){
	TOGGLE_IF_MORE_TO_COME(i); // But don't forget to make it, so that it can be reversed later.
	continue; 
      }

      int change = IS_OUTEDGE(tail, head) ? -1 : +1;
      int taildiff = (od[tail]+id[tail] + change == deg)-(od[tail]+id[tail] == deg);
      int headdiff = (od[head]+id[head] + change == deg)-(od[head]+id[head] == deg);

      Edge e;
      Vertex head1, tail1;

      // This kludge is necessary because ElapsedTime will return the wrong result for a timestamp that hasn't been updated, which it hasn't been yet.
#define TMP_SET_ET(t,h)							\
      unsigned int just_formed = FALSE;					\
      for(unsigned int l=0; l<i; l++){					\
	if(t==TAIL(l) && h==HEAD(l)){					\
	  just_formed = TRUE;						\
	  break;							\
	}								\
      }									\
      int et = just_formed ? 0 : ElapsedTime(t,h,nwp);
      
      switch(taildiff * (tailattr==testattr)){ // If tailattr!=testattr, it'll look for case 0, i.e., do nothing.
      case -1: // tail was previously counted, but is no longer
	STEP_THROUGH_OUTEDGES(tail, e, head1){
	  TMP_SET_ET(tail,head1);
	  s1 -= et+1;
	  e1--;
	}
	STEP_THROUGH_INEDGES(tail, e, head1){
	  TMP_SET_ET(head1,tail);
	  s1 -= et+1;
	  e1--;
	}
	break;
	// We don't need to do anything special for the focus dyad here:
	// if it's formed, then it wasn't counted in the first place;
	// if it's dissolved, then it will have been subtracted off by the previous two loops.

      case +1: // tail was previously not counted, but is now
	STEP_THROUGH_OUTEDGES(tail, e, head1){
	  TMP_SET_ET(tail,head1);
	  s1 += et+1;
	  e1++;
	}
	STEP_THROUGH_INEDGES(tail, e, head1){
	  TMP_SET_ET(head1,tail);
	  s1 += et+1;
	  e1++;
	}
	// Here, we need to handle the focus dyad:
	if(change==+1){// if it's formed, add 1 to s1
	  s1 += 1;
	  e1++;
	}else{// if it's dissolved, it had been counted in the previous two loops, and it should be subtracted
	  int et = ElapsedTime(tail,head,nwp);
	  s1 -= et+1;
	  e1--;
	}
	break;
      }

      switch(headdiff * (headattr==testattr)){ // If headattr!=testattr, it'll look for case 0, i.e., do nothing.
      case -1: // head was previously counted, but is no longer
	STEP_THROUGH_OUTEDGES(head, e, tail1){
	  TMP_SET_ET(head,tail1);
	  s1 -= et+1;
	  e1--;
	}
	STEP_THROUGH_INEDGES(head, e, tail1){
	  TMP_SET_ET(tail1,head);
	  s1 -= et+1;
	  e1--;
	}
	break;
	// We don't need to do anything special for the focus dyad here:
	// if it's formed, then it wasn't counted in the first place;
	// if it's dissolved, then it will have been subtracted off by the previous two loops.

      case +1: // head was previously not counted, but is now
	STEP_THROUGH_OUTEDGES(head, e, tail1){
	  TMP_SET_ET(head,tail1);
	  s1 += et+1;
	  e1++;
	}
	STEP_THROUGH_INEDGES(head, e, tail1){
	  TMP_SET_ET(tail1,head);
	  s1 += et+1;
	  e1++;
	}
	// Here, we need to handle the focus dyad:
	if(change==+1){// if it's formed, add 1 to s1
	  s1 += 1;
	  e1++;
	}else{// if it's dissolved, it had been counted in the previous two loops, and it should be subtracted
	  int et = ElapsedTime(tail,head,nwp);
	  s1 -= et+1;
	  e1--;
	}
	break;
      }
#undef TMP_SET_ET
      TOGGLE_IF_MORE_TO_COME(i);
    }

    UNDO_PREVIOUS_TOGGLES(i);
  
    CHANGE_STAT[j]=(e1==0?zeroval:s1/e1)-(e0==0?zeroval:s0/e0);
  }
}

S_CHANGESTAT_FN(s_degree_by_attr_mean_age_mon){
  int i;
  Vertex *id=IN_DEG, *od=OUT_DEG;
  double zeroval = INPUT_PARAM[0];

  ZERO_ALL_CHANGESTATS(i);

  for(unsigned int j = 0; j < N_CHANGE_STATS; j++){
    Edge e=0;

    for (Edge k=1; k <= N_EDGES; k++){
      Vertex tail, head;
      FindithEdge(&tail, &head, k, nwp);
      Vertex taildeg = od[tail]+id[tail], headdeg = od[head]+id[head];
      int tailattr = INPUT_PARAM[2*N_CHANGE_STATS + tail]; 
      int headattr = INPUT_PARAM[2*N_CHANGE_STATS + head]; 
    
      Vertex deg = INPUT_PARAM[2*j+1];
      int testattr = INPUT_PARAM[2*j+2];

      unsigned int w = ((taildeg==deg && tailattr==testattr) ? 1 : 0) +
	((headdeg==deg && headattr==testattr) ? 1 : 0);
      
      if(w){
	CHANGE_STAT[j] += ElapsedTime(tail,head,nwp)*w;
	e+=w;
      }
    }
    
    if(e) CHANGE_STAT[j] /= e;
    else CHANGE_STAT[j] = zeroval;
  }
}
