#include "MHproposals_DynMLE.h"
#include "edgelist.h"

/* Shorthand. */
#define Mtail (MHp->toggletail)
#define Mhead (MHp->togglehead)

/********************
   void MH_FormationMLE
   Propose ONLY edges not in the reference graph
***********************/
void MH_FormationMLE (MHproposal *MHp, Network *nwp) 
{  
  static Vertex nnodes;
  static Edge ndyads, nedges0;

  if(MHp->ntoggles == 0) { /* Initialize */
    MHp->ntoggles=1;
    nnodes = nwp[0].nnodes;
    ndyads = DYADCOUNT(nnodes, 0, nwp[0].directed_flag);
    nedges0 = MHp->inputs[0];
    return;
  }
  
  if(nedges0==ndyads){ /* Attempting formation on a complete graph. */
    Mtail[0]=MH_FAILED;
    Mhead[0]=MH_IMPOSSIBLE;
    return;
  }

  BD_COND_LOOP({
    /* Keep trying dyads until a one that is not an edge in the reference network is found. */
    /* Generate. */
      GetRandDyad(Mtail, Mhead, nwp);
    }, !dEdgeListSearch(Mtail[0],Mhead[0],MHp->inputs), 2);
}

/********************
   void MH_FormationMLE
   Propose ONLY edges not in the reference graph
***********************/
void MH_FormationMLETNT(MHproposal *MHp, Network *nwp) 
{  
  static Vertex nnodes;
  Vertex tail,head;
  static Edge ndyads;
  static double comp=0.5, odds;
  static Network discord;

  if(MHp->ntoggles == 0) { /* Initialize */
    MHp->ntoggles=1;
    nnodes = nwp->nnodes;
    odds = comp/(1.0-comp);
    ndyads = DYADCOUNT(nnodes, nwp->bipartite, nwp[0].directed_flag);

    Edge nedges0 = MHp->inputs[0];
    MHp->discord = (Network**) calloc(2,sizeof(Network*)); // A space for the sentinel NULL pointer.
    MHp->discord[0] = &discord;
    discord = NetworkInitializeD(MHp->inputs+1, MHp->inputs+1+nedges0, nedges0, nnodes, nwp->directed_flag, nwp->bipartite, 0, 0, NULL);
   
    for(Edge i=0; i<nwp->nedges; i++){
      FindithEdge(&tail, &head, i+1, nwp);
      ToggleEdge(tail,head, &discord);
    }
    
    return;
  }
  

  Edge nedges = nwp->nedges, ndedges = discord.nedges, nempty = ndyads-nedges;

  if(nempty==0 && ndedges==0){ /* Attempting formation on a complete graph. */
    Mtail[0]=MH_FAILED;
    Mhead[0]=MH_IMPOSSIBLE;
    return;
  }

  BD_LOOP({
    if(ndedges != 0 && unif_rand() < comp) { /* Select a discordant dyad at random */
      GetRandEdge(Mtail, Mhead, &discord);
      
      MHp->logratio += log(ndedges / ((double)nempty + 1)/comp * ((ndedges==1)? 1 : (1-comp)));
    }else{ /* select an empty dyad in nwp[0] at random */
      GetRandNonedge(Mtail, Mhead, nwp);

      if(ndedges==0){
        MHp->logratio += log(nempty*comp);
      }else{
        MHp->logratio += log(((double)nempty)/(ndedges+1) *odds);
      }
    }
    });
}

/********************
   void MH_DissolutionMLE
   Propose ONLY edges in the reference graph
***********************/
void MH_DissolutionMLETNT(MHproposal *MHp, Network *nwp) 
{  
  static Vertex nnodes;
  Vertex tail,head;
  static Edge ndyads;
  static double comp=0.5, odds;
  static Network discord;

  if(MHp->ntoggles == 0) { /* Initialize */
    MHp->ntoggles=1;
    nnodes = nwp->nnodes;
    odds = comp/(1.0-comp);
    ndyads = DYADCOUNT(nnodes, nwp->bipartite, nwp[0].directed_flag);

    Edge nedges0 = MHp->inputs[0];
    MHp->discord = (Network**) calloc(2,sizeof(Network*)); // A space for the sentinel NULL pointer.
    MHp->discord[0] = &discord;
    discord = NetworkInitializeD(MHp->inputs+1, MHp->inputs+1+nedges0, nedges0, nnodes, nwp->directed_flag, nwp->bipartite, 0, 0, NULL);
   
    for(Edge i=0; i<nwp->nedges; i++){
      FindithEdge(&tail, &head, i+1, nwp);
      ToggleEdge(tail,head, &discord);
    }
    
    return;
  }
  

  Edge nedges = nwp->nedges, ndedges = discord.nedges;

  if(nedges==0 && ndedges==0){ /* Attempting dissolution on an empty graph. */
    Mtail[0]=MH_FAILED;
    Mhead[0]=MH_IMPOSSIBLE;
    return;
  }

  BD_LOOP({
    if(ndedges != 0 && unif_rand() < comp) { /* Select a discordant dyad at random */
      GetRandEdge(Mtail, Mhead, &discord);
      
      MHp->logratio += log(ndedges / ((double)nedges + 1)/comp * ((ndedges==1)? 1 : (1-comp)));
    }else{ /* select an edge in nwp[0] at random */
      GetRandEdge(Mtail, Mhead, nwp);

      if(ndedges==0){
        MHp->logratio += log(nedges*comp);
      }else{
        MHp->logratio += log(((double)nedges)/(ndedges+1) *odds);
      }
    }
    });
}
