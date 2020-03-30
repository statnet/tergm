#include "ergm_MHproposal.h"
#include "ergm_edgelist.h"
#include "ergm_changestat.h"
#include "ergm_rlebdm.h"
#include "ergm_MHstorage.h"
#include "ergm_unsorted_edgelist.h"
#include "changestats_lasttoggle.h"
#include "tergm_model.h"

typedef struct {
  UnsrtEL *discordantEdges;
  UnsrtEL *discordantNonEdges;
  
  int in_discord;
  
} discordTNTStorage; 

MH_I_FN(Mi_discordTNT) {
  MHp->ntoggles = 1;
  
  ALLOC_STORAGE(1, discordTNTStorage, sto);
  sto->discordantEdges = UnsrtELInitialize(0, NULL, NULL, FALSE);
  sto->discordantNonEdges = UnsrtELInitialize(0, NULL, NULL, FALSE);
  // we ignore discord for this initialization (assuming a TICK will precede any proposals)
}

MH_X_FN(Mx_discordTNT) {
  GET_STORAGE(discordTNTStorage, sto);
  
  if(type == TICK) {
    // "clear" the discordant dyads
    sto->discordantEdges->nedges = 0;
    sto->discordantNonEdges->nedges = 0;
  }
}

MH_P_FN(MH_discordTNT) {
  GET_AUX_STORAGE(StoreTimeAndLasttoggle, dur_inf);
  
  GET_STORAGE(discordTNTStorage, sto);
  
  int in_discord;
  int in_network;
  
  int nedges = EDGECOUNT(nwp);
  int nddyads = kh_size(dur_inf->discord);
  
  if(nddyads == 0 || unif_rand() < 0.5) {
    // propose from network
    if(nedges == 0 || unif_rand() < 0.5) {
      // propose toggling a random dyad
      GetRandDyad(Mtail, Mhead, nwp);
      in_network = IS_OUTEDGE(Mtail[0], Mhead[0]);
      in_discord = kh_get(DyadMapInt, dur_inf->discord, THKey(dur_inf->discord, Mtail[0], Mhead[0])) != kh_none;
    } else {
      // propose toggling off an edge in network
      GetRandEdge(Mtail, Mhead, nwp);
      in_network = TRUE;
      in_discord = kh_get(DyadMapInt, dur_inf->discord, THKey(dur_inf->discord, Mtail[0], Mhead[0])) != kh_none;      
    }

    if(in_discord) {
      // need to resample to know index
      if(in_network) {
        UnsrtELGetRand(Mtail, Mhead, sto->discordantEdges);
      } else {
        UnsrtELGetRand(Mtail, Mhead, sto->discordantNonEdges);
      }
    }
  } else {
    // propose from discord
    if(unif_rand() < sto->discordantEdges->nedges/((double) nddyads)) {
      UnsrtELGetRand(Mtail, Mhead, sto->discordantEdges);
      in_network = TRUE;
    } else {
      UnsrtELGetRand(Mtail, Mhead, sto->discordantNonEdges);
      in_network = FALSE;
    }
    in_discord = TRUE;
  }
  
  // compute logratio
  
  int ndyads = DYADCOUNT(nwp);
  
  // these ignore overall factor of 1/2 which cancels out in ratio
  double forward_discord = in_discord ? 1.0/nddyads : 0;
  double backward_discord = in_discord ? 0 : 1.0/(1 + nddyads);
  
  double forward_network = in_network ? (0.5/nedges + 0.5/ndyads) : (nedges == 0 ? 1.0/ndyads : 0.5/ndyads);
  double backward_network = in_network ? (nedges == 1 ? 1.0/ndyads : 0.5/ndyads) : (0.5/(nedges + 1) + 0.5/ndyads);
  
  if(nddyads == 0) forward_network *= 2;
  if(nddyads == 1 && in_discord) backward_network *= 2;
  
  double forward = forward_discord + forward_network;
  double backward = backward_discord + backward_network;

  MHp->logratio = log(backward/forward);

  sto->in_discord = in_discord;
}

MH_U_FN(Mu_discordTNT) {
  // add or remove the dyad from the appropriate discordance edgelist
  
  GET_STORAGE(discordTNTStorage, sto);
  
  if(sto->in_discord) {
    // currently in discord; take it out
    if(edgeflag) {        
      UnsrtELDelete(tail, head, sto->discordantEdges);
    } else {
      UnsrtELDelete(tail, head, sto->discordantNonEdges);  
    }
  } else {
    // not currently in discord; add it in
    if(edgeflag) {
      UnsrtELInsert(tail, head, sto->discordantNonEdges);
    } else {
      UnsrtELInsert(tail, head, sto->discordantEdges);        
    }
  }
}

MH_F_FN(Mf_discordTNT) {
  GET_STORAGE(discordTNTStorage, sto);
  
  UnsrtELDestroy(sto->discordantNonEdges);
  UnsrtELDestroy(sto->discordantEdges);
}



/********************
   discordStratTNT
********************/

typedef struct {
  UnsrtEL **nonDiscordantELs;
  UnsrtEL **discordantELs;
  UnsrtEL **discordantNonELs;
  
  double **nodesbycode;
  
  double *pmat;
  double *nodecountsbycode;
  double *dyadcounts;
  
  double *tailtypes;
  double *headtypes;
  
  int mixingtype;
  
  int nmixtypes;
  
  int in_discord;
} discordStratTNTStorage; 


MH_I_FN(Mi_discordStratTNT) {
  // process the inputs and initialize all the edgelists in storage; set MHp->ntoggles to 1
  MHp->ntoggles = 1;
  
  int nmixtypes = MHp->inputs[0];
    
  int nattrcodes = MHp->inputs[1 + 3*nmixtypes];
  
  double *vattr = MHp->inputs + 1 + 3*nmixtypes + 1 + nattrcodes + N_NODES;
  
  ALLOC_STORAGE(1, discordStratTNTStorage, sto);
  
  sto->nonDiscordantELs = (UnsrtEL **)Calloc(nmixtypes, UnsrtEL *);
  sto->discordantELs = (UnsrtEL **)Calloc(nmixtypes, UnsrtEL *);
  sto->discordantNonELs = (UnsrtEL **)Calloc(nmixtypes, UnsrtEL *);
  
  for(int i = 0; i < nmixtypes; i++) {
    sto->nonDiscordantELs[i] = UnsrtELInitialize(0, NULL, NULL, FALSE);
    sto->discordantELs[i] = UnsrtELInitialize(0, NULL, NULL, FALSE);
    sto->discordantNonELs[i] = UnsrtELInitialize(0, NULL, NULL, FALSE);
  }
    
  double *inputindmat = MHp->inputs + 1 + 3*nmixtypes + 1 + nattrcodes + N_NODES + N_NODES + nmixtypes;  
  
  double **indmat = (double **)Calloc(nattrcodes, double *);
  indmat[0] = inputindmat;
  for(int i = 1; i < nattrcodes; i++) {
    indmat[i] = indmat[i - 1] + nattrcodes;
  }
  
  // we are treating all edges as nondiscordant, 
  // assuming a TICK will precede any proposals
  Vertex head;
  Edge e;
  for(Vertex tail = 1; tail <= N_NODES; tail++) {
    STEP_THROUGH_OUTEDGES(tail, e, head) {
      int index = indmat[(int)vattr[tail - 1]][(int)vattr[head - 1]];
      if(index >= 0) {
        UnsrtELInsert(tail, head, sto->nonDiscordantELs[index]);
      }
    }
  }
  Free(indmat);
  
  sto->nodecountsbycode = MHp->inputs + 1 + 3*nmixtypes + 1;
  
  sto->nodesbycode = (double **)Calloc(nattrcodes, double *);
  sto->nodesbycode[0] = MHp->inputs + 1 + 3*nmixtypes + 1 + nattrcodes;
  for(int i = 1; i < nattrcodes; i++) {
    sto->nodesbycode[i] = sto->nodesbycode[i - 1] + (int)sto->nodecountsbycode[i - 1];
  }
  
  sto->nmixtypes = nmixtypes;
  sto->pmat = MHp->inputs + 1 + 2*nmixtypes;
  
  sto->dyadcounts = MHp->inputs + 1 + 3*nmixtypes + 1 + nattrcodes + N_NODES + N_NODES;
  
  sto->tailtypes = MHp->inputs + 1;
  sto->headtypes = MHp->inputs + 1 + nmixtypes;
} 

MH_X_FN(Mx_discordStratTNT) {
  GET_STORAGE(discordStratTNTStorage, sto);
    
  if(type == TICK) {
    // transfer discordant edges to nondiscordant edges
    // clear discordant edges and discordant nonedges
    for(int i = 0; i < sto->nmixtypes; i++) {
      Vertex *tails = sto->discordantELs[i]->tails;
      Vertex *heads = sto->discordantELs[i]->heads;
      int nedges = sto->discordantELs[i]->nedges;
      
      // UnsrtELs start at index 1 here
      for(int j = 1; j <= nedges; j++) {
        UnsrtELInsert(tails[j], heads[j], sto->nonDiscordantELs[i]);
      }
      
      sto->discordantELs[i]->nedges = 0;
      sto->discordantNonELs[i]->nedges = 0;      
    }
  }
}

MH_P_FN(MH_discordStratTNT) {
  GET_AUX_STORAGE(StoreTimeAndLasttoggle, dur_inf);
  GET_STORAGE(discordStratTNTStorage, sto);
      
  double ur = unif_rand();
  
  // find the first mixing type i with (cumulative) probability larger than ur
  int i = 0;
  while(ur > sto->pmat[i]) {
    i++;
  }
  
  // record the mixing type of the toggle, in case it's needed in the U function later
  sto->mixingtype = i;    
  
  int tailtype = sto->tailtypes[i];
  int headtype = sto->headtypes[i];
    
  // number of edges of this mixing type
  int nedgestype = sto->nonDiscordantELs[i]->nedges + sto->discordantELs[i]->nedges;

  // number of dyads of this mixing type
  int ndyadstype = sto->dyadcounts[i];
  
  // number of discordant dyads of this mixing type
  int nddyadstype = sto->discordantNonELs[i]->nedges + sto->discordantELs[i]->nedges;
  
  // flags
  int in_discord;
  int in_network;
  
  if(nddyadstype == 0 || unif_rand() < 0.5) {
    // propose from network
    if(nedgestype == 0 || unif_rand() < 0.5) {
      // propose toggling a random dyad of the specified mixing type
      
      int tailindex = sto->nodecountsbycode[tailtype]*unif_rand();
      int headindex;
      if(tailtype == headtype) {
        // need to avoid sampling a loop
        headindex = (sto->nodecountsbycode[headtype] - 1)*unif_rand();
        if(headindex == tailindex) {
          headindex = sto->nodecountsbycode[headtype] - 1;
        }
      } else {
        // any old head will do
        headindex = sto->nodecountsbycode[headtype]*unif_rand();
      }
            
      Vertex tail = sto->nodesbycode[tailtype][tailindex];
      Vertex head = sto->nodesbycode[headtype][headindex];
      
      if(tail > head && !DIRECTED) {
        Vertex tmp = tail;
        tail = head;
        head = tmp;
      }
      
      Mtail[0] = tail;
      Mhead[0] = head;

      in_network = IS_OUTEDGE(Mtail[0], Mhead[0]);
      in_discord = kh_get(DyadMapInt, dur_inf->discord, THKey(dur_inf->discord, Mtail[0], Mhead[0])) != kh_none;
      
      // if it resides in any of the edgelists we store, we need to resample
      if(in_network) {
        if(in_discord) {
          UnsrtELGetRand(Mtail, Mhead, sto->discordantELs[i]);
        } else {
          UnsrtELGetRand(Mtail, Mhead, sto->nonDiscordantELs[i]);
        }
      } else if(in_discord) {
        UnsrtELGetRand(Mtail, Mhead, sto->discordantNonELs[i]);
      }
    } else {
      // propose toggling off an edge of the specified mixing type
      if(unif_rand() < sto->nonDiscordantELs[i]->nedges/((double) nedgestype)) {
        UnsrtELGetRand(Mtail, Mhead, sto->nonDiscordantELs[i]);
        in_discord = FALSE;
      } else {
        UnsrtELGetRand(Mtail, Mhead, sto->discordantELs[i]);
        in_discord = TRUE;
      }
      in_network = TRUE;
    }
  } else {
    // propose from discord
    if(unif_rand() < sto->discordantELs[i]->nedges/((double) nddyadstype)) {
      UnsrtELGetRand(Mtail, Mhead, sto->discordantELs[i]);
      in_network = TRUE;
    } else {
      UnsrtELGetRand(Mtail, Mhead, sto->discordantNonELs[i]);
      in_network = FALSE;
    }
    in_discord = TRUE;
  }
  
  // compute logratio
  
  // these ignore overall factor of 1/2 which cancels out in ratio
  // and overall factor of prob(mixing type), which also cancels out in ratio
  double forward_discord = in_discord ? 1.0/nddyadstype : 0;
  double backward_discord = in_discord ? 0 : 1.0/(1 + nddyadstype);
  
  double forward_network = in_network ? (0.5/nedgestype + 0.5/ndyadstype) : (nedgestype == 0 ? 1.0/ndyadstype : 0.5/ndyadstype);
  double backward_network = in_network ? (nedgestype == 1 ? 1.0/ndyadstype : 0.5/ndyadstype) : (0.5/(nedgestype + 1) + 0.5/ndyadstype);
  
  if(nddyadstype == 0) forward_network *= 2;
  if(nddyadstype == 1 && in_discord) backward_network *= 2;
  
  double forward = forward_discord + forward_network;
  double backward = backward_discord + backward_network;

  MHp->logratio = log(backward/forward);

  sto->in_discord = in_discord;
}

MH_U_FN(Mu_discordStratTNT) {
  // add or remove edge from appropriate edgelist
  GET_STORAGE(discordStratTNTStorage, sto);
  
  if(edgeflag) {
    // we are removing an existing edge
    if(sto->in_discord) {
      UnsrtELDelete(tail, head, sto->discordantELs[sto->mixingtype]);
    } else {
      UnsrtELDelete(tail, head, sto->nonDiscordantELs[sto->mixingtype]);
      UnsrtELInsert(tail, head, sto->discordantNonELs[sto->mixingtype]);
    }
  } else {
    // we are adding a new edge
    if(sto->in_discord) {
      UnsrtELInsert(tail, head, sto->nonDiscordantELs[sto->mixingtype]);
      UnsrtELDelete(tail, head, sto->discordantNonELs[sto->mixingtype]);
    } else {
      UnsrtELInsert(tail, head, sto->discordantELs[sto->mixingtype]);
    }
  }
}

MH_F_FN(Mf_discordStratTNT) {
  // Free all the things
  GET_STORAGE(discordStratTNTStorage, sto);
  
  for(int i = 0; i < sto->nmixtypes; i++) {
    UnsrtELDestroy(sto->nonDiscordantELs[i]);
    UnsrtELDestroy(sto->discordantELs[i]);
    UnsrtELDestroy(sto->discordantNonELs[i]);    
  }

  Free(sto->nonDiscordantELs);
  Free(sto->discordantELs);
  Free(sto->discordantNonELs);

  Free(sto->nodesbycode);

  // MHp->storage itself should be Freed by MHProposalDestroy
}



/********************
    MH_discordBDTNT
********************/

typedef struct {
  int *attrcounts;
  Vertex **nodesvec;
  int *nodepos;
  
  int tailtype;
  int tailmaxl;
  
  int headtype;
  int headmaxl;
  
  int currentdyads;
  int proposeddyads;
  
  int bound;
  int nmixtypes;
  double *vattr;
  double *tailtypes;
  double *headtypes;
  
  Network *BDTDNE;
  Network *nonBDTDNE;
  
  UnsrtEL *discordantEdges;
  
  int in_discord;
  
  UnsrtEL *transferEL; // for moving edges between BDTDNE and nonBDTDNE
} discordBDTNTStorage;

MH_I_FN(Mi_discordBDTNT) {
  // process the inputs and initialize all the node lists in storage; set MHp->ntoggles to 1
  MHp->ntoggles = 1;
  
  int bound = MHp->inputs[0];
  int nlevels = MHp->inputs[1];
  
  double *nodecountsbycode = MHp->inputs + 2;
  
  int nmixtypes = MHp->inputs[2 + nlevels];
  
  double *tailtypes = MHp->inputs + 3 + nlevels;
  double *headtypes = tailtypes + nmixtypes;
  
  double *vattr = headtypes + nmixtypes;
    
  Vertex **nodesvec = (Vertex **)Calloc(nlevels, Vertex *);
  
  int *attrcounts = (int *)Calloc(nlevels, int);
    
  for(int i = 0; i < nlevels; i++) {
    // make room for maximum number of nodes of each type
    nodesvec[i] = (Vertex *)Calloc((int)nodecountsbycode[i], Vertex);
  }

  // nodepos[i - 1] is node i's position within nodesvec[(int)vattr[i - 1]] if node i is 
  // present within that list, and is undefined otherwise
  int *nodepos = Calloc(N_NODES, int);

  for(Vertex vertex = 1; vertex <= N_NODES; vertex++) {
    if(IN_DEG[vertex] + OUT_DEG[vertex] < bound) {
      // add vertex to the submaximal list corresponding to its attribute type
      nodesvec[(int)vattr[vertex - 1]][attrcounts[(int)vattr[vertex - 1]]] = vertex;
      nodepos[vertex - 1] = attrcounts[(int)vattr[vertex - 1]];
      attrcounts[(int)vattr[vertex - 1]]++;
    }
  }
  
  // count number of "BD-toggleable" dyads in current network
  int currentdyads = 0;    
  for(int i = 0; i < nmixtypes; i++) {
    if(tailtypes[i] == headtypes[i]) {
      currentdyads += attrcounts[(int)tailtypes[i]]*(attrcounts[(int)headtypes[i]] - 1)/2;
    } else {
      currentdyads += attrcounts[(int)tailtypes[i]]*attrcounts[(int)headtypes[i]];
    }
  }

  // if we cannot toggle any edges or dyads, error
  if(EDGECOUNT(nwp) == 0 && currentdyads == 0) {
    MHp->ntoggles = MH_FAILED;
    return;
  }  
  
  ALLOC_STORAGE(1, discordBDTNTStorage, sto);
  
  sto->attrcounts = attrcounts;
  sto->nodesvec = nodesvec;
  sto->currentdyads = currentdyads;
  sto->bound = bound;
  sto->nmixtypes = nmixtypes;
  sto->vattr = vattr;
  sto->tailtypes = tailtypes;
  sto->headtypes = headtypes;
  
  sto->nodepos = nodepos;
  
  sto->BDTDNE = NetworkInitialize(NULL, NULL, 0, N_NODES, DIRECTED, BIPARTITE, FALSE, 0, NULL);
  sto->nonBDTDNE = NetworkInitialize(NULL, NULL, 0, N_NODES, DIRECTED, BIPARTITE, FALSE, 0, NULL);
  
  sto->discordantEdges = UnsrtELInitialize(0, NULL, NULL, FALSE);
  sto->transferEL = UnsrtELInitialize(0, NULL, NULL, FALSE);
}

MH_X_FN(Mx_discordBDTNT) {
  if(type == TICK) {
    // clear all the discordance information
    
    GET_STORAGE(discordBDTNTStorage, sto);
    
    Vertex tail;
    Vertex head;
    
    for(Edge i = EDGECOUNT(sto->BDTDNE); i > 0; i--) {
      FindithEdge(&tail, &head, i, sto->BDTDNE);
      ToggleKnownEdge(tail, head, sto->BDTDNE, TRUE);
    }
    
    for(Edge i = EDGECOUNT(sto->nonBDTDNE); i > 0; i--) {
      FindithEdge(&tail, &head, i, sto->nonBDTDNE);
      ToggleKnownEdge(tail, head, sto->nonBDTDNE, TRUE);
    }
    
    sto->discordantEdges->nedges = 0;
  }
}

MH_P_FN(MH_discordBDTNT) {    
  GET_STORAGE(discordBDTNTStorage, sto);
  GET_AUX_STORAGE(StoreTimeAndLasttoggle, dur_inf);

  int nedges = EDGECOUNT(nwp);
  
  int in_network;
  int in_discord;

  // number of *BD-toggleable* discordant dyads
  int nddyads = sto->discordantEdges->nedges + EDGECOUNT(sto->BDTDNE);

  if(nddyads == 0 || unif_rand() < 0.5) {
    // propose from network

    // if currentdyads == 0, we *must* propose toggling off an existing edge;
    // the case nedges == 0 && currentdyads == 0 was excluded during initialization,
    // and we cannot end up in that case if we don't start in that case
    // (assuming the initial network is valid)  
    if((unif_rand() < 0.5 && nedges > 0) || (sto->currentdyads == 0)) {
      // select an existing edge at random, and propose toggling it off
      GetRandEdge(Mtail, Mhead, nwp);
            
      in_network = TRUE;
    } else {
      // select a BD-toggleable dyad and propose toggling it
      // doubling here allows more efficient calculation in 
      // the case tailtypes[i] == headtypes[i]
      
      int dyadindex = 2*sto->currentdyads*unif_rand();
              
      Vertex tail;
      Vertex head;
      
      // this rather ugly block of code is just finding the dyad that corresponds
      // to the dyadindex we drew above, and then setting the info for
      // tail and head appropriately
      for(int i = 0; i < sto->nmixtypes; i++) {
        int dyadstype;
        if(sto->tailtypes[i] == sto->headtypes[i]) {
          dyadstype = sto->attrcounts[(int)sto->tailtypes[i]]*(sto->attrcounts[(int)sto->headtypes[i]] - 1)/2;
        } else {
          dyadstype = sto->attrcounts[(int)sto->tailtypes[i]]*sto->attrcounts[(int)sto->headtypes[i]];
        }
        
        if(dyadindex < 2*dyadstype) {
          int tailindex;
          int headindex;
          
          if(sto->tailtypes[i] == sto->headtypes[i]) {
            tailindex = dyadindex / sto->attrcounts[(int)sto->headtypes[i]];
            headindex = dyadindex % (sto->attrcounts[(int)sto->headtypes[i]] - 1);
            if(tailindex == headindex) {
              headindex = sto->attrcounts[(int)sto->headtypes[i]] - 1;
            }
                      
            tail = sto->nodesvec[(int)sto->tailtypes[i]][tailindex];
            head = sto->nodesvec[(int)sto->headtypes[i]][headindex];
            
          } else {
            dyadindex /= 2;
            tailindex = dyadindex / sto->attrcounts[(int)sto->headtypes[i]];
            headindex = dyadindex % sto->attrcounts[(int)sto->headtypes[i]];
            
            tail = sto->nodesvec[(int)sto->tailtypes[i]][tailindex];
            head = sto->nodesvec[(int)sto->headtypes[i]][headindex];
          }
          
          if(tail > head) {
            Mtail[0] = head;
            Mhead[0] = tail;
          } else {
            Mtail[0] = tail;
            Mhead[0] = head;
          }
          
          break;
        } else {
          dyadindex -= 2*dyadstype;
        }
      }
          
      in_network = IS_OUTEDGE(Mtail[0],Mhead[0]);
    }
    
    in_discord = kh_get(DyadMapInt, dur_inf->discord, THKey(dur_inf->discord, Mtail[0], Mhead[0])) != kh_none;

    if(in_network && in_discord) {
      // need to resample
      UnsrtELGetRand(Mtail, Mhead, sto->discordantEdges);
    }
  } else {
    // propose from (the BD-toggleable subset of) discord
    if(unif_rand() < ((double) sto->discordantEdges->nedges)/nddyads) {
      UnsrtELGetRand(Mtail, Mhead, sto->discordantEdges);
      in_network = TRUE;
    } else {
      GetRandEdge(Mtail, Mhead, sto->BDTDNE);
      in_network = FALSE;
    }
    
    in_discord = TRUE;
  }
  
  sto->in_discord = in_discord;
  
  sto->tailtype = sto->vattr[Mtail[0] - 1];
  sto->headtype = sto->vattr[Mhead[0] - 1];      
  
  if(in_network) {
    sto->tailmaxl = IN_DEG[Mtail[0]] + OUT_DEG[Mtail[0]] == sto->bound;
    sto->headmaxl = IN_DEG[Mhead[0]] + OUT_DEG[Mhead[0]] == sto->bound;   
  } else {
    sto->tailmaxl = IN_DEG[Mtail[0]] + OUT_DEG[Mtail[0]] == sto->bound - 1;
    sto->headmaxl = IN_DEG[Mhead[0]] + OUT_DEG[Mhead[0]] == sto->bound - 1;
  }
        
  // the count of dyads that can be toggled in the "GetRandBDDyad" branch,
  // in the proposed network
  sto->proposeddyads = 0;
  
  for(int i = 0; i < sto->nmixtypes; i++) {
    int proposedtailadjustment = (sto->vattr[Mtail[0] - 1] == sto->tailtypes[i] && sto->tailmaxl) + (sto->vattr[Mhead[0] - 1] == sto->tailtypes[i] && sto->headmaxl);
    int proposedheadadjustment = (sto->vattr[Mtail[0] - 1] == sto->headtypes[i] && sto->tailmaxl) + (sto->vattr[Mhead[0] - 1] == sto->headtypes[i] && sto->headmaxl);
    
    if(!in_network) {
      proposedtailadjustment = -proposedtailadjustment;
      proposedheadadjustment = -proposedheadadjustment;
    }
    
    if(sto->tailtypes[i] == sto->headtypes[i]) {
      sto->proposeddyads += (sto->attrcounts[(int)sto->tailtypes[i]] + proposedtailadjustment)*(sto->attrcounts[(int)sto->headtypes[i]] + proposedheadadjustment - 1)/2;
    } else {
      sto->proposeddyads += (sto->attrcounts[(int)sto->tailtypes[i]] + proposedtailadjustment)*(sto->attrcounts[(int)sto->headtypes[i]] + proposedheadadjustment);
    }
  }

  // how nddyads can change:
  // the dyad we toggle can add/remove one discordant dyad
  // discordant nonedges can change toggleability status based on the toggle we make,
  // but only if they have Mtail or Mhead as an endpoint
  // so check BDTDNE and/or nonBDTDNE neighbors of Mtail and Mhead, other than Mtail -> Mhead itself
  // which should be taken into account at the outset
  int propnddyads = nddyads;
  
  // direct effect of the toggle we are proposing
  if(in_discord) {
    propnddyads--;
  } else {
    propnddyads++;
  }
  
  // indirect effect of causing other discordant nonedges to change toggleability status
  if(!in_network) {
    // may reduce propnddyads; subtract in_discord to avoid repeatedly counting the Mtail[0] -> Mhead[0] edge
    if(sto->tailmaxl) {
      propnddyads -= sto->BDTDNE->indegree[Mtail[0]] + sto->BDTDNE->outdegree[Mtail[0]] - in_discord;
    }
    
    if(sto->headmaxl) {
      propnddyads -= sto->BDTDNE->indegree[Mhead[0]] + sto->BDTDNE->outdegree[Mhead[0]] - in_discord;
    }
  } else {
    // may increase propnddyads, but only if other endpoint is also submaximal
    if(sto->tailmaxl) {
      if(sto->nonBDTDNE->indegree[Mtail[0]] > 0) {
        Edge e = Mtail[0];
        do { // NB the neighbor cannot be head
          if(IN_DEG[sto->nonBDTDNE->inedges[e].value] + OUT_DEG[sto->nonBDTDNE->inedges[e].value] < sto->bound) {
            propnddyads++;
          }
        } while((e = EdgetreePreSuccessor(sto->nonBDTDNE->inedges, e)) != 0);       
      }
      
      if(sto->nonBDTDNE->outdegree[Mtail[0]] > 0) {
        Edge e = Mtail[0];
        do { // NB the neighbor cannot be head
          if(IN_DEG[sto->nonBDTDNE->outedges[e].value] + OUT_DEG[sto->nonBDTDNE->outedges[e].value] < sto->bound) {
            propnddyads++;
          }
        } while((e = EdgetreePreSuccessor(sto->nonBDTDNE->outedges, e)) != 0);               
      }
    }
    
    if(sto->headmaxl) {
      if(sto->nonBDTDNE->indegree[Mhead[0]] > 0) {
        Edge e = Mhead[0];
        do { // NB the neighbor cannot be head
          if(IN_DEG[sto->nonBDTDNE->inedges[e].value] + OUT_DEG[sto->nonBDTDNE->inedges[e].value] < sto->bound) {
            propnddyads++;
          }
        } while((e = EdgetreePreSuccessor(sto->nonBDTDNE->inedges, e)) != 0);       
      }
      
      if(sto->nonBDTDNE->outdegree[Mhead[0]] > 0) {
        Edge e = Mhead[0];
        do { // NB the neighbor cannot be head
          if(IN_DEG[sto->nonBDTDNE->outedges[e].value] + OUT_DEG[sto->nonBDTDNE->outedges[e].value] < sto->bound) {
            propnddyads++;
          }
        } while((e = EdgetreePreSuccessor(sto->nonBDTDNE->outedges, e)) != 0);               
      }
    }
  }
  
  double forward_network = in_network ? ((sto->currentdyads == 0 ? 1.0 : 0.5)/nedges + (sto->tailmaxl || sto->headmaxl ? 0 : 0.5/sto->currentdyads)) : ((nedges == 0 ? 1.0 : 0.5)/sto->currentdyads);
  
  double forward_discord = in_discord ? 1.0/nddyads : 0;
  
  double backward_network = in_network ? ((nedges == 1 ? 1.0 : 0.5)/sto->proposeddyads) : ((sto->proposeddyads == 0 ? 1.0 : 0.5)/(nedges + 1) + (sto->tailmaxl || sto->headmaxl ? 0 : 0.5/sto->proposeddyads));
  
  double backward_discord = in_discord ? 0 : 1.0/propnddyads;
  
  if(nddyads == 0) forward_network *= 2;
  if(propnddyads == 0) backward_network *= 2;
  
  double forward = forward_discord + forward_network;
  double backward = backward_discord + backward_network;

  MHp->logratio = log(backward/forward);
}

// this U_FN is called *before* the toggle is made in the network
MH_U_FN(Mu_discordBDTNT) {  
  GET_STORAGE(discordBDTNTStorage, sto);
  
  // add or remove the dyad being toggled from the relevant edge set
  if(edgeflag) {
    if(sto->in_discord) {
      // discordant edge
      UnsrtELDelete(tail, head, sto->discordantEdges);
    } else {
      // nondiscordant edge; will become BDTDNE
      ToggleKnownEdge(tail, head, sto->BDTDNE, FALSE);
    }
  } else {
    if(sto->in_discord) {
      // discordant non-edge; evidently BD toggleable
      ToggleKnownEdge(tail, head, sto->BDTDNE, TRUE); // it's an edge in the (BDT) discordant non-edge network :)        
    } else {
      // nondiscordant nonedge; will become discordantEdge
      UnsrtELInsert(tail, head, sto->discordantEdges);        
    }
  }
  
  // update current dyad count
  sto->currentdyads = sto->proposeddyads;  
  
  // if (sub)maximality/BD toggleability has changed for any nodes/dyads, update storage accordingly
  if(edgeflag) {
    if(sto->tailmaxl) {
      // tail will be newly submaxl after toggle, so add it to the appropriate node list
      sto->nodesvec[sto->tailtype][sto->attrcounts[sto->tailtype]] = tail;
      sto->nodepos[tail - 1] = sto->attrcounts[sto->tailtype];
      sto->attrcounts[sto->tailtype]++;

      sto->transferEL->nedges = 0;
      
      // iterate over all nonBDTDNEs with tail as an endpoint; if other endpoint is also
      // submaxl, move this edge from nonBDTDNE to BDTDNE
      if(sto->nonBDTDNE->indegree[tail] > 0) {
        Edge e = tail;
        do { // NB the neighbor cannot be head
          if(IN_DEG[sto->nonBDTDNE->inedges[e].value] + OUT_DEG[sto->nonBDTDNE->inedges[e].value] < sto->bound) {
            UnsrtELInsert(sto->nonBDTDNE->inedges[e].value, tail, sto->transferEL);
          }
        } while((e = EdgetreePreSuccessor(sto->nonBDTDNE->inedges, e)) != 0);
      }
      
      if(sto->nonBDTDNE->outdegree[tail] > 0) {
        Edge e = tail;
        do {
          if(IN_DEG[sto->nonBDTDNE->outedges[e].value] + OUT_DEG[sto->nonBDTDNE->outedges[e].value] < sto->bound) {
            UnsrtELInsert(tail, sto->nonBDTDNE->outedges[e].value, sto->transferEL);
          }
        } while((e = EdgetreePreSuccessor(sto->nonBDTDNE->outedges, e)) != 0);
      }                

      for(int i = 1; i <= sto->transferEL->nedges; i++) {
        ToggleKnownEdge(sto->transferEL->tails[i], sto->transferEL->heads[i], sto->nonBDTDNE, TRUE);            
        ToggleKnownEdge(sto->transferEL->tails[i], sto->transferEL->heads[i], sto->BDTDNE, FALSE);
      }        
    }
    
    if(sto->headmaxl) {
      // head will be newly submaxl after toggle, so add it to the appropriate node list    
      sto->nodesvec[sto->headtype][sto->attrcounts[sto->headtype]] = head;
      sto->nodepos[head - 1] = sto->attrcounts[sto->headtype];
      sto->attrcounts[sto->headtype]++;

      sto->transferEL->nedges = 0;

      // iterate over all nonBDTDNEs with head as an endpoint; if other endpoint is also
      // submaxl, move this edge from nonBDTDNE to BDTDNE
      if(sto->nonBDTDNE->indegree[head] > 0) {
        Edge e = head;
        do { // NB the neighbor cannot be tail
          if(IN_DEG[sto->nonBDTDNE->inedges[e].value] + OUT_DEG[sto->nonBDTDNE->inedges[e].value] < sto->bound) {
            UnsrtELInsert(sto->nonBDTDNE->inedges[e].value, head, sto->transferEL);
          }
        } while((e = EdgetreePreSuccessor(sto->nonBDTDNE->inedges, e)) != 0);
      }
      
      if(sto->nonBDTDNE->outdegree[head] > 0) {
        Edge e = head;
        do {
          if(IN_DEG[sto->nonBDTDNE->outedges[e].value] + OUT_DEG[sto->nonBDTDNE->outedges[e].value] < sto->bound) {
            UnsrtELInsert(head, sto->nonBDTDNE->outedges[e].value, sto->transferEL);
          }
        } while((e = EdgetreePreSuccessor(sto->nonBDTDNE->outedges, e)) != 0);
      }                

      for(int i = 1; i <= sto->transferEL->nedges; i++) {
        ToggleKnownEdge(sto->transferEL->tails[i], sto->transferEL->heads[i], sto->nonBDTDNE, TRUE);            
        ToggleKnownEdge(sto->transferEL->tails[i], sto->transferEL->heads[i], sto->BDTDNE, FALSE);
      }
    }
  } else {
    if(sto->tailmaxl) {
      // tail will be newly maxl after toggle, so remove it from the appropriate node list, updating nodepos
      // for whatever node is currently at the end of sto->nodesvec[sto->tailtype], since that will be
      // moved into tail's position
      sto->nodepos[sto->nodesvec[sto->tailtype][sto->attrcounts[sto->tailtype] - 1] - 1] = sto->nodepos[tail - 1];
      sto->nodesvec[sto->tailtype][sto->nodepos[tail - 1]] = sto->nodesvec[sto->tailtype][sto->attrcounts[sto->tailtype] - 1];
      sto->attrcounts[sto->tailtype]--;

      sto->transferEL->nedges = 0;
      
      // transfer all BDTDNEs with tail as an endpoint to nonBDTDNE
      if(sto->BDTDNE->indegree[tail] > 0) {
        Edge e = tail;
        do {
          UnsrtELInsert(sto->BDTDNE->inedges[e].value, tail, sto->transferEL);
        } while((e = EdgetreePreSuccessor(sto->BDTDNE->inedges, e)) != 0);
      }
      
      if(sto->BDTDNE->outdegree[tail] > 0) {
        Edge e = tail;
        do {
          UnsrtELInsert(sto->BDTDNE->outedges[e].value, tail, sto->transferEL);
        } while((e = EdgetreePreSuccessor(sto->BDTDNE->outedges, e)) != 0);
      }                

      for(int i = 1; i <= sto->transferEL->nedges; i++) {
        ToggleKnownEdge(sto->transferEL->tails[i], sto->transferEL->heads[i], sto->BDTDNE, TRUE);
        ToggleKnownEdge(sto->transferEL->tails[i], sto->transferEL->heads[i], sto->nonBDTDNE, FALSE);            
      }
    }
    
    if(sto->headmaxl) {
      // head will be newly maxl after toggle, so remove it from the appropriate node list, updating nodepos
      // for whatever node is currently at the end of sto->nodesvec[sto->headtype], since that will be
      // moved into head's position
      sto->nodepos[sto->nodesvec[sto->headtype][sto->attrcounts[sto->headtype] - 1] - 1] = sto->nodepos[head - 1];
      sto->nodesvec[sto->headtype][sto->nodepos[head - 1]] = sto->nodesvec[sto->headtype][sto->attrcounts[sto->headtype] - 1];
      sto->attrcounts[sto->headtype]--;

      sto->transferEL->nedges = 0;
      
      // transfer all BDTDNEs with head as an endpoint to nonBDTDNE
      if(sto->BDTDNE->indegree[head] > 0) {
        Edge e = head;
        do {
          UnsrtELInsert(sto->BDTDNE->inedges[e].value, head, sto->transferEL);
        } while((e = EdgetreePreSuccessor(sto->BDTDNE->inedges, e)) != 0);
      }
      
      if(sto->BDTDNE->outdegree[head] > 0) {
        Edge e = head;
        do {
          UnsrtELInsert(sto->BDTDNE->outedges[e].value, head, sto->transferEL);
        } while((e = EdgetreePreSuccessor(sto->BDTDNE->outedges, e)) != 0);
      }                

      for(int i = 1; i <= sto->transferEL->nedges; i++) {
        ToggleKnownEdge(sto->transferEL->tails[i], sto->transferEL->heads[i], sto->BDTDNE, TRUE);
        ToggleKnownEdge(sto->transferEL->tails[i], sto->transferEL->heads[i], sto->nonBDTDNE, FALSE);            
      }
    }      
  }
}

MH_F_FN(Mf_discordBDTNT) {
  // Free all the things
  GET_STORAGE(discordBDTNTStorage, sto);
  
  int nlevels = MHp->inputs[1];

  for(int i = 0; i < nlevels; i++) {
    Free(sto->nodesvec[i]);
  }

  Free(sto->nodesvec);
  Free(sto->attrcounts);
  Free(sto->nodepos);

  NetworkDestroy(sto->BDTDNE);
  NetworkDestroy(sto->nonBDTDNE);
  UnsrtELDestroy(sto->discordantEdges);
  
  // MHp->storage itself should be Freed by MHProposalDestroy
}



/********************
    MH_discordBDStratTNT
********************/

typedef struct {  
  Vertex ***nodesvec;
  int **attrcounts;
  int *nodepos;
  
  int strattailtype;
  int bdtailtype;
  int tailmaxl;
  
  int stratheadtype;
  int bdheadtype;  
  int headmaxl;
  
  int stratmixingtype;
  
  int *currentdyads;
  int *proposeddyads;
  
  double *originalprobvec;
  double *currentprobvec;
  double *proposedprobvec;
  
  int bound;
  int nmixtypes;
  
  double *strat_vattr;
  double *bd_vattr;
  
  double *BDtypesbyStrattype;
  double *BDtailsbyStrattype;
  double *BDheadsbyStrattype;
  
  double *strattailtypes;
  double *stratheadtypes;
  
  Network **BDTDNE;
  Network **nonBDTDNE;
  
  UnsrtEL **discordantEdges;
  UnsrtEL **nonDiscordantEdges;
  
  int in_discord;
  
  UnsrtEL *transferEL; // for moving edges between BDTDNE and nonBDTDNE  
} discordBDStratTNTStorage;

MH_I_FN(Mi_discordBDStratTNT) {
  // process the inputs and initialize all the edgelists in storage; set MHp->ntoggles to 1
  MHp->ntoggles = 1;
  
  int nmixtypes = MHp->inputs[0];
  
  double *probvec = MHp->inputs + 1 + 2*nmixtypes;
  
  double *strattailattrs = MHp->inputs + 1;
  double *stratheadattrs = MHp->inputs + 1 + nmixtypes;
  
  int nattrcodes = MHp->inputs[1 + 3*nmixtypes];
  
  double *strat_vattr = MHp->inputs + 1 + 3*nmixtypes + 1;
  
  UnsrtEL **nonDiscordantEdges = (UnsrtEL **)Calloc(nmixtypes, UnsrtEL *);
      
  for(int i = 0; i < nmixtypes; i++) {
    nonDiscordantEdges[i] = UnsrtELInitialize(0, NULL, NULL, FALSE);
  }
  
  double *inputindmat = MHp->inputs + 1 + 3*nmixtypes + 1 + N_NODES;  
  
  double **indmat = (double **)Calloc(nattrcodes, double *);
  indmat[0] = inputindmat;
  for(int i = 1; i < nattrcodes; i++) {
    indmat[i] = indmat[i - 1] + nattrcodes;
  }
  
  Vertex head;
  Edge e;
  for(Vertex tail = 1; tail <= N_NODES; tail++) {
    STEP_THROUGH_OUTEDGES(tail, e, head) {
      int index = indmat[(int)strat_vattr[tail - 1]][(int)strat_vattr[head - 1]];
      if(index >= 0) {
        UnsrtELInsert(tail, head, nonDiscordantEdges[index]);
      }
    }
  }
  Free(indmat);
  
  // above handles initialization of edgelists according to the "strat" part of BDStratTNT
  
  int npairings = MHp->inputs[1 + 3*nmixtypes + 1 + N_NODES + nattrcodes*nattrcodes];
  
  int bound = MHp->inputs[1 + 3*nmixtypes + 1 + N_NODES + nattrcodes*nattrcodes + 1 + npairings];
  
  int bdlevels = MHp->inputs[1 + 3*nmixtypes + 1 + N_NODES + nattrcodes*nattrcodes + 1 + npairings + 1];
  
  int bdmixtypes = MHp->inputs[1 + 3*nmixtypes + 1 + N_NODES + nattrcodes*nattrcodes + 1 + npairings + 1 + 1];
  
  double *bd_vattr = MHp->inputs + 1 + 3*nmixtypes + 1 + N_NODES + nattrcodes*nattrcodes + 1 + npairings + 1 + 1 + 1 + 2*bdmixtypes;
  
  double *nodecountsbypairedcode = MHp->inputs + 1 + 3*nmixtypes + 1 + N_NODES + nattrcodes*nattrcodes + 1;
  
  
  double **nodecountsmat = (double **)Calloc(nattrcodes, double *);
  nodecountsmat[0] = nodecountsbypairedcode;
  for(int i = 1; i < nattrcodes; i++) {
    nodecountsmat[i] = nodecountsmat[i - 1] + bdlevels;
  }

  Vertex ***nodesvec = (Vertex ***)Calloc(nattrcodes, Vertex **);
  for(int i = 0; i < nattrcodes; i++) {
    nodesvec[i] = (Vertex **)Calloc(bdlevels, Vertex *);
    for(int j = 0; j < bdlevels; j++) {
      nodesvec[i][j] = (Vertex *)Calloc(nodecountsmat[i][j], Vertex);
    }
  }
  Free(nodecountsmat);


  int **attrcounts = (int **)Calloc(nattrcodes, int *);
  for(int i = 0; i < nattrcodes; i++) {
    attrcounts[i] = (int *)Calloc(bdlevels, int);
  }
  
  int *nodepos = Calloc(N_NODES, int);
  for(Vertex vertex = 1; vertex <= N_NODES; vertex++) {
    if(IN_DEG[vertex] + OUT_DEG[vertex] < bound) {
      // add vertex to the submaximal list corresponding to its attribute type
      nodesvec[(int)strat_vattr[vertex - 1]][(int)bd_vattr[vertex - 1]][attrcounts[(int)strat_vattr[vertex - 1]][(int)bd_vattr[vertex - 1]]] = vertex;
      nodepos[vertex - 1] = attrcounts[(int)strat_vattr[vertex - 1]][(int)bd_vattr[vertex - 1]];
      attrcounts[(int)strat_vattr[vertex - 1]][(int)bd_vattr[vertex - 1]]++;
    }
  }

  double *BDtypesbyStrattype = MHp->inputs + 1 + 3*nmixtypes + 1 + N_NODES + nattrcodes*nattrcodes + 1 + npairings + 1 + 1 + 1 + 2*bdmixtypes + N_NODES;
  
  int sumtypes = MHp->inputs[1 + 3*nmixtypes + 1 + N_NODES + nattrcodes*nattrcodes + 1 + npairings + 1 + 1 + 1 + 2*bdmixtypes + N_NODES + nmixtypes];
  
  double *BDtailsbyStrattype = MHp->inputs + 1 + 3*nmixtypes + 1 + N_NODES + nattrcodes*nattrcodes + 1 + npairings + 1 + 1 + 1 + 2*bdmixtypes + N_NODES + nmixtypes + 1;
  
  double *BDheadsbyStrattype = MHp->inputs + 1 + 3*nmixtypes + 1 + N_NODES + nattrcodes*nattrcodes + 1 + npairings + 1 + 1 + 1 + 2*bdmixtypes + N_NODES + nmixtypes + 1 + sumtypes;
  
  int sumcurrentdyads = 0;
  int *currentdyads = Calloc(nmixtypes, int);
  
  double *currentprobvec = Calloc(nmixtypes, double);
  double sumprobs = 0;
  
  for(int i = 0; i < nmixtypes; i++) {
    for(int j = 0; j < BDtypesbyStrattype[i]; j++) {
      int tailcounts = attrcounts[(int)strattailattrs[i]][(int)BDtailsbyStrattype[j]];
      int headcounts = attrcounts[(int)stratheadattrs[i]][(int)BDheadsbyStrattype[j]];
      
      if(strattailattrs[i] != stratheadattrs[i] || BDtailsbyStrattype[j] != BDheadsbyStrattype[j]) {
        currentdyads[i] += tailcounts*headcounts;
      } else {
        currentdyads[i] += tailcounts*(headcounts - 1)/2;
      }
    }
    
    BDtailsbyStrattype += (int)BDtypesbyStrattype[i];
    BDheadsbyStrattype += (int)BDtypesbyStrattype[i];
    sumcurrentdyads += currentdyads[i];
    
    if(currentdyads[i] > 0 || nonDiscordantEdges[i]->nedges > 0) {
      currentprobvec[i] = probvec[i];
      sumprobs += probvec[i];
    } // else it's already 0
  }
  
  // if we cannot toggle any edges or dyads, error
  if(EDGECOUNT(nwp) == 0 && sumcurrentdyads == 0) {
    MHp->ntoggles = MH_FAILED;
    return;
  }

  for(int i = 0; i < nmixtypes; i++) {
    currentprobvec[i] /= sumprobs;
  }

  ALLOC_STORAGE(1, discordBDStratTNTStorage, sto);
  
  sto->nonDiscordantEdges = nonDiscordantEdges;
  sto->nodesvec = nodesvec;
  sto->attrcounts = attrcounts;
  
  sto->currentdyads = currentdyads;
  sto->proposeddyads = Calloc(nmixtypes, int);
  
  sto->originalprobvec = probvec;
  sto->currentprobvec = currentprobvec;
  sto->proposedprobvec = Calloc(nmixtypes, double);
  
  sto->bound = bound;
  sto->nmixtypes = nmixtypes;
  
  sto->strat_vattr = strat_vattr;
  sto->bd_vattr = bd_vattr;
  
  sto->BDtypesbyStrattype = BDtypesbyStrattype;
  sto->BDtailsbyStrattype = MHp->inputs + 1 + 3*nmixtypes + 1 + N_NODES + nattrcodes*nattrcodes + 1 + npairings + 1 + 1 + 1 + 2*bdmixtypes + N_NODES + nmixtypes + 1;
  sto->BDheadsbyStrattype = MHp->inputs + 1 + 3*nmixtypes + 1 + N_NODES + nattrcodes*nattrcodes + 1 + npairings + 1 + 1 + 1 + 2*bdmixtypes + N_NODES + nmixtypes + 1 + sumtypes;
  
  sto->strattailtypes = MHp->inputs + 1;
  sto->stratheadtypes = MHp->inputs + 1 + nmixtypes;
  
  sto->nodepos = nodepos;
  
  sto->BDTDNE = Calloc(nmixtypes, Network *);
  sto->nonBDTDNE = Calloc(nmixtypes, Network *);
  sto->discordantEdges = Calloc(nmixtypes, UnsrtEL *);
  
  for(int i = 0; i < nmixtypes; i++) {
    sto->BDTDNE[i] = NetworkInitialize(NULL, NULL, 0, N_NODES, DIRECTED, BIPARTITE, FALSE, 0, NULL);
    sto->nonBDTDNE[i] = NetworkInitialize(NULL, NULL, 0, N_NODES, DIRECTED, BIPARTITE, FALSE, 0, NULL);
    
    sto->discordantEdges[i] = UnsrtELInitialize(0, NULL, NULL, FALSE);
  }
  
  sto->transferEL = UnsrtELInitialize(0, NULL, NULL, FALSE);
}

MH_X_FN(Mx_discordBDStratTNT) {
  if(type == TICK) {
    GET_STORAGE(discordBDStratTNTStorage, sto);

    Vertex tail;
    Vertex head;

    for(int i = 0; i < sto->nmixtypes; i++) {
      // clear all the discordance information
            
      for(Edge j = EDGECOUNT(sto->BDTDNE[i]); j > 0; j--) {
        FindithEdge(&tail, &head, j, sto->BDTDNE[i]);
        ToggleKnownEdge(tail, head, sto->BDTDNE[i], TRUE);
      }
      
      for(Edge j = EDGECOUNT(sto->nonBDTDNE[i]); j > 0; j--) {
        FindithEdge(&tail, &head, j, sto->nonBDTDNE[i]);
        ToggleKnownEdge(tail, head, sto->nonBDTDNE[i], TRUE);
      }
      
      sto->discordantEdges[i]->nedges = 0;
    }
  }
}

MH_P_FN(MH_discordBDStratTNT) {
  GET_STORAGE(discordBDStratTNTStorage, sto);
  GET_AUX_STORAGE(StoreTimeAndLasttoggle, dur_inf);

  double ur = unif_rand();
  
  // find the first mixing type strat_i with (cumulative) probability larger than ur
  int strat_i = 0;
  while(ur > sto->currentprobvec[strat_i]) {
    ur -= sto->currentprobvec[strat_i];
    strat_i++;
  }

  // at this point strat mixing type strat_i must have a toggleable dyad

  // record the mixing type of the toggle, in case it's needed in the U function later
  sto->stratmixingtype = strat_i;    

  int strattailtype = sto->strattailtypes[strat_i];
  int stratheadtype = sto->stratheadtypes[strat_i];
    
  // number of edges of this mixing type
  int nedgestype = sto->nonDiscordantEdges[strat_i]->nedges + sto->discordantEdges[strat_i]->nedges;
  int ndyadstype = sto->currentdyads[strat_i];
  
  int nddyadstype = sto->discordantEdges[strat_i]->nedges + EDGECOUNT(sto->BDTDNE[strat_i]);
  
  int in_network;
  int in_discord;
  
  if(nddyadstype == 0 || unif_rand() < 0.5) {
    // propose from network
    if((unif_rand() < 0.5 && nedgestype > 0) || ndyadstype == 0) {
      // propose toggling off an existing edge of strat mixing type strat_i
      if(unif_rand() < (((double) sto->nonDiscordantEdges[strat_i]->nedges) / nedgestype)) {
        UnsrtELGetRand(Mtail, Mhead, sto->nonDiscordantEdges[strat_i]);
        in_discord = FALSE;
      } else {
        UnsrtELGetRand(Mtail, Mhead, sto->discordantEdges[strat_i]);
        in_discord = TRUE;
      }
      in_network = TRUE;
    } else {
      // select a random BD toggleable dyad of strat mixing type strat_i and propose toggling it
      int dyadindex = 2*sto->currentdyads[strat_i]*unif_rand();
  
      Vertex head;
      Vertex tail;
  
      double *BDtailsbyStrattype = sto->BDtailsbyStrattype;
      double *BDheadsbyStrattype = sto->BDheadsbyStrattype;
  
      // need to increment BDtails/headsbyStrattype    
      for(int j = 0; j < strat_i; j++) {
        BDtailsbyStrattype += (int)sto->BDtypesbyStrattype[j];
        BDheadsbyStrattype += (int)sto->BDtypesbyStrattype[j];
      }
      
      // this rather ugly block of code is just finding the dyad that corresponds
      // to the dyadindex we drew above, and then setting the info for
      // tail and head appropriately
      for(int j = 0; j < sto->BDtypesbyStrattype[strat_i]; j++) {
        int dyadstype;
  
        int tailcounts = sto->attrcounts[(int)sto->strattailtypes[strat_i]][(int)BDtailsbyStrattype[j]];
        int headcounts = sto->attrcounts[(int)sto->stratheadtypes[strat_i]][(int)BDheadsbyStrattype[j]];
  
        if(sto->strattailtypes[strat_i] != sto->stratheadtypes[strat_i] || BDtailsbyStrattype[j] != BDheadsbyStrattype[j]) {
          dyadstype = tailcounts*headcounts;
        } else {
          dyadstype = tailcounts*(headcounts - 1)/2;
        }
        
        if(dyadindex < 2*dyadstype) {
          int tailindex;
          int headindex;
          
          if(sto->strattailtypes[strat_i] == sto->stratheadtypes[strat_i] && BDtailsbyStrattype[j] == BDheadsbyStrattype[j]) {
            tailindex = dyadindex / tailcounts;
            headindex = dyadindex % (headcounts - 1);
            if(tailindex == headindex) {
              headindex = headcounts - 1;
            }
                      
            tail = sto->nodesvec[strattailtype][(int)BDtailsbyStrattype[j]][tailindex];
            head = sto->nodesvec[stratheadtype][(int)BDheadsbyStrattype[j]][headindex];
          } else {
            dyadindex /= 2;
            tailindex = dyadindex / headcounts;
            headindex = dyadindex % headcounts;
            
            tail = sto->nodesvec[strattailtype][(int)BDtailsbyStrattype[j]][tailindex];
            head = sto->nodesvec[stratheadtype][(int)BDheadsbyStrattype[j]][headindex];
          }
              
          if(tail > head) {
            Mtail[0] = head;
            Mhead[0] = tail;
          } else {
            Mtail[0] = tail;
            Mhead[0] = head;
          }
          
          break;
        } else {
          dyadindex -= 2*dyadstype;
        }
      }
      
      
      in_discord = kh_get(DyadMapInt, dur_inf->discord, THKey(dur_inf->discord, Mtail[0], Mhead[0])) != kh_none;    
      in_network = IS_OUTEDGE(Mtail[0],Mhead[0]);
      // now check if the dyad we drew is already an edge or not
      if(in_network) {
        // must resample to know edge index; we will fix strat and bd types below;
        // indices within submaximal lists will not be used in the U_FN since this is an off-toggle,
        // so if they're wrong, that's fine
        if(in_discord) {
          UnsrtELGetRand(Mtail, Mhead, sto->discordantEdges[strat_i]);
        } else {
          UnsrtELGetRand(Mtail, Mhead, sto->nonDiscordantEdges[strat_i]);
        }
      }
    }
  } else {
    // propose from discord
    if(unif_rand() < ((double) sto->discordantEdges[strat_i]->nedges)/nddyadstype) {
      UnsrtELGetRand(Mtail, Mhead, sto->discordantEdges[strat_i]);
      in_network = TRUE;
    } else {
      GetRandEdge(Mtail, Mhead, sto->BDTDNE[strat_i]);
      in_network = FALSE;
    }
    
    in_discord = TRUE;
  }

  sto->in_discord = in_discord;

  sto->strattailtype = sto->strat_vattr[Mtail[0] - 1];
  sto->stratheadtype = sto->strat_vattr[Mhead[0] - 1];
    
  sto->bdtailtype = sto->bd_vattr[Mtail[0] - 1];
  sto->bdheadtype = sto->bd_vattr[Mhead[0] - 1];

  if(in_network) {    
    sto->tailmaxl = IN_DEG[Mtail[0]] + OUT_DEG[Mtail[0]] == sto->bound;
    sto->headmaxl = IN_DEG[Mhead[0]] + OUT_DEG[Mhead[0]] == sto->bound;
  } else {
    sto->tailmaxl = IN_DEG[Mtail[0]] + OUT_DEG[Mtail[0]] == sto->bound - 1;
    sto->headmaxl = IN_DEG[Mhead[0]] + OUT_DEG[Mhead[0]] == sto->bound - 1;
  }
  
  // next need to figure out which mixing types will have > 0 toggleable dyads in the proposed network
  // this can be due to having toggleable edges or nonedges, discordant or nondiscordant
  // note however that anything toggleable in discord must also be toggleable through the ``network'' branch of the above code
  // so e.g. we don't need to separately consider BDTDNEs because those must show up in the proposeddyads calculation
  
  double proposedcumprob = 0;
  
  double *BDtailsbyStrattype = sto->BDtailsbyStrattype;
  double *BDheadsbyStrattype = sto->BDheadsbyStrattype;
  
  for(int i = 0; i < sto->nmixtypes; i++) {
    sto->proposeddyads[i] = 0;        
    for(int j = 0; j < sto->BDtypesbyStrattype[i]; j++) {
      // adjustments
      int proposedtailadjustment = (sto->strattailtype == sto->strattailtypes[i] && sto->bdtailtype == BDtailsbyStrattype[j] && sto->tailmaxl) + (sto->stratheadtype == sto->strattailtypes[i] && sto->bdheadtype == BDtailsbyStrattype[j] && sto->headmaxl);
      int proposedheadadjustment = (sto->strattailtype == sto->stratheadtypes[i] && sto->bdtailtype == BDheadsbyStrattype[j] && sto->tailmaxl) + (sto->stratheadtype == sto->stratheadtypes[i] && sto->bdheadtype == BDheadsbyStrattype[j] && sto->headmaxl);
      
      int tailcounts = sto->attrcounts[(int)sto->strattailtypes[i]][(int)BDtailsbyStrattype[j]];
      int headcounts = sto->attrcounts[(int)sto->stratheadtypes[i]][(int)BDheadsbyStrattype[j]];
      
      if(!in_network) {
        proposedtailadjustment = -proposedtailadjustment;
        proposedheadadjustment = -proposedheadadjustment;
      }
      
      if(sto->strattailtypes[i] != sto->stratheadtypes[i] || BDtailsbyStrattype[j] != BDheadsbyStrattype[j]) {
        sto->proposeddyads[i] += (tailcounts + proposedtailadjustment)*(headcounts + proposedheadadjustment);
      } else {
        sto->proposeddyads[i] += (tailcounts + proposedtailadjustment)*(headcounts + proposedheadadjustment - 1)/2;
      }
    }
  
    BDtailsbyStrattype += (int)sto->BDtypesbyStrattype[i];
    BDheadsbyStrattype += (int)sto->BDtypesbyStrattype[i];
    
    if(sto->proposeddyads[i] + sto->nonDiscordantEdges[i]->nedges + sto->discordantEdges[i]->nedges + (in_network ? -(i == strat_i) : (i == strat_i)) > 0) {
      sto->proposedprobvec[i] = sto->originalprobvec[i];
      proposedcumprob += sto->originalprobvec[i];
    } else {
      sto->proposedprobvec[i] = 0;
    }
  }

  for(int i = 0; i < sto->nmixtypes; i++) {
    sto->proposedprobvec[i] /= proposedcumprob;
  }
  
  // now we need statistics for the current mixing type so we can compute the logratio, including discordant dyad counts

  // how nddyads can change:
  // the dyad we toggle can add/remove one discordant dyad
  // discordant nonedges can change toggleability status based on the toggle we make,
  // but only if they have Mtail or Mhead as an endpoint
  // so check BDTDNE and/or nonBDTDNE neighbors of Mtail and Mhead, other than Mtail -> Mhead itself
  // which should be taken into account at the outset
  int propnddyadstype = nddyadstype;
  
  // direct effect of the toggle we are proposing
  if(in_discord) {
    propnddyadstype--;
  } else {
    propnddyadstype++;
  }
  
  // indirect effect of causing other discordant nonedges to change toggleability status
  if(!in_network) {
    // may reduce propnddyadstype; subtract in_discord to avoid repeatedly counting the Mtail[0] -> Mhead[0] edge
    if(sto->tailmaxl) {
      propnddyadstype -= sto->BDTDNE[strat_i]->indegree[Mtail[0]] + sto->BDTDNE[strat_i]->outdegree[Mtail[0]] - in_discord;
    }
    
    if(sto->headmaxl) {
      propnddyadstype -= sto->BDTDNE[strat_i]->indegree[Mhead[0]] + sto->BDTDNE[strat_i]->outdegree[Mhead[0]] - in_discord;
    }
  } else {
    // may increase propnddyads, but only if other endpoint is also submaximal
    if(sto->tailmaxl) {
      if(sto->nonBDTDNE[strat_i]->indegree[Mtail[0]] > 0) {
        Edge e = Mtail[0];
        do { // NB the neighbor cannot be head
          if(IN_DEG[sto->nonBDTDNE[strat_i]->inedges[e].value] + OUT_DEG[sto->nonBDTDNE[strat_i]->inedges[e].value] < sto->bound) {
            propnddyadstype++;
          }
        } while((e = EdgetreePreSuccessor(sto->nonBDTDNE[strat_i]->inedges, e)) != 0);       
      }
      
      if(sto->nonBDTDNE[strat_i]->outdegree[Mtail[0]] > 0) {
        Edge e = Mtail[0];
        do { // NB the neighbor cannot be head
          if(IN_DEG[sto->nonBDTDNE[strat_i]->outedges[e].value] + OUT_DEG[sto->nonBDTDNE[strat_i]->outedges[e].value] < sto->bound) {
            propnddyadstype++;
          }
        } while((e = EdgetreePreSuccessor(sto->nonBDTDNE[strat_i]->outedges, e)) != 0);               
      }
    }
    
    if(sto->headmaxl) {
      if(sto->nonBDTDNE[strat_i]->indegree[Mhead[0]] > 0) {
        Edge e = Mhead[0];
        do { // NB the neighbor cannot be head
          if(IN_DEG[sto->nonBDTDNE[strat_i]->inedges[e].value] + OUT_DEG[sto->nonBDTDNE[strat_i]->inedges[e].value] < sto->bound) {
            propnddyadstype++;
          }
        } while((e = EdgetreePreSuccessor(sto->nonBDTDNE[strat_i]->inedges, e)) != 0);       
      }
      
      if(sto->nonBDTDNE[strat_i]->outdegree[Mhead[0]] > 0) {
        Edge e = Mhead[0];
        do { // NB the neighbor cannot be head
          if(IN_DEG[sto->nonBDTDNE[strat_i]->outedges[e].value] + OUT_DEG[sto->nonBDTDNE[strat_i]->outedges[e].value] < sto->bound) {
            propnddyadstype++;
          }
        } while((e = EdgetreePreSuccessor(sto->nonBDTDNE[strat_i]->outedges, e)) != 0);               
      }
    }
  }

  double forward_weight = sto->currentprobvec[strat_i];
  double backward_weight = sto->proposedprobvec[strat_i];

  double forward_network = in_network ? ((sto->currentdyads[strat_i] == 0 ? 1.0 : 0.5)/nedgestype + (sto->tailmaxl || sto->headmaxl ? 0 : 0.5/sto->currentdyads[strat_i])) : ((nedgestype == 0 ? 1.0 : 0.5)/sto->currentdyads[strat_i]);
  
  double forward_discord = in_discord ? 1.0/nddyadstype : 0;
  
  double backward_network = in_network ? ((nedgestype == 1 ? 1.0 : 0.5)/sto->proposeddyads[strat_i]) : ((sto->proposeddyads[strat_i] == 0 ? 1.0 : 0.5)/(nedgestype + 1) + (sto->tailmaxl || sto->headmaxl ? 0 : 0.5/sto->proposeddyads[strat_i]));
  
  double backward_discord = in_discord ? 0 : 1.0/propnddyadstype;
  
  if(nddyadstype == 0) forward_network *= 2;
  if(propnddyadstype == 0) backward_network *= 2;
  
  double forward = forward_discord + forward_network;
  double backward = backward_discord + backward_network;

  MHp->logratio = log((backward_weight*backward)/(forward_weight*forward));
}

MH_U_FN(Mu_discordBDStratTNT) {
  GET_STORAGE(discordBDStratTNTStorage, sto);

  // add or remove the dyad being toggled from the relevant edge set
  if(edgeflag) {
    if(sto->in_discord) {
      // discordant edge
      UnsrtELDelete(tail, head, sto->discordantEdges[sto->stratmixingtype]);
    } else {
      // nondiscordant edge; will become BDTDNE
      UnsrtELDelete(tail, head, sto->nonDiscordantEdges[sto->stratmixingtype]);      
      ToggleKnownEdge(tail, head, sto->BDTDNE[sto->stratmixingtype], FALSE);
    }
  } else {
    if(sto->in_discord) {
      // discordant non-edge; evidently BD toggleable
      ToggleKnownEdge(tail, head, sto->BDTDNE[sto->stratmixingtype], TRUE); // it's an edge in the (BDT) discordant non-edge network :)    
      UnsrtELInsert(tail, head, sto->nonDiscordantEdges[sto->stratmixingtype]);      
    } else {
      // nondiscordant nonedge; will become discordantEdge
      UnsrtELInsert(tail, head, sto->discordantEdges[sto->stratmixingtype]);        
    }
  }
  
  // update dyad counts and prob vec
  memcpy(sto->currentdyads, sto->proposeddyads, sto->nmixtypes*sizeof(int));
  memcpy(sto->currentprobvec, sto->proposedprobvec, sto->nmixtypes*sizeof(double));
  
  
  // if (sub)maximality/BD toggleability has changed for any nodes/dyads, update storage accordingly
  if(edgeflag) {
    if(sto->tailmaxl) {
      // tail will be newly submaxl after toggle, so add it to the appropriate node list
      sto->nodesvec[sto->strattailtype][sto->bdtailtype][sto->attrcounts[sto->strattailtype][sto->bdtailtype]] = tail;
      sto->nodepos[tail - 1] = sto->attrcounts[sto->strattailtype][sto->bdtailtype];
      sto->attrcounts[sto->strattailtype][sto->bdtailtype]++;
      
      
      for(int i = 0; i < sto->nmixtypes; i++) {
        sto->transferEL->nedges = 0;
        
        // iterate over all nonBDTDNEs with tail as an endpoint; if other endpoint is also
        // submaxl, move this edge from nonBDTDNE to BDTDNE
        if(sto->nonBDTDNE[i]->indegree[tail] > 0) {
          Edge e = tail;
          do { // NB the neighbor cannot be head
            if(IN_DEG[sto->nonBDTDNE[i]->inedges[e].value] + OUT_DEG[sto->nonBDTDNE[i]->inedges[e].value] < sto->bound) {
              UnsrtELInsert(sto->nonBDTDNE[i]->inedges[e].value, tail, sto->transferEL);
            }
          } while((e = EdgetreePreSuccessor(sto->nonBDTDNE[i]->inedges, e)) != 0);
        }
        
        if(sto->nonBDTDNE[i]->outdegree[tail] > 0) {
          Edge e = tail;
          do {
            if(IN_DEG[sto->nonBDTDNE[i]->outedges[e].value] + OUT_DEG[sto->nonBDTDNE[i]->outedges[e].value] < sto->bound) {
              UnsrtELInsert(tail, sto->nonBDTDNE[i]->outedges[e].value, sto->transferEL);
            }
          } while((e = EdgetreePreSuccessor(sto->nonBDTDNE[i]->outedges, e)) != 0);
        }                
  
        for(int j = 1; j <= sto->transferEL->nedges; j++) {
          ToggleKnownEdge(sto->transferEL->tails[j], sto->transferEL->heads[j], sto->nonBDTDNE[i], TRUE);            
          ToggleKnownEdge(sto->transferEL->tails[j], sto->transferEL->heads[j], sto->BDTDNE[i], FALSE);
        }
      }
    }
    
    if(sto->headmaxl) {
      // head will be newly submaxl after toggle, so add it to the appropriate node list    
      sto->nodesvec[sto->stratheadtype][sto->bdheadtype][sto->attrcounts[sto->stratheadtype][sto->bdheadtype]] = head;
      sto->nodepos[head - 1] = sto->attrcounts[sto->stratheadtype][sto->bdheadtype];
      sto->attrcounts[sto->stratheadtype][sto->bdheadtype]++;

      for(int i = 0; i < sto->nmixtypes; i++) {
        sto->transferEL->nedges = 0;
  
        // iterate over all nonBDTDNEs with head as an endpoint; if other endpoint is also
        // submaxl, move this edge from nonBDTDNE to BDTDNE
        if(sto->nonBDTDNE[i]->indegree[head] > 0) {
          Edge e = head;
          do { // NB the neighbor cannot be tail
            if(IN_DEG[sto->nonBDTDNE[i]->inedges[e].value] + OUT_DEG[sto->nonBDTDNE[i]->inedges[e].value] < sto->bound) {
              UnsrtELInsert(sto->nonBDTDNE[i]->inedges[e].value, head, sto->transferEL);
            }
          } while((e = EdgetreePreSuccessor(sto->nonBDTDNE[i]->inedges, e)) != 0);
        }
        
        if(sto->nonBDTDNE[i]->outdegree[head] > 0) {
          Edge e = head;
          do {
            if(IN_DEG[sto->nonBDTDNE[i]->outedges[e].value] + OUT_DEG[sto->nonBDTDNE[i]->outedges[e].value] < sto->bound) {
              UnsrtELInsert(head, sto->nonBDTDNE[i]->outedges[e].value, sto->transferEL);
            }
          } while((e = EdgetreePreSuccessor(sto->nonBDTDNE[i]->outedges, e)) != 0);
        }                
  
        for(int j = 1; j <= sto->transferEL->nedges; j++) {
          ToggleKnownEdge(sto->transferEL->tails[j], sto->transferEL->heads[j], sto->nonBDTDNE[i], TRUE);            
          ToggleKnownEdge(sto->transferEL->tails[j], sto->transferEL->heads[j], sto->BDTDNE[i], FALSE);
        }
      }
    }
  } else {
    if(sto->tailmaxl) {
      // tail will be newly maxl after toggle, so remove it from the appropriate node list
      sto->nodesvec[sto->strattailtype][sto->bdtailtype][sto->nodepos[tail - 1]] = sto->nodesvec[sto->strattailtype][sto->bdtailtype][sto->attrcounts[sto->strattailtype][sto->bdtailtype] - 1];
      sto->nodepos[sto->nodesvec[sto->strattailtype][sto->bdtailtype][sto->nodepos[tail - 1]]] = sto->nodepos[tail - 1];
      sto->attrcounts[sto->strattailtype][sto->bdtailtype]--;

      for(int i = 0; i < sto->nmixtypes; i++) {
        sto->transferEL->nedges = 0;
        
        // transfer all BDTDNEs with tail as an endpoint to nonBDTDNE
        if(sto->BDTDNE[i]->indegree[tail] > 0) {
          Edge e = tail;
          do {
            UnsrtELInsert(sto->BDTDNE[i]->inedges[e].value, tail, sto->transferEL);
          } while((e = EdgetreePreSuccessor(sto->BDTDNE[i]->inedges, e)) != 0);
        }
        
        if(sto->BDTDNE[i]->outdegree[tail] > 0) {
          Edge e = tail;
          do {
            UnsrtELInsert(sto->BDTDNE[i]->outedges[e].value, tail, sto->transferEL);
          } while((e = EdgetreePreSuccessor(sto->BDTDNE[i]->outedges, e)) != 0);
        }                
  
        for(int j = 1; j <= sto->transferEL->nedges; j++) {
          ToggleKnownEdge(sto->transferEL->tails[j], sto->transferEL->heads[j], sto->BDTDNE[i], TRUE);
          ToggleKnownEdge(sto->transferEL->tails[j], sto->transferEL->heads[j], sto->nonBDTDNE[i], FALSE);            
        }
      }
    }
    
    if(sto->headmaxl) {
      // head will be newly maxl after toggle, so remove it from the appropriate node list
      sto->nodesvec[sto->stratheadtype][sto->bdheadtype][sto->nodepos[head - 1]] = sto->nodesvec[sto->stratheadtype][sto->bdheadtype][sto->attrcounts[sto->stratheadtype][sto->bdheadtype] - 1];
      sto->nodepos[sto->nodesvec[sto->stratheadtype][sto->bdheadtype][sto->nodepos[head - 1]]] = sto->nodepos[head - 1];
      sto->attrcounts[sto->stratheadtype][sto->bdheadtype]--;

      for(int i = 0; i < sto->nmixtypes; i++) {
        sto->transferEL->nedges = 0;
        
        // transfer all BDTDNEs with head as an endpoint to nonBDTDNE
        if(sto->BDTDNE[i]->indegree[head] > 0) {
          Edge e = head;
          do {
            UnsrtELInsert(sto->BDTDNE[i]->inedges[e].value, head, sto->transferEL);
          } while((e = EdgetreePreSuccessor(sto->BDTDNE[i]->inedges, e)) != 0);
        }
        
        if(sto->BDTDNE[i]->outdegree[head] > 0) {
          Edge e = head;
          do {
            UnsrtELInsert(sto->BDTDNE[i]->outedges[e].value, head, sto->transferEL);
          } while((e = EdgetreePreSuccessor(sto->BDTDNE[i]->outedges, e)) != 0);
        }                
  
        for(int j = 1; j <= sto->transferEL->nedges; j++) {
          ToggleKnownEdge(sto->transferEL->tails[j], sto->transferEL->heads[j], sto->BDTDNE[i], TRUE);
          ToggleKnownEdge(sto->transferEL->tails[j], sto->transferEL->heads[j], sto->nonBDTDNE[i], FALSE);            
        }
      }
    }      
  }
}

MH_F_FN(Mf_discordBDStratTNT) {
  // Free all the things
  GET_STORAGE(discordBDStratTNTStorage, sto);

  int nattrcodes = MHp->inputs[1 + 3*sto->nmixtypes];
  int npairings = MHp->inputs[1 + 3*sto->nmixtypes + 1 + N_NODES + nattrcodes*nattrcodes];
  int bdlevels = MHp->inputs[1 + 3*sto->nmixtypes + 1 + N_NODES + nattrcodes*nattrcodes + 1 + npairings + 1];

  for(int i = 0; i < nattrcodes; i++) {
    Free(sto->attrcounts[i]);
  }
  Free(sto->attrcounts);

  for(int i = 0; i < nattrcodes; i++) {
    for(int j = 0; j < bdlevels; j++) {
      Free(sto->nodesvec[i][j]);
    }
    Free(sto->nodesvec[i]);
  }
  Free(sto->nodesvec);
  
  // free dyad counts and probvecs
  Free(sto->proposeddyads);
  Free(sto->currentdyads);

  Free(sto->proposedprobvec);
  Free(sto->currentprobvec);
  
  Free(sto->nodepos);
  
  UnsrtELDestroy(sto->transferEL);
  
  for(int i = 0; i < sto->nmixtypes; i++) {
    NetworkDestroy(sto->BDTDNE[i]);
    NetworkDestroy(sto->nonBDTDNE[i]);
    
    UnsrtELDestroy(sto->discordantEdges[i]);
    UnsrtELDestroy(sto->nonDiscordantEdges[i]);
  }
  
  Free(sto->BDTDNE);
  Free(sto->nonBDTDNE);
  Free(sto->discordantEdges);
  Free(sto->nonDiscordantEdges);
  // MHp->storage itself should be Freed by MHProposalDestroy
}
