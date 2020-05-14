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
  Dyad *dyadcounts;
  
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
    
  double *inputindmat = MHp->inputs + 1 + 3*nmixtypes + 1 + nattrcodes + N_NODES + N_NODES;  
  
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
  sto->pmat = Calloc(nmixtypes, double);

  int empirical_flag = MHp->inputs[1 + 3*nmixtypes + 1 + nattrcodes + N_NODES + N_NODES + nattrcodes*nattrcodes];
  if(empirical_flag) {
    sto->pmat[0] = sto->nonDiscordantELs[0]->nedges;
    for(int i = 1; i < nmixtypes; i++) {
      sto->pmat[i] = sto->pmat[i - 1] + sto->nonDiscordantELs[i]->nedges;
    }

    // empirical_flag with no edges is an error
    if(sto->pmat[nmixtypes - 1] == 0) {
      MHp->ntoggles = MH_FAILED;
      return;
    }

    for(int i = 0; i < nmixtypes; i++) {
      sto->pmat[i] /= sto->pmat[nmixtypes - 1];
    }
  } else {
    memcpy(sto->pmat, MHp->inputs + 1 + 2*nmixtypes, nmixtypes*sizeof(double));
  }
  
  
  sto->tailtypes = MHp->inputs + 1;
  sto->headtypes = MHp->inputs + 1 + nmixtypes;
  
  sto->dyadcounts = Calloc(nmixtypes, Dyad);
  for(int i = 0; i < nmixtypes; i++) {
    int tailtype = sto->tailtypes[i];
    int headtype = sto->headtypes[i];
    
    int tailcounts = sto->nodecountsbycode[tailtype];
    int headcounts = sto->nodecountsbycode[headtype];
    
    if(tailtype == headtype) {
      if(DIRECTED) {
        sto->dyadcounts[i] = (Dyad)tailcounts*(headcounts - 1);
      } else {
        sto->dyadcounts[i] = (Dyad)tailcounts*(headcounts - 1)/2;
      }
    } else {
      sto->dyadcounts[i] = (Dyad)tailcounts*headcounts;
    }
  }  
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
  Dyad ndyadstype = sto->dyadcounts[i];
  
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

  Free(sto->pmat);
  Free(sto->dyadcounts);

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
  
  Dyad currentdyads;
  Dyad proposeddyads;
  
  int bound;
  int nmixtypes;
  double *vattr;
  double *tailtypes;
  double *headtypes;

  UnsrtEL *BDTDNE;
  UnsrtEL *nonBDTDNE;

  UnsrtEL *nonDiscordantEdges;
  UnsrtEL *discordantEdges;
  
  int *BDTDNE_deg;
  int *nonBDTDNE_deg;
  
  int in_discord;  
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
  Dyad currentdyads = 0;    
  for(int i = 0; i < nmixtypes; i++) {
    if(tailtypes[i] == headtypes[i]) {
      currentdyads += (Dyad)attrcounts[(int)tailtypes[i]]*(attrcounts[(int)headtypes[i]] - 1)/2;
    } else {
      currentdyads += (Dyad)attrcounts[(int)tailtypes[i]]*attrcounts[(int)headtypes[i]];
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

  sto->nonDiscordantEdges = UnsrtELInitialize(0, NULL, NULL, FALSE);
  Vertex head;
  Edge e;
  for(Vertex tail = 1; tail <= N_NODES; tail++) {
    STEP_THROUGH_OUTEDGES(tail, e, head) {
      UnsrtELInsert(tail, head, sto->nonDiscordantEdges);
    }
  }
  
  sto->nodepos = nodepos;
  
  sto->discordantEdges = UnsrtELInitialize(0, NULL, NULL, FALSE);
  
  sto->BDTDNE = UnsrtELInitialize(0, NULL, NULL, FALSE);
  sto->nonBDTDNE = UnsrtELInitialize(0, NULL, NULL, FALSE);
  sto->BDTDNE_deg = Calloc(N_NODES + 1, int);
  sto->nonBDTDNE_deg = Calloc(N_NODES + 1, int);
}

MH_X_FN(Mx_discordBDTNT) {
  if(type == TICK) {
    // clear all the discordance information
    
    GET_STORAGE(discordBDTNTStorage, sto);
        
    sto->BDTDNE->nedges = 0;
    sto->nonBDTDNE->nedges = 0;
    
    memset(sto->BDTDNE_deg, 0, (N_NODES + 1)*sizeof(int));
    memset(sto->nonBDTDNE_deg, 0, (N_NODES + 1)*sizeof(int));
    
    for(int i = 1; i <= sto->discordantEdges->nedges; i++) {
      UnsrtELInsert(sto->discordantEdges->tails[i], sto->discordantEdges->heads[i], sto->nonDiscordantEdges);
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
  
  int nddyads = sto->discordantEdges->nedges + sto->BDTDNE->nedges;
  
  if(nddyads == 0 || unif_rand() < 0.5) {
    // propose from network
    
    // if currentdyads == 0, we *must* propose toggling off an existing edge;
    // the case nedges == 0 && currentdyads == 0 was excluded during initialization,
    // and we cannot end up in that case if we don't start in that case
    // (assuming the initial network is valid)  
    if((unif_rand() < 0.5 && nedges > 0) || (sto->currentdyads == 0)) {
      // select an existing edge at random, and propose toggling it off
      if(unif_rand() < ((double) sto->nonDiscordantEdges->nedges)/nedges) {
        UnsrtELGetRand(Mtail, Mhead, sto->nonDiscordantEdges);
        in_discord = FALSE;
      } else {
        UnsrtELGetRand(Mtail, Mhead, sto->discordantEdges);
        in_discord = TRUE;
      }
          
      in_network = TRUE;
    } else {
      // select a BD-toggleable dyad and propose toggling it
      // doubling here allows more efficient calculation in 
      // the case tailtypes[i] == headtypes[i]
      
      Dyad dyadindex = 2*sto->currentdyads*unif_rand();
              
      Vertex head;
      Vertex tail;
      
      // this rather ugly block of code is just finding the dyad that corresponds
      // to the dyadindex we drew above, and then setting the info for
      // tail and head appropriately
      for(int i = 0; i < sto->nmixtypes; i++) {
        Dyad dyadstype;
        if(sto->tailtypes[i] == sto->headtypes[i]) {
          dyadstype = (Dyad)sto->attrcounts[(int)sto->tailtypes[i]]*(sto->attrcounts[(int)sto->headtypes[i]] - 1)/2;
        } else {
          dyadstype = (Dyad)sto->attrcounts[(int)sto->tailtypes[i]]*sto->attrcounts[(int)sto->headtypes[i]];
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
      in_discord = kh_get(DyadMapInt, dur_inf->discord, THKey(dur_inf->discord, Mtail[0], Mhead[0])) != kh_none;

      if(in_discord) {
        // need to resample
        if(in_network) {
          UnsrtELGetRand(Mtail, Mhead, sto->discordantEdges);
        } else {
          UnsrtELGetRand(Mtail, Mhead, sto->BDTDNE);
        }
      } else if(in_network) {
        UnsrtELGetRand(Mtail, Mhead, sto->nonDiscordantEdges);
      }
    }
  } else {
    // propose from discord      
    if(unif_rand() < ((double) sto->discordantEdges->nedges)/nddyads) {
      UnsrtELGetRand(Mtail, Mhead, sto->discordantEdges);
      in_network = TRUE;
    } else {
      UnsrtELGetRand(Mtail, Mhead, sto->BDTDNE);
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
  sto->proposeddyads = sto->currentdyads;
  
  int delta = in_network ? +1 : -1;
  
  for(int i = 0; i < sto->nmixtypes; i++) {
    int corr = 0;
    int ha = 0;

    if(sto->tailtype == sto->headtypes[i] && sto->tailmaxl) {
      ha += delta;
      corr += sto->attrcounts[(int)sto->tailtypes[i]];
    }
      
    if(sto->headtype == sto->headtypes[i] && sto->headmaxl) {
      ha += delta;
      corr += sto->attrcounts[(int)sto->tailtypes[i]];        
    }
      
    if(sto->tailtype == sto->tailtypes[i] && sto->tailmaxl) {
      corr += sto->attrcounts[(int)sto->headtypes[i]] + ha;
    }
      
    if(sto->headtype == sto->tailtypes[i] && sto->headmaxl) {
      corr += sto->attrcounts[(int)sto->headtypes[i]] + ha;
    }
      
    if(sto->tailtypes[i] == sto->headtypes[i]) {
      if(in_network) {
        corr -= ha;
      } else {
        corr += ha;
      }
      corr /= 2;
    }
    
    if(in_network) {
      sto->proposeddyads += corr;
    } else {
      sto->proposeddyads -= corr;
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
      propnddyads -= sto->BDTDNE_deg[Mtail[0]] - in_discord;
    }
    
    if(sto->headmaxl) {
      propnddyads -= sto->BDTDNE_deg[Mhead[0]] - in_discord;
    }
  } else {
    // may increase propnddyads, but only if other endpoint is also submaximal
    if(sto->tailmaxl) {
      if(sto->nonBDTDNE_deg[Mtail[0]] > 0) {
        for(int i = 1; i <= sto->nonBDTDNE->nedges; i++) {
          if((sto->nonBDTDNE->tails[i] == Mtail[0] && IN_DEG[sto->nonBDTDNE->heads[i]] + OUT_DEG[sto->nonBDTDNE->heads[i]] < sto->bound) ||
             (sto->nonBDTDNE->heads[i] == Mtail[0] && IN_DEG[sto->nonBDTDNE->tails[i]] + OUT_DEG[sto->nonBDTDNE->tails[i]] < sto->bound)) {
            propnddyads++;
          }
        }
      }
    }
    
    if(sto->headmaxl) {
      if(sto->nonBDTDNE_deg[Mhead[0]] > 0) {
        for(int i = 1; i <= sto->nonBDTDNE->nedges; i++) {
          if((sto->nonBDTDNE->tails[i] == Mhead[0] && IN_DEG[sto->nonBDTDNE->heads[i]] + OUT_DEG[sto->nonBDTDNE->heads[i]] < sto->bound) ||
             (sto->nonBDTDNE->heads[i] == Mhead[0] && IN_DEG[sto->nonBDTDNE->tails[i]] + OUT_DEG[sto->nonBDTDNE->tails[i]] < sto->bound)) {
            propnddyads++;
          }
        }
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
      UnsrtELDelete(tail, head, sto->nonDiscordantEdges);
      UnsrtELInsert(tail, head, sto->BDTDNE);

      sto->BDTDNE_deg[tail]++;
      sto->BDTDNE_deg[head]++;
    }
  } else {
    if(sto->in_discord) {
      // discordant non-edge; evidently BD toggleable
      UnsrtELDelete(tail, head, sto->BDTDNE);
      UnsrtELInsert(tail, head, sto->nonDiscordantEdges);

      sto->BDTDNE_deg[tail]--;
      sto->BDTDNE_deg[head]--;      
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
      
      // iterate over all nonBDTDNEs with tail as an endpoint; if other endpoint is also
      // submaxl, move this edge from nonBDTDNE to BDTDNE
      if(sto->nonBDTDNE_deg[tail] > 0) {
        for(int j = sto->nonBDTDNE->nedges; j >= 1; j--) {
          Vertex ttail = sto->nonBDTDNE->tails[j];
          Vertex thead = sto->nonBDTDNE->heads[j];
          if((tail == ttail && IN_DEG[thead] + OUT_DEG[thead] < sto->bound) || (tail == thead && IN_DEG[ttail] + OUT_DEG[ttail] < sto->bound)) {
            sto->nonBDTDNE->lindex = j;
            sto->nonBDTDNE->ltail = ttail;
            sto->nonBDTDNE->lhead = thead;
            UnsrtELDelete(ttail, thead, sto->nonBDTDNE);
            UnsrtELInsert(ttail, thead, sto->BDTDNE);
            sto->nonBDTDNE_deg[ttail]--;
            sto->nonBDTDNE_deg[thead]--;
            sto->BDTDNE_deg[ttail]++;
            sto->BDTDNE_deg[thead]++;            
          }
        }
      }
    }
    
    if(sto->headmaxl) {
      // head will be newly submaxl after toggle, so add it to the appropriate node list    
      sto->nodesvec[sto->headtype][sto->attrcounts[sto->headtype]] = head;
      sto->nodepos[head - 1] = sto->attrcounts[sto->headtype];
      sto->attrcounts[sto->headtype]++;

      // iterate over all nonBDTDNEs with head as an endpoint; if other endpoint is also
      // submaxl, move this edge from nonBDTDNE to BDTDNE
      if(sto->nonBDTDNE_deg[head] > 0) {
        for(int j = sto->nonBDTDNE->nedges; j >= 1; j--) {
          Vertex ttail = sto->nonBDTDNE->tails[j];
          Vertex thead = sto->nonBDTDNE->heads[j];
          if((head == ttail && IN_DEG[thead] + OUT_DEG[thead] < sto->bound) || (head == thead && IN_DEG[ttail] + OUT_DEG[ttail] < sto->bound)) {
            sto->nonBDTDNE->lindex = j;
            sto->nonBDTDNE->ltail = ttail;
            sto->nonBDTDNE->lhead = thead;
            UnsrtELDelete(ttail, thead, sto->nonBDTDNE);
            UnsrtELInsert(ttail, thead, sto->BDTDNE);
            sto->nonBDTDNE_deg[ttail]--;
            sto->nonBDTDNE_deg[thead]--;
            sto->BDTDNE_deg[ttail]++;
            sto->BDTDNE_deg[thead]++;            
          }
        }
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

      // transfer all BDTDNEs with tail as an endpoint to nonBDTDNE
      if(sto->BDTDNE_deg[tail] > 0) {
        for(int j = sto->BDTDNE->nedges; j >= 1; j--) {
          Vertex ttail = sto->BDTDNE->tails[j];
          Vertex thead = sto->BDTDNE->heads[j];
          if(tail == ttail || tail == thead) {
            sto->BDTDNE->lindex = j;
            sto->BDTDNE->ltail = ttail;
            sto->BDTDNE->lhead = thead;
            UnsrtELDelete(ttail, thead, sto->BDTDNE);
            UnsrtELInsert(ttail, thead, sto->nonBDTDNE);
            sto->nonBDTDNE_deg[ttail]++;
            sto->nonBDTDNE_deg[thead]++;
            sto->BDTDNE_deg[ttail]--;
            sto->BDTDNE_deg[thead]--;            
          }
        }        
      }          
    }
    
    if(sto->headmaxl) {
      // head will be newly maxl after toggle, so remove it from the appropriate node list, updating nodepos
      // for whatever node is currently at the end of sto->nodesvec[sto->headtype], since that will be
      // moved into head's position
      sto->nodepos[sto->nodesvec[sto->headtype][sto->attrcounts[sto->headtype] - 1] - 1] = sto->nodepos[head - 1];
      sto->nodesvec[sto->headtype][sto->nodepos[head - 1]] = sto->nodesvec[sto->headtype][sto->attrcounts[sto->headtype] - 1];
      sto->attrcounts[sto->headtype]--;

      if(sto->BDTDNE_deg[head] > 0) {
        for(int j = sto->BDTDNE->nedges; j >= 1; j--) {
          Vertex ttail = sto->BDTDNE->tails[j];
          Vertex thead = sto->BDTDNE->heads[j];
          if(head == ttail || head == thead) {
            sto->BDTDNE->lindex = j;
            sto->BDTDNE->ltail = ttail;
            sto->BDTDNE->lhead = thead;
            UnsrtELDelete(ttail, thead, sto->BDTDNE);
            UnsrtELInsert(ttail, thead, sto->nonBDTDNE);
            sto->nonBDTDNE_deg[ttail]++;
            sto->nonBDTDNE_deg[thead]++;
            sto->BDTDNE_deg[ttail]--;
            sto->BDTDNE_deg[thead]--;                        
          }
        }        
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

  UnsrtELDestroy(sto->BDTDNE);
  UnsrtELDestroy(sto->nonBDTDNE);
  UnsrtELDestroy(sto->discordantEdges);
  UnsrtELDestroy(sto->nonDiscordantEdges);
  
  Free(sto->BDTDNE_deg);
  Free(sto->nonBDTDNE_deg);
  
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
  
  double currentcumprob;
  double proposedcumprob;
  
  double *originalprobvec;
  double *currentprobvec;
  double *proposedprobvec;
  
  int bound;
  int nmixtypes;
  
  double *strat_vattr;
  double *bd_vattr;
  
  double *BDtypesbyStrattype;
  double **BDtailsbyStrattype;
  double **BDheadsbyStrattype;
  
  double *strattailtypes;
  double *stratheadtypes;

  int *totBDTDNE_deg;
  int *totnonBDTDNE_deg;

  int **BDTDNE_deg;
  int **nonBDTDNE_deg;

  UnsrtEL **BDTDNE;
  UnsrtEL **nonBDTDNE;
  
  UnsrtEL **discordantEdges;
  UnsrtEL **nonDiscordantEdges;
  
  int in_discord;
} discordBDStratTNTStorage;

MH_I_FN(Mi_discordBDStratTNT) {
  // process the inputs and initialize all the edgelists in storage; set MHp->ntoggles to 1
  MHp->ntoggles = 1;

  ALLOC_STORAGE(1, discordBDStratTNTStorage, sto);
  
  int nmixtypes = MHp->inputs[0];
  
  double *probvec = MHp->inputs + 1 + 2*nmixtypes;
  
  double *strattailattrs = MHp->inputs + 1;
  double *stratheadattrs = MHp->inputs + 1 + nmixtypes;
  
  int nattrcodes = MHp->inputs[1 + 3*nmixtypes];
  
  double *strat_vattr = MHp->inputs + 1 + 3*nmixtypes + 1;
  
  UnsrtEL **els = (UnsrtEL **)Calloc(nmixtypes, UnsrtEL *);
      
  for(int i = 0; i < nmixtypes; i++) {
    els[i] = UnsrtELInitialize(0, NULL, NULL, FALSE);
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
        UnsrtELInsert(tail, head, els[index]);
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
  
  sto->nodepos = nodepos;

  sto->BDtypesbyStrattype = MHp->inputs + 1 + 3*nmixtypes + 1 + N_NODES + nattrcodes*nattrcodes + 1 + npairings + 1 + 1 + 1 + 2*bdmixtypes + N_NODES;
  
  int sumtypes = MHp->inputs[1 + 3*nmixtypes + 1 + N_NODES + nattrcodes*nattrcodes + 1 + npairings + 1 + 1 + 1 + 2*bdmixtypes + N_NODES + nmixtypes];
  
  sto->BDtailsbyStrattype = Calloc(nmixtypes, double *);
  sto->BDheadsbyStrattype = Calloc(nmixtypes, double *);
  
  sto->BDtailsbyStrattype[0] = MHp->inputs + 1 + 3*nmixtypes + 1 + N_NODES + nattrcodes*nattrcodes + 1 + npairings + 1 + 1 + 1 + 2*bdmixtypes + N_NODES + nmixtypes + 1;
  sto->BDheadsbyStrattype[0] = MHp->inputs + 1 + 3*nmixtypes + 1 + N_NODES + nattrcodes*nattrcodes + 1 + npairings + 1 + 1 + 1 + 2*bdmixtypes + N_NODES + nmixtypes + 1 + sumtypes;
  
  for(int i = 1; i < nmixtypes; i++) {
    sto->BDtailsbyStrattype[i] = sto->BDtailsbyStrattype[i - 1] + (int)sto->BDtypesbyStrattype[i - 1];
    sto->BDheadsbyStrattype[i] = sto->BDheadsbyStrattype[i - 1] + (int)sto->BDtypesbyStrattype[i - 1];
  }
  
  int empirical_flag = MHp->inputs[1 + 3*nmixtypes + 1 + N_NODES + nattrcodes*nattrcodes + 1 + npairings + 1 + 1 + 1 + 2*bdmixtypes + N_NODES + nmixtypes + 1 + sumtypes + sumtypes];
  sto->originalprobvec = Calloc(nmixtypes, double);
  if(empirical_flag) {
    int sumedges = 0;
    for(int i = 0; i < nmixtypes; i++) {
      sto->originalprobvec[i] = els[i]->nedges;
      sumedges += els[i]->nedges;
    }
    
    // empirical_flag with no edges is an error
    if(sumedges == 0) {
      MHp->ntoggles = MH_FAILED;
      return;
    }
    
    for(int i = 0; i < nmixtypes; i++) {
      sto->originalprobvec[i] /= sumedges;
    }
  } else {
    memcpy(sto->originalprobvec, probvec, nmixtypes*sizeof(double));
  }
  
  Dyad sumcurrentdyads = 0;
  
  double *currentprobvec = Calloc(nmixtypes, double);
  double sumprobs = 0;
  
  for(int i = 0; i < nmixtypes; i++) {
    Dyad currentdyads = 0;
    for(int j = 0; j < (int)sto->BDtypesbyStrattype[i]; j++) {
      int tailcounts = attrcounts[(int)strattailattrs[i]][(int)sto->BDtailsbyStrattype[i][j]];
      int headcounts = attrcounts[(int)stratheadattrs[i]][(int)sto->BDheadsbyStrattype[i][j]];
      
      if(strattailattrs[i] != stratheadattrs[i] || sto->BDtailsbyStrattype[i][j] != sto->BDheadsbyStrattype[i][j]) {
        currentdyads += (Dyad)tailcounts*headcounts;
      } else {
        currentdyads += (Dyad)tailcounts*(headcounts - 1)/2;
      }
    }
    
    sumcurrentdyads += currentdyads;
    
    if(currentdyads > 0 || els[i]->nedges > 0) {
      currentprobvec[i] = sto->originalprobvec[i];
      sumprobs += sto->originalprobvec[i];
    } // else it's already 0
  }
  
  // if we cannot toggle any edges or dyads, error
  if(EDGECOUNT(nwp) == 0 && sumcurrentdyads == 0) {
    MHp->ntoggles = MH_FAILED;
    return;
  }
  
  sto->nonDiscordantEdges = els;
  sto->nodesvec = nodesvec;
  sto->attrcounts = attrcounts;
    
  sto->currentprobvec = currentprobvec;
  sto->proposedprobvec = Calloc(nmixtypes, double);
  
  sto->currentcumprob = sumprobs;
  
  // NB: the next two lines reference the logic in the P function
  sto->proposedcumprob = 1;
  memcpy(sto->proposedprobvec, sto->originalprobvec, nmixtypes*sizeof(double));
  
  sto->bound = bound;
  sto->nmixtypes = nmixtypes;
  
  sto->strat_vattr = strat_vattr;
  sto->bd_vattr = bd_vattr;
    
  sto->strattailtypes = MHp->inputs + 1;
  sto->stratheadtypes = MHp->inputs + 1 + nmixtypes;

  sto->BDTDNE_deg = Calloc(nmixtypes, int *);
  sto->nonBDTDNE_deg = Calloc(nmixtypes, int *);

  sto->totBDTDNE_deg = Calloc(N_NODES + 1, int);
  sto->totnonBDTDNE_deg = Calloc(N_NODES + 1, int);

  sto->BDTDNE = Calloc(nmixtypes, UnsrtEL *);
  sto->nonBDTDNE = Calloc(nmixtypes, UnsrtEL *);
  sto->discordantEdges = Calloc(nmixtypes, UnsrtEL *);
  
  // also initialize discord stuff
  for(int i = 0; i < nmixtypes; i++) {
    sto->BDTDNE_deg[i] = Calloc(N_NODES + 1, int);
    sto->nonBDTDNE_deg[i] = Calloc(N_NODES + 1, int);
    
    sto->BDTDNE[i] = UnsrtELInitialize(0, NULL, NULL, FALSE);
    sto->nonBDTDNE[i] = UnsrtELInitialize(0, NULL, NULL, FALSE);
    
    sto->discordantEdges[i] = UnsrtELInitialize(0, NULL, NULL, FALSE);
  }
}

MH_X_FN(Mx_discordBDStratTNT) {
  if(type == TICK) {
    GET_STORAGE(discordBDStratTNTStorage, sto);

    memset(sto->totBDTDNE_deg, 0, (N_NODES + 1)*sizeof(int));
    memset(sto->totnonBDTDNE_deg, 0, (N_NODES + 1)*sizeof(int));


    for(int i = 0; i < sto->nmixtypes; i++) {
      // clear all the discordance information

      if(sto->BDTDNE[i]->nedges > 0) {
        memset(sto->BDTDNE_deg[i], 0, (N_NODES + 1)*sizeof(int));
        sto->BDTDNE[i]->nedges = 0;
      }
      
      if(sto->nonBDTDNE[i]->nedges > 0) {
        memset(sto->nonBDTDNE_deg[i], 0, (N_NODES + 1)*sizeof(int));
        sto->nonBDTDNE[i]->nedges = 0;
      }
      
      // transfer discordantEdges to nonDiscordantEdges
      Vertex *tails = sto->discordantEdges[i]->tails;
      Vertex *heads = sto->discordantEdges[i]->heads;
      for(int j = 1; j <= sto->discordantEdges[i]->nedges; j++) {
        UnsrtELInsert(tails[j], heads[j], sto->nonDiscordantEdges[i]);
      }
      sto->discordantEdges[i]->nedges = 0;
    }
  }
}

MH_P_FN(MH_discordBDStratTNT) {
  GET_STORAGE(discordBDStratTNTStorage, sto);
  GET_AUX_STORAGE(StoreTimeAndLasttoggle, dur_inf);

  double ur = unif_rand()*sto->currentcumprob;
  
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
  
  Dyad ndyadstype = 0;
  for(int j = 0; j < (int)sto->BDtypesbyStrattype[strat_i]; j++) {
    if(strattailtype == stratheadtype && sto->BDtailsbyStrattype[strat_i][j] == sto->BDheadsbyStrattype[strat_i][j]) {
      ndyadstype += (Dyad)sto->attrcounts[strattailtype][(int)sto->BDtailsbyStrattype[strat_i][j]]*(sto->attrcounts[stratheadtype][(int)sto->BDheadsbyStrattype[strat_i][j]] - 1)/2;
    } else {
      ndyadstype += (Dyad)sto->attrcounts[strattailtype][(int)sto->BDtailsbyStrattype[strat_i][j]]*sto->attrcounts[stratheadtype][(int)sto->BDheadsbyStrattype[strat_i][j]];
    }
  }
  
  int nddyadstype = sto->discordantEdges[strat_i]->nedges + sto->BDTDNE[strat_i]->nedges;
  
  int in_network;
  int in_discord;
  
  if(nddyadstype == 0 || unif_rand() < 0.5) {
    // propose from network
    if((unif_rand() < 0.5 && nedgestype > 0) || ndyadstype == 0) {
      // propose toggling off an existing edge of strat mixing type strat_i
      if(unif_rand() < ((double) sto->nonDiscordantEdges[strat_i]->nedges)/nedgestype) {
        UnsrtELGetRand(Mtail, Mhead, sto->nonDiscordantEdges[strat_i]);
        in_discord = FALSE;
      } else {
        UnsrtELGetRand(Mtail, Mhead, sto->discordantEdges[strat_i]);
        in_discord = TRUE;
      }
      in_network = TRUE;
    } else {
      // select a random BD toggleable dyad of strat mixing type strat_i and propose toggling it
      Dyad dyadindex = 2*ndyadstype*unif_rand();
  
      Vertex head;
      Vertex tail;
      
      // this rather ugly block of code is just finding the dyad that corresponds
      // to the dyadindex we drew above, and then setting the info for
      // tail and head appropriately
      for(int j = 0; j < (int)sto->BDtypesbyStrattype[strat_i]; j++) {
        Dyad dyadsthistype;
  
        int tailcounts = sto->attrcounts[strattailtype][(int)sto->BDtailsbyStrattype[strat_i][j]];
        int headcounts = sto->attrcounts[stratheadtype][(int)sto->BDheadsbyStrattype[strat_i][j]];
  
        if(strattailtype != stratheadtype || sto->BDtailsbyStrattype[strat_i][j] != sto->BDheadsbyStrattype[strat_i][j]) {
          dyadsthistype = (Dyad)tailcounts*headcounts;
        } else {
          dyadsthistype = (Dyad)tailcounts*(headcounts - 1)/2;
        }
        
        if(dyadindex < 2*dyadsthistype) {
          int tailindex;
          int headindex;
          
          if(strattailtype == stratheadtype && sto->BDtailsbyStrattype[strat_i][j] == sto->BDheadsbyStrattype[strat_i][j]) {
            tailindex = dyadindex / tailcounts;
            headindex = dyadindex % (headcounts - 1);
            if(tailindex == headindex) {
              headindex = headcounts - 1;
            }
                      
            tail = sto->nodesvec[strattailtype][(int)sto->BDtailsbyStrattype[strat_i][j]][tailindex];
            head = sto->nodesvec[stratheadtype][(int)sto->BDheadsbyStrattype[strat_i][j]][headindex];
          } else {
            dyadindex /= 2;
            tailindex = dyadindex / headcounts;
            headindex = dyadindex % headcounts;
            
            tail = sto->nodesvec[strattailtype][(int)sto->BDtailsbyStrattype[strat_i][j]][tailindex];
            head = sto->nodesvec[stratheadtype][(int)sto->BDheadsbyStrattype[strat_i][j]][headindex];
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
          dyadindex -= 2*dyadsthistype;
        }
      }
         
      in_network = IS_OUTEDGE(Mtail[0],Mhead[0]);
      in_discord = kh_get(DyadMapInt, dur_inf->discord, THKey(dur_inf->discord, Mtail[0], Mhead[0])) != kh_none;

      if(in_discord) {
        // need to resample
        if(in_network) {
          UnsrtELGetRand(Mtail, Mhead, sto->discordantEdges[strat_i]);
        } else {
          UnsrtELGetRand(Mtail, Mhead, sto->BDTDNE[strat_i]);
        }
      } else if(in_network) {
        UnsrtELGetRand(Mtail, Mhead, sto->nonDiscordantEdges[strat_i]);
      }
    }
  } else {
    // propose from discord
    if(unif_rand() < ((double) sto->discordantEdges[strat_i]->nedges)/nddyadstype) {
      UnsrtELGetRand(Mtail, Mhead, sto->discordantEdges[strat_i]);
      in_network = TRUE;
    } else {
      UnsrtELGetRand(Mtail, Mhead, sto->BDTDNE[strat_i]);
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
    // strat and bd types already set above
    sto->tailmaxl = IN_DEG[Mtail[0]] + OUT_DEG[Mtail[0]] == sto->bound - 1;
    sto->headmaxl = IN_DEG[Mhead[0]] + OUT_DEG[Mhead[0]] == sto->bound - 1;
  }
    
  if(sto->proposedcumprob < 1) {
    sto->proposedcumprob = 1;
    memcpy(sto->proposedprobvec, sto->originalprobvec, sto->nmixtypes*sizeof(double));
  } // else it's already initialized appropriately for us
  
  for(int i = 0; i < sto->nmixtypes; i++) {
    if(sto->nonDiscordantEdges[i]->nedges > 0 || sto->discordantEdges[i]->nedges > 0 || i == strat_i) {
      continue; // something is obviously toggleable in the proposed network
    }
    
    // if no edges and not strat_i, need to check for any toggleable dyads in the proposed network
    // note BDTDNEs are always counted in the ordinary BD toggleable dyads so we don't need to
    // consider those separately
    int anytoggleable = FALSE;
    
    for(int j = 0; j < (int)sto->BDtypesbyStrattype[i]; j++) {
      // adjustments
      int proposedtailadjustment = (sto->strattailtype == sto->strattailtypes[i] && sto->bdtailtype == sto->BDtailsbyStrattype[i][j] && sto->tailmaxl) + (sto->stratheadtype == sto->strattailtypes[i] && sto->bdheadtype == sto->BDtailsbyStrattype[i][j] && sto->headmaxl);
      int proposedheadadjustment = (sto->strattailtype == sto->stratheadtypes[i] && sto->bdtailtype == sto->BDheadsbyStrattype[i][j] && sto->tailmaxl) + (sto->stratheadtype == sto->stratheadtypes[i] && sto->bdheadtype == sto->BDheadsbyStrattype[i][j] && sto->headmaxl);
      
      int tailcounts = sto->attrcounts[(int)sto->strattailtypes[i]][(int)sto->BDtailsbyStrattype[i][j]];
      int headcounts = sto->attrcounts[(int)sto->stratheadtypes[i]][(int)sto->BDheadsbyStrattype[i][j]];
      
      if(in_network) {
        proposedtailadjustment = -proposedtailadjustment;
        proposedheadadjustment = -proposedheadadjustment;
      }
      
      if(tailcounts > proposedtailadjustment && headcounts > proposedheadadjustment + (sto->strattailtypes[i] == sto->stratheadtypes[i] && sto->BDtailsbyStrattype[i][j] == sto->BDheadsbyStrattype[i][j])) {
        anytoggleable = TRUE;
        break;
      }      
    }
  
    // if no toggleable dyads in the proposed network, must update proposed probs accordingly
    if(!anytoggleable) {
      sto->proposedcumprob -= sto->originalprobvec[i];
      sto->proposedprobvec[i] = 0;
    }
  }

  // need to compute proposed dyad count for current mixing type (only)
  Dyad proposeddyadstype = ndyadstype;

  int delta = in_network ? +1 : -1;

  for(int j = 0; j < (int)sto->BDtypesbyStrattype[strat_i]; j++) {
    int corr = 0;
    int ha = 0;

    if(sto->strattailtype == sto->stratheadtypes[strat_i] && sto->bdtailtype == sto->BDheadsbyStrattype[strat_i][j] && sto->tailmaxl) {
      ha += delta;
      corr += sto->attrcounts[(int)sto->strattailtypes[strat_i]][(int)sto->BDtailsbyStrattype[strat_i][j]];
    }
      
    if(sto->stratheadtype == sto->stratheadtypes[strat_i] && sto->bdheadtype == sto->BDheadsbyStrattype[strat_i][j] && sto->headmaxl) {
      ha += delta;
      corr += sto->attrcounts[(int)sto->strattailtypes[strat_i]][(int)sto->BDtailsbyStrattype[strat_i][j]];        
    }
      
    if(sto->strattailtype == sto->strattailtypes[strat_i] && sto->bdtailtype == sto->BDtailsbyStrattype[strat_i][j] && sto->tailmaxl) {
      corr += sto->attrcounts[(int)sto->stratheadtypes[strat_i]][(int)sto->BDheadsbyStrattype[strat_i][j]] + ha;
    }
      
    if(sto->stratheadtype == sto->strattailtypes[strat_i] && sto->bdheadtype == sto->BDtailsbyStrattype[strat_i][j] && sto->headmaxl) {
      corr += sto->attrcounts[(int)sto->stratheadtypes[strat_i]][(int)sto->BDheadsbyStrattype[strat_i][j]] + ha;
    }
      
    if(sto->strattailtypes[strat_i] == sto->stratheadtypes[strat_i] && sto->BDtailsbyStrattype[strat_i][j] == sto->BDheadsbyStrattype[strat_i][j]) {
      if(in_network) {
        corr -= ha;
      } else {
        corr += ha;
      }
      corr /= 2;
    }
    
    if(in_network) {
      proposeddyadstype += corr;
    } else {
      proposeddyadstype -= corr;
    }      
  }
  
  // need propnddyads count
  
  int propnddyadstype = nddyadstype;
  
  // direct effect of the toggle we are proposing
  if(in_discord) {
    propnddyadstype--;
  } else {
    propnddyadstype++;
  }
  
  // indirect effect of causing other discordant nonedges to change toggleability status
  if(!in_network) {
    // may reduce propnddyads; subtract in_discord to avoid repeatedly counting the Mtail[0] -> Mhead[0] edge
    if(sto->tailmaxl) {
      propnddyadstype -= sto->BDTDNE_deg[strat_i][Mtail[0]] - in_discord;
    }
    
    if(sto->headmaxl) {
      propnddyadstype -= sto->BDTDNE_deg[strat_i][Mhead[0]] - in_discord;
    }
  } else {
    // may increase propnddyads, but only if other endpoint is also submaximal
    if(sto->tailmaxl) {
      if(sto->nonBDTDNE_deg[strat_i][Mtail[0]] > 0) {
        for(int i = 1; i <= sto->nonBDTDNE[strat_i]->nedges; i++) {
          if((sto->nonBDTDNE[strat_i]->tails[i] == Mtail[0] && IN_DEG[sto->nonBDTDNE[strat_i]->heads[i]] + OUT_DEG[sto->nonBDTDNE[strat_i]->heads[i]] < sto->bound) ||
             (sto->nonBDTDNE[strat_i]->heads[i] == Mtail[0] && IN_DEG[sto->nonBDTDNE[strat_i]->tails[i]] + OUT_DEG[sto->nonBDTDNE[strat_i]->tails[i]] < sto->bound)) {
            propnddyadstype++;
          }
        }
      }
    }
    
    if(sto->headmaxl) {
      if(sto->nonBDTDNE_deg[strat_i][Mhead[0]] > 0) {
        for(int i = 1; i <= sto->nonBDTDNE[strat_i]->nedges; i++) {
          if((sto->nonBDTDNE[strat_i]->tails[i] == Mhead[0] && IN_DEG[sto->nonBDTDNE[strat_i]->heads[i]] + OUT_DEG[sto->nonBDTDNE[strat_i]->heads[i]] < sto->bound) ||
             (sto->nonBDTDNE[strat_i]->heads[i] == Mhead[0] && IN_DEG[sto->nonBDTDNE[strat_i]->tails[i]] + OUT_DEG[sto->nonBDTDNE[strat_i]->tails[i]] < sto->bound)) {
            propnddyadstype++;
          }
        }
      }
    }
  }
    
  double prob_weight = sto->currentcumprob/sto->proposedcumprob;
  
  double forward_network = in_network ? ((ndyadstype == 0 ? 1.0 : 0.5)/nedgestype + (sto->tailmaxl || sto->headmaxl ? 0 : 0.5/ndyadstype)) : ((nedgestype == 0 ? 1.0 : 0.5)/ndyadstype);
  
  double forward_discord = in_discord ? 1.0/nddyadstype : 0;
  
  double backward_network = in_network ? ((nedgestype == 1 ? 1.0 : 0.5)/proposeddyadstype) : ((proposeddyadstype == 0 ? 1.0 : 0.5)/(nedgestype + 1) + (sto->tailmaxl || sto->headmaxl ? 0 : 0.5/proposeddyadstype));
  
  double backward_discord = in_discord ? 0 : 1.0/propnddyadstype;
  
  if(nddyadstype == 0) forward_network *= 2;
  if(propnddyadstype == 0) backward_network *= 2;
  
  double forward = forward_discord + forward_network;
  double backward = backward_discord + backward_network;

  MHp->logratio = log(prob_weight*backward/forward);
  
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
      UnsrtELInsert(tail, head, sto->BDTDNE[sto->stratmixingtype]);
      sto->BDTDNE_deg[sto->stratmixingtype][tail]++;
      sto->BDTDNE_deg[sto->stratmixingtype][head]++;
      sto->totBDTDNE_deg[tail]++;
      sto->totBDTDNE_deg[head]++;      
    }
  } else {
    if(sto->in_discord) {
      // discordant non-edge; evidently BD toggleable
      UnsrtELDelete(tail, head, sto->BDTDNE[sto->stratmixingtype]);
      UnsrtELInsert(tail, head, sto->nonDiscordantEdges[sto->stratmixingtype]);      
      sto->BDTDNE_deg[sto->stratmixingtype][tail]--;
      sto->BDTDNE_deg[sto->stratmixingtype][head]--;
      sto->totBDTDNE_deg[tail]--;
      sto->totBDTDNE_deg[head]--;      
    } else {
      // nondiscordant nonedge; will become discordantEdge
      UnsrtELInsert(tail, head, sto->discordantEdges[sto->stratmixingtype]);        
    }
  }
  
  // avoid copying in the common case where all strat types are toggleable in both the current and proposed networks
  if(sto->currentcumprob != 1 || sto->proposedcumprob != 1) {
    sto->currentcumprob = sto->proposedcumprob;
    memcpy(sto->currentprobvec, sto->proposedprobvec, sto->nmixtypes*sizeof(double));
  }
  
  // if (sub)maximality/BD toggleability has changed for any nodes/dyads, update storage accordingly
  if(edgeflag) {
    if(sto->tailmaxl) {
      // tail will be newly submaxl after toggle, so add it to the appropriate node list
      sto->nodesvec[sto->strattailtype][sto->bdtailtype][sto->attrcounts[sto->strattailtype][sto->bdtailtype]] = tail;
      sto->nodepos[tail - 1] = sto->attrcounts[sto->strattailtype][sto->bdtailtype];
      sto->attrcounts[sto->strattailtype][sto->bdtailtype]++;
      
      if(sto->totnonBDTDNE_deg[tail] > 0) {
        for(int i = 0; i < sto->nmixtypes; i++) {
          // iterate over all nonBDTDNEs with tail as an endpoint; if other endpoint is also
          // submaxl, move this edge from nonBDTDNE to BDTDNE
          if(sto->nonBDTDNE_deg[i][tail] > 0) {
            for(int j = sto->nonBDTDNE[i]->nedges; j >= 1; j--) {
              Vertex ttail = sto->nonBDTDNE[i]->tails[j];
              Vertex thead = sto->nonBDTDNE[i]->heads[j];
              if((tail == ttail && IN_DEG[thead] + OUT_DEG[thead] < sto->bound) || (tail == thead && IN_DEG[ttail] + OUT_DEG[ttail] < sto->bound)) {
                sto->nonBDTDNE[i]->lindex = j;
                sto->nonBDTDNE[i]->ltail = ttail;
                sto->nonBDTDNE[i]->lhead = thead;
                UnsrtELDelete(ttail, thead, sto->nonBDTDNE[i]);
                UnsrtELInsert(ttail, thead, sto->BDTDNE[i]);
                sto->nonBDTDNE_deg[i][ttail]--;
                sto->nonBDTDNE_deg[i][thead]--;
                sto->BDTDNE_deg[i][ttail]++;
                sto->BDTDNE_deg[i][thead]++;            
  
                sto->totnonBDTDNE_deg[ttail]--;
                sto->totnonBDTDNE_deg[thead]--;
                sto->totBDTDNE_deg[ttail]++;
                sto->totBDTDNE_deg[thead]++;            
              }
            }
          }
        }
      }
    }
    
    if(sto->headmaxl) {
      // head will be newly submaxl after toggle, so add it to the appropriate node list    
      sto->nodesvec[sto->stratheadtype][sto->bdheadtype][sto->attrcounts[sto->stratheadtype][sto->bdheadtype]] = head;
      sto->nodepos[head - 1] = sto->attrcounts[sto->stratheadtype][sto->bdheadtype];
      sto->attrcounts[sto->stratheadtype][sto->bdheadtype]++;

      if(sto->totnonBDTDNE_deg[head] > 0) {
        for(int i = 0; i < sto->nmixtypes; i++) {
          // iterate over all nonBDTDNEs with head as an endpoint; if other endpoint is also
          // submaxl, move this edge from nonBDTDNE to BDTDNE
          if(sto->nonBDTDNE_deg[i][head] > 0) {
            for(int j = sto->nonBDTDNE[i]->nedges; j >= 1; j--) {
              Vertex ttail = sto->nonBDTDNE[i]->tails[j];
              Vertex thead = sto->nonBDTDNE[i]->heads[j];
              if((head == ttail && IN_DEG[thead] + OUT_DEG[thead] < sto->bound) || (head == thead && IN_DEG[ttail] + OUT_DEG[ttail] < sto->bound)) {
                sto->nonBDTDNE[i]->lindex = j;
                sto->nonBDTDNE[i]->ltail = ttail;
                sto->nonBDTDNE[i]->lhead = thead;
                UnsrtELDelete(ttail, thead, sto->nonBDTDNE[i]);
                UnsrtELInsert(ttail, thead, sto->BDTDNE[i]);
                sto->nonBDTDNE_deg[i][ttail]--;
                sto->nonBDTDNE_deg[i][thead]--;
                sto->BDTDNE_deg[i][ttail]++;
                sto->BDTDNE_deg[i][thead]++;            
  
                sto->totnonBDTDNE_deg[ttail]--;
                sto->totnonBDTDNE_deg[thead]--;
                sto->totBDTDNE_deg[ttail]++;
                sto->totBDTDNE_deg[thead]++;
              }
            }
          }
        }
      }
    }
  } else {
    if(sto->tailmaxl) {
      // tail will be newly maxl after toggle, so remove it from the appropriate node list
      sto->nodesvec[sto->strattailtype][sto->bdtailtype][sto->nodepos[tail - 1]] = sto->nodesvec[sto->strattailtype][sto->bdtailtype][sto->attrcounts[sto->strattailtype][sto->bdtailtype] - 1];
      sto->nodepos[sto->nodesvec[sto->strattailtype][sto->bdtailtype][sto->nodepos[tail - 1]] - 1] = sto->nodepos[tail - 1];
      sto->attrcounts[sto->strattailtype][sto->bdtailtype]--;

      if(sto->totBDTDNE_deg[tail] > 0) {
        for(int i = 0; i < sto->nmixtypes; i++) {
          // transfer all BDTDNEs with tail as an endpoint to nonBDTDNE
          if(sto->BDTDNE_deg[i][tail] > 0) {
            for(int j = sto->BDTDNE[i]->nedges; j >= 1; j--) {
              Vertex ttail = sto->BDTDNE[i]->tails[j];
              Vertex thead = sto->BDTDNE[i]->heads[j];
              if(tail == ttail || tail == thead) {
                sto->BDTDNE[i]->lindex = j;
                sto->BDTDNE[i]->ltail = ttail;
                sto->BDTDNE[i]->lhead = thead;
                UnsrtELDelete(ttail, thead, sto->BDTDNE[i]);
                UnsrtELInsert(ttail, thead, sto->nonBDTDNE[i]);
                sto->nonBDTDNE_deg[i][ttail]++;
                sto->nonBDTDNE_deg[i][thead]++;
                sto->BDTDNE_deg[i][ttail]--;
                sto->BDTDNE_deg[i][thead]--;            
  
                sto->totnonBDTDNE_deg[ttail]++;
                sto->totnonBDTDNE_deg[thead]++;
                sto->totBDTDNE_deg[ttail]--;
                sto->totBDTDNE_deg[thead]--;            
              }
            }        
          }
        }          
      }
    }
    
    if(sto->headmaxl) {
      // head will be newly maxl after toggle, so remove it from the appropriate node list
      sto->nodesvec[sto->stratheadtype][sto->bdheadtype][sto->nodepos[head - 1]] = sto->nodesvec[sto->stratheadtype][sto->bdheadtype][sto->attrcounts[sto->stratheadtype][sto->bdheadtype] - 1];
      sto->nodepos[sto->nodesvec[sto->stratheadtype][sto->bdheadtype][sto->nodepos[head - 1]] - 1] = sto->nodepos[head - 1];
      sto->attrcounts[sto->stratheadtype][sto->bdheadtype]--;

      if(sto->totBDTDNE_deg[head] > 0) {
        for(int i = 0; i < sto->nmixtypes; i++) {
          if(sto->BDTDNE_deg[i][head] > 0) {
            for(int j = sto->BDTDNE[i]->nedges; j >= 1; j--) {
              Vertex ttail = sto->BDTDNE[i]->tails[j];
              Vertex thead = sto->BDTDNE[i]->heads[j];
              if(head == ttail || head == thead) {
                sto->BDTDNE[i]->lindex = j;
                sto->BDTDNE[i]->ltail = ttail;
                sto->BDTDNE[i]->lhead = thead;
                UnsrtELDelete(ttail, thead, sto->BDTDNE[i]);
                UnsrtELInsert(ttail, thead, sto->nonBDTDNE[i]);
                sto->nonBDTDNE_deg[i][ttail]++;
                sto->nonBDTDNE_deg[i][thead]++;
                sto->BDTDNE_deg[i][ttail]--;
                sto->BDTDNE_deg[i][thead]--;                        
  
                sto->totnonBDTDNE_deg[ttail]++;
                sto->totnonBDTDNE_deg[thead]++;
                sto->totBDTDNE_deg[ttail]--;
                sto->totBDTDNE_deg[thead]--;                        
              }
            }        
          }
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
  
  Free(sto->proposedprobvec);
  Free(sto->currentprobvec);
  
  Free(sto->nodepos);
    
  for(int i = 0; i < sto->nmixtypes; i++) {
    UnsrtELDestroy(sto->BDTDNE[i]);
    UnsrtELDestroy(sto->nonBDTDNE[i]);
    
    UnsrtELDestroy(sto->discordantEdges[i]);
    UnsrtELDestroy(sto->nonDiscordantEdges[i]);
    
    Free(sto->BDTDNE_deg[i]);
    Free(sto->nonBDTDNE_deg[i]);
  }
  
  Free(sto->BDTDNE);
  Free(sto->nonBDTDNE);
  Free(sto->discordantEdges);
  Free(sto->nonDiscordantEdges);
  Free(sto->BDTDNE_deg);
  Free(sto->nonBDTDNE_deg);
  Free(sto->totBDTDNE_deg);
  Free(sto->totnonBDTDNE_deg);
  // MHp->storage itself should be Freed by MHProposalDestroy
}
