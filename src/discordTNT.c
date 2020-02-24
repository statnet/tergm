#include "ergm_MHproposal.h"
#include "ergm_edgelist.h"
#include "ergm_changestat.h"
#include "ergm_rlebdm.h"
#include "ergm_MHstorage.h"
#include "ergm_unsorted_edgelist.h"
#include "changestats_lasttoggle.h"
#include "tergm_model.h"

MH_I_FN(Mi_discordTNT) {
  MHp->ntoggles = 1;
  
  ALLOC_STORAGE(3, void *, sto);
  sto[0] = UnsrtELInitialize(0, NULL, NULL, FALSE);
  sto[1] = UnsrtELInitialize(0, NULL, NULL, FALSE);
  sto[2] = Calloc(1, int);
  
  // we ignore discord for this initialization (assuming a TICK will precede any proposals)
}

MH_X_FN(Mx_discordTNT) {
  GET_STORAGE(void *, sto);
  
  if(type == TICK) {
    // "clear" the discordant dyads
    ((UnsrtEL *)sto[0])->nedges = 0;
    ((UnsrtEL *)sto[1])->nedges = 0;
  }
}

MH_P_FN(MH_discordTNT) {
  GET_AUX_STORAGE(StoreTimeAndLasttoggle, dur_inf);
  
  GET_STORAGE(void *, sto);
  UnsrtEL *discordantEdges = sto[0];
  UnsrtEL *discordantNonEdges = sto[1];
  
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
        UnsrtELGetRand(Mtail, Mhead, discordantEdges);
      } else {
        UnsrtELGetRand(Mtail, Mhead, discordantNonEdges);
      }
    }
  } else {
    // propose from discord
    if(unif_rand() < discordantEdges->nedges/((double) nddyads)) {
      UnsrtELGetRand(Mtail, Mhead, discordantEdges);
      in_network = TRUE;
    } else {
      UnsrtELGetRand(Mtail, Mhead, discordantNonEdges);
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

  ((int *)sto[2])[0] = in_discord;
}

MH_U_FN(Mu_discordTNT) {
  // add or remove the dyad from the appropriate discordance edgelist
  
  GET_STORAGE(void *, sto);
  
  if(((int *)sto[2])[0]) {
    // currently in discord; take it out
    if(edgeflag) {        
      UnsrtELDelete(tail, head, sto[0]);
    } else {
      UnsrtELDelete(tail, head, sto[1]);  
    }
  } else {
    // not currently in discord; add it in
    if(edgeflag) {
      UnsrtELInsert(tail, head, sto[1]);
    } else {
      UnsrtELInsert(tail, head, sto[0]);        
    }
  }
}

MH_F_FN(Mf_discordTNT) {
  GET_STORAGE(void *, sto);
  
  Free(sto[2]);
  UnsrtELDestroy(sto[1]);
  UnsrtELDestroy(sto[0]);
}



/********************
   discordStratTNT
********************/

MH_I_FN(Mi_discordStratTNT) {
  // process the inputs and initialize all the edgelists in storage; set MHp->ntoggles to 1
  MHp->ntoggles = 1;
  
  int nmixtypes = MHp->inputs[0];
    
  int nattrcodes = MHp->inputs[1 + 3*nmixtypes];
  
  double *vattr = MHp->inputs + 1 + 3*nmixtypes + 1 + nattrcodes + N_NODES;
  
  ALLOC_STORAGE(6, void *, sto);
  
  UnsrtEL **nonDiscordantELs = (UnsrtEL **)Calloc(nmixtypes, UnsrtEL *);
  UnsrtEL **discordantELs = (UnsrtEL **)Calloc(nmixtypes, UnsrtEL *);
  UnsrtEL **discordantNonELs = (UnsrtEL **)Calloc(nmixtypes, UnsrtEL *);
  
  for(int i = 0; i < nmixtypes; i++) {
    nonDiscordantELs[i] = UnsrtELInitialize(0, NULL, NULL, FALSE);
    discordantELs[i] = UnsrtELInitialize(0, NULL, NULL, FALSE);
    discordantNonELs[i] = UnsrtELInitialize(0, NULL, NULL, FALSE);
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
        UnsrtELInsert(tail, head, nonDiscordantELs[index]);
      }
    }
  }
  Free(indmat);
  
  // assign edgelists to storage
  sto[0] = nonDiscordantELs;
  sto[1] = discordantELs;
  sto[2] = discordantNonELs;
 
  // will track mixing type of current proposal
  sto[3] = Calloc(1, int);

  // will track if current dyad is discordant
  sto[4] = Calloc(1, int);
  
  double *nodecountsbycode = MHp->inputs + 1 + 3*nmixtypes + 1;
  
  double **nodesbycode = (double **)Calloc(nattrcodes, double *);
  nodesbycode[0] = MHp->inputs + 1 + 3*nmixtypes + 1 + nattrcodes;
  for(int i = 1; i < nattrcodes; i++) {
    nodesbycode[i] = nodesbycode[i - 1] + (int)nodecountsbycode[i - 1];
  }
  
  sto[5] = nodesbycode;
} 

MH_X_FN(Mx_discordStratTNT) {
  GET_STORAGE(void *, sto);
    
  if(type == TICK) {
    // transfer discordant edges to nondiscordant edges
    // clear discordant edges and discordant nonedges
    int nmixtypes = MHp->inputs[0];
    
    UnsrtEL **nonDiscordantELs = sto[0];
    UnsrtEL **discordantELs = sto[1];
    UnsrtEL **discordantNonELs = sto[2];
    
    for(int i = 0; i < nmixtypes; i++) {
      Vertex *tails = discordantELs[i]->tails;
      Vertex *heads = discordantELs[i]->heads;
      int nedges = discordantELs[i]->nedges;
      
      for(int j = 0; j < nedges; j++) {
        UnsrtELInsert(tails[j], heads[j], nonDiscordantELs[i]);
      }
      
      discordantELs[i]->nedges = 0;
      discordantNonELs[i]->nedges = 0;      
    }
  }
}

MH_P_FN(MH_discordStratTNT) {
  GET_AUX_STORAGE(StoreTimeAndLasttoggle, dur_inf);
  GET_STORAGE(void *, sto);
  
  int nmixtypes = MHp->inputs[0];
  int nattrcodes = MHp->inputs[1 + 3*nmixtypes];
  
  double *pmat = MHp->inputs + 1 + 2*nmixtypes;
  double *nodecountsbycode = MHp->inputs + 1 + 3*nmixtypes + 1;  
  
  UnsrtEL **nonDiscordantELs = sto[0];
  UnsrtEL **discordantELs = sto[1];
  UnsrtEL **discordantNonELs = sto[2];

  double **nodesbycode = sto[5];
    
  double ur = unif_rand();
  
  // find the first mixing type i with (cumulative) probability larger than ur
  int i = 0;
  while(ur > pmat[i]) {
    i++;
  }
  
  // record the mixing type of the toggle, in case it's needed in the U function later
  ((int *)sto[3])[0] = i;    
  
  int tailtype = MHp->inputs[1 + i];
  int headtype = MHp->inputs[1 + nmixtypes + i];
    
  // number of edges of this mixing type
  int nedgestype = nonDiscordantELs[i]->nedges + discordantELs[i]->nedges;

  // number of dyads of this mixing type
  int ndyadstype = MHp->inputs[1 + 3*nmixtypes + 1 + nattrcodes + N_NODES + N_NODES + i];
  
  // number of discordant dyads of this mixing type
  int nddyadstype = discordantNonELs[i]->nedges + discordantELs[i]->nedges;
  
  // flags
  int in_discord;
  int in_network;
  
  if(nddyadstype == 0 || unif_rand() < 0.5) {
    // propose from network
    if(nedgestype == 0 || unif_rand() < 0.5) {
      // propose toggling a random dyad of the specified mixing type
      
      int tailindex = nodecountsbycode[tailtype]*unif_rand();
      int headindex;
      if(tailtype == headtype) {
        // need to avoid sampling a loop
        headindex = (nodecountsbycode[headtype] - 1)*unif_rand();
        if(headindex == tailindex) {
          headindex = nodecountsbycode[headtype] - 1;
        }
      } else {
        // any old head will do
        headindex = nodecountsbycode[headtype]*unif_rand();
      }
            
      Vertex tail = nodesbycode[tailtype][tailindex];
      Vertex head = nodesbycode[headtype][headindex];
      
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
          UnsrtELGetRand(Mtail, Mhead, discordantELs[i]);
        } else {
          UnsrtELGetRand(Mtail, Mhead, nonDiscordantELs[i]);
        }
      } else if(in_discord) {
        UnsrtELGetRand(Mtail, Mhead, discordantNonELs[i]);
      }
    } else {
      // propose toggling off an edge of the specified mixing type
      if(unif_rand() < nonDiscordantELs[i]->nedges/((double) nedgestype)) {
        UnsrtELGetRand(Mtail, Mhead, nonDiscordantELs[i]);
        in_discord = FALSE;
      } else {
        UnsrtELGetRand(Mtail, Mhead, discordantELs[i]);
        in_discord = TRUE;
      }
      in_network = TRUE;
    }
  } else {
    // propose from discord
    if(unif_rand() < discordantELs[i]->nedges/((double) nddyadstype)) {
      UnsrtELGetRand(Mtail, Mhead, discordantELs[i]);
      in_network = FALSE;
    } else {
      UnsrtELGetRand(Mtail, Mhead, discordantNonELs[i]);
      in_network = TRUE;
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

  ((int *)sto[4])[0] = in_discord;
}

MH_U_FN(Mu_discordStratTNT) {
  // add or remove edge from appropriate edgelist
  GET_STORAGE(void *, sto);
  
  UnsrtEL **nonDiscordantELs = sto[0];
  UnsrtEL **discordantELs = sto[1];
  UnsrtEL **discordantNonELs = sto[2];
  
  int i = ((int *)sto[3])[0];
  int in_discord = ((int *)sto[4])[0];
  
  if(edgeflag) {
    // we are removing an existing edge
    if(in_discord) {
      UnsrtELDelete(tail, head, discordantELs[i]);
    } else {
      UnsrtELDelete(tail, head, nonDiscordantELs[i]);
      UnsrtELInsert(tail, head, discordantNonELs[i]);
    }
  } else {
    // we are adding a new edge
    if(in_discord) {
      UnsrtELInsert(tail, head, nonDiscordantELs[i]);
      UnsrtELDelete(tail, head, discordantNonELs[i]);
    } else {
      UnsrtELInsert(tail, head, discordantELs[i]);
    }
  }
}

MH_F_FN(Mf_discordStratTNT) {
  // Free all the things
  int nmixtypes = MHp->inputs[0];

  GET_STORAGE(void *, sto);
  
  UnsrtEL **nonDiscordantELs = sto[0];
  UnsrtEL **discordantELs = sto[1];
  UnsrtEL **discordantNonELs = sto[2];  
  
  Free(sto[5]);
  Free(sto[4]);
  Free(sto[3]);
  
  for(int i = 0; i < nmixtypes; i++) {
    UnsrtELDestroy(nonDiscordantELs[i]);
    UnsrtELDestroy(discordantELs[i]);
    UnsrtELDestroy(discordantNonELs[i]);    
  }

  Free(nonDiscordantELs);
  Free(discordantELs);
  Free(discordantNonELs);

  // MHp->storage itself should be Freed by MHProposalDestroy
}


