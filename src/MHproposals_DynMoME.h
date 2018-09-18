/*  File src/MHproposals_DynMoME.h in package tergm, part of the Statnet suite
 *  of packages for network analysis, http://statnet.org .
 *
 *  This software is distributed under the GPL-3 license.  It is free,
 *  open source, and has the attribution requirements (GPL Section 7) at
 *  http://statnet.org/attribution
 *
 *  Copyright 2008-2018 Statnet Commons
 */
#ifndef MHproposals_DynMoME_H
#define MHproposals_DynMoME_H

#include "edgetree.h"
#include "changestat.h"
#include "model.h"
#include "MHproposal.h"

void MH_Formation(MHProposal *MHp, Network *nwp);
void MH_FormationTNT(MHProposal *MHp, Network *nwp);
void MH_Dissolution(MHProposal *MHp, Network *nwp);
void MH_DissolutionTNT(MHProposal *MHp, Network *nwp);

#endif 

