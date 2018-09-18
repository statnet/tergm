/*  File src/MHproposals_DynMLE.h in package tergm, part of the Statnet suite
 *  of packages for network analysis, http://statnet.org .
 *
 *  This software is distributed under the GPL-3 license.  It is free,
 *  open source, and has the attribution requirements (GPL Section 7) at
 *  http://statnet.org/attribution
 *
 *  Copyright 2008-2018 Statnet Commons
 */
#ifndef MHproposals_DynMLE_H
#define MHproposals_DynMLE_H

#include "edgetree.h"
#include "changestat.h"
#include "model.h"
#include "MHproposal.h"

void MH_FormationMLE(MHProposal *MHp, Network *nwp);
void MH_FormationMLETNT(MHProposal *MHp, Network *nwp);
void MH_DissolutionMLETNT(MHProposal *MHp, Network *nwp);

#endif 

