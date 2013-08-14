#ifndef MHproposals_DynMLE_H
#define MHproposals_DynMLE_H

#include "edgetree.h"
#include "changestat.h"
#include "model.h"
#include "MHproposal.h"

void MH_FormationMLE(MHproposal *MHp, Network *nwp);
void MH_FormationMLETNT(MHproposal *MHp, Network *nwp);
void MH_DissolutionMLETNT(MHproposal *MHp, Network *nwp);

#endif 

