#ifndef MHproposals_DynMLE_H
#define MHproposals_DynMLE_H

#include "edgetree.h"
#include "changestat.h"
#include "model.h"
#include "MHproposal.h"

void MH_FormationMLEblockdiag(MHproposal *MHp, Network *nwp);
void MH_FormationMLEblockdiagTNT(MHproposal *MHp, Network *nwp);

#endif 

