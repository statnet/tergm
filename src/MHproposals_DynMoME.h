#ifndef MHproposals_DynMoME_H
#define MHproposals_DynMoME_H

#include "edgetree.h"
#include "changestat.h"
#include "model.h"
#include "MHproposal.h"

void MH_Formation(MHproposal *MHp, Network *nwp);
void MH_FormationTNT(MHproposal *MHp, Network *nwp);
void MH_Dissolution(MHproposal *MHp, Network *nwp);
void MH_DissolutionTNT(MHproposal *MHp, Network *nwp);

#endif 

