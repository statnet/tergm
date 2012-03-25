#include "MHproposal.h"
#include "R_ext/Rdynload.h"
#include "edgetree.h"
#include "R_ext/Rdynload.h"

DegreeBound* DegreeBoundInitialize(int *attribs, int *maxout, int *maxin,int *minout, int *minin, int condAllDegExact,int attriblength, Network *nwp){
static DegreeBound* (*fun)(int *,int *,int *,int *,int *,int,int,Network *) = NULL;
if(fun==NULL) fun = (DegreeBound* (*)(int *,int *,int *,int *,int *,int,int,Network *)) R_FindSymbol("DegreeBoundInitialize", "ergm", NULL);
return fun(attribs,maxout,maxin,minout,minin,condAllDegExact,attriblength,nwp);
}
void DegreeBoundDestroy(DegreeBound *bd){
static void (*fun)(DegreeBound *) = NULL;
if(fun==NULL) fun = (void (*)(DegreeBound *)) R_FindSymbol("DegreeBoundDestroy", "ergm", NULL);
fun(bd);
}
void MH_init(MHproposal *MHp,char *MHproposaltype, char *MHproposalpackage,double *inputs,int fVerbose,Network *nwp,int *attribs, int *maxout, int *maxin,int *minout, int *minin, int condAllDegExact,int attriblength){
static void (*fun)(MHproposal *,char *,char *,double *,int,Network *,int *,int *,int *,int *,int *,int,int) = NULL;
if(fun==NULL) fun = (void (*)(MHproposal *,char *,char *,double *,int,Network *,int *,int *,int *,int *,int *,int,int)) R_FindSymbol("MH_init", "ergm", NULL);
fun(MHp,MHproposaltype,MHproposalpackage,inputs,fVerbose,nwp,attribs,maxout,maxin,minout,minin,condAllDegExact,attriblength);
}
void MH_free(MHproposal *MHp){
static void (*fun)(MHproposal *) = NULL;
if(fun==NULL) fun = (void (*)(MHproposal *)) R_FindSymbol("MH_free", "ergm", NULL);
fun(MHp);
}
int CheckTogglesValid(MHproposal *MHp, Network *nwp){
static int (*fun)(MHproposal *,Network *) = NULL;
if(fun==NULL) fun = (int (*)(MHproposal *,Network *)) R_FindSymbol("CheckTogglesValid", "ergm", NULL);
return fun(MHp,nwp);
}
int CheckConstrainedTogglesValid(MHproposal *MHp, Network *nwp){
static int (*fun)(MHproposal *,Network *) = NULL;
if(fun==NULL) fun = (int (*)(MHproposal *,Network *)) R_FindSymbol("CheckConstrainedTogglesValid", "ergm", NULL);
return fun(MHp,nwp);
}
