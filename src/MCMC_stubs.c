#include "MCMC.h"
#include "R_ext/Rdynload.h"
#include "edgetree.h"
#include "changestat.h"
#include "MHproposal.h"
#include "model.h"

void MCMC_wrapper(int *dnumnets, int *dnedges,int *tails, int *heads,int *dn, int *dflag, int *bipartite,int *nterms, char **funnames,char **sonames,char **MHproposaltype, char **MHproposalpackage,double *inputs, double *theta0, int *samplesize,double *sample, int *burnin, int *interval,int *newnetworktails,int *newnetworkheads,int *fVerbose,int *attribs, int *maxout, int *maxin, int *minout,int *minin, int *condAllDegExact, int *attriblength,int *maxedges,int *status){
static void (*fun)(int *,int *,int *,int *,int *,int *,int *,int *,char **,char **,char **,char **,double *,double *,int *,double *,int *,int *,int *,int *,int *,int *,int *,int *,int *,int *,int *,int *,int *,int *) = NULL;
if(fun==NULL) fun = (void (*)(int *,int *,int *,int *,int *,int *,int *,int *,char **,char **,char **,char **,double *,double *,int *,double *,int *,int *,int *,int *,int *,int *,int *,int *,int *,int *,int *,int *,int *,int *)) R_FindSymbol("MCMC_wrapper", "ergm", NULL);
fun(dnumnets,dnedges,tails,heads,dn,dflag,bipartite,nterms,funnames,sonames,MHproposaltype,MHproposalpackage,inputs,theta0,samplesize,sample,burnin,interval,newnetworktails,newnetworkheads,fVerbose,attribs,maxout,maxin,minout,minin,condAllDegExact,attriblength,maxedges,status);
}
MCMCStatus MCMCSample(MHproposal *MHp,double *theta, double *networkstatistics,int samplesize, int burnin,int interval, int fVerbose, int nmax,Network *nwp, Model *m){
static MCMCStatus (*fun)(MHproposal *,double *,double *,int,int,int,int,int,Network *,Model *) = NULL;
if(fun==NULL) fun = (MCMCStatus (*)(MHproposal *,double *,double *,int,int,int,int,int,Network *,Model *)) R_FindSymbol("MCMCSample", "ergm", NULL);
return fun(MHp,theta,networkstatistics,samplesize,burnin,interval,fVerbose,nmax,nwp,m);
}
MCMCStatus MetropolisHastings(MHproposal *MHp,double *theta, double *statistics,int nsteps, int *staken,int fVerbose,Network *nwp, Model *m){
static MCMCStatus (*fun)(MHproposal *,double *,double *,int,int *,int,Network *,Model *) = NULL;
if(fun==NULL) fun = (MCMCStatus (*)(MHproposal *,double *,double *,int,int *,int,Network *,Model *)) R_FindSymbol("MetropolisHastings", "ergm", NULL);
return fun(MHp,theta,statistics,nsteps,staken,fVerbose,nwp,m);
}
void MCMCPhase12(int *tails, int *heads, int *dnedges,int *dn, int *dflag, int *bipartite,int *nterms, char **funnames,char **sonames,char **MHproposaltype, char **MHproposalpackage,double *inputs,double *theta0, int *samplesize,double *gain, double *meanstats, int *phase1, int *nsub,double *sample, int *burnin, int *interval,int *newnetworktails,int *newnetworkheads,int *fVerbose,int *attribs, int *maxout, int *maxin, int *minout,int *minin, int *condAllDegExact, int *attriblength,int *maxedges,int *mtails, int *mheads, int *mdnedges){
static void (*fun)(int *,int *,int *,int *,int *,int *,int *,char **,char **,char **,char **,double *,double *,int *,double *,double *,int *,int *,double *,int *,int *,int *,int *,int *,int *,int *,int *,int *,int *,int *,int *,int *,int *,int *,int *) = NULL;
if(fun==NULL) fun = (void (*)(int *,int *,int *,int *,int *,int *,int *,char **,char **,char **,char **,double *,double *,int *,double *,double *,int *,int *,double *,int *,int *,int *,int *,int *,int *,int *,int *,int *,int *,int *,int *,int *,int *,int *,int *)) R_FindSymbol("MCMCPhase12", "ergm", NULL);
fun(tails,heads,dnedges,dn,dflag,bipartite,nterms,funnames,sonames,MHproposaltype,MHproposalpackage,inputs,theta0,samplesize,gain,meanstats,phase1,nsub,sample,burnin,interval,newnetworktails,newnetworkheads,fVerbose,attribs,maxout,maxin,minout,minin,condAllDegExact,attriblength,maxedges,mtails,mheads,mdnedges);
}
void MCMCSamplePhase12(MHproposal *MH,double *theta, double gain, double *meanstats,int nphase1, int nsubphases, double *networkstatistics,int samplesize, int burnin,int interval, int fVerbose,Network *nwp, Model *m){
static void (*fun)(MHproposal *,double *,double,double *,int,int,double *,int,int,int,int,Network *,Model *) = NULL;
if(fun==NULL) fun = (void (*)(MHproposal *,double *,double,double *,int,int,double *,int,int,int,int,Network *,Model *)) R_FindSymbol("MCMCSamplePhase12", "ergm", NULL);
fun(MH,theta,gain,meanstats,nphase1,nsubphases,networkstatistics,samplesize,burnin,interval,fVerbose,nwp,m);
}
