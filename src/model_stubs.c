#include "model.h"
#include "R_ext/Rdynload.h"
Model* ModelInitialize(char *fnames, char *sonames, double **inputs,int n_terms){
static Model* (*fun)(char *,char *,double **,int) = NULL;
if(fun==NULL) fun = (Model* (*)(char *,char *,double **,int)) R_FindSymbol("ModelInitialize", "ergm", NULL);
return fun(fnames,sonames,inputs,n_terms);
}
void ModelDestroy(Model *m){
static void (*fun)(Model *) = NULL;
if(fun==NULL) fun = (void (*)(Model *)) R_FindSymbol("ModelDestroy", "ergm", NULL);
fun(m);
}
int GetIndexForAttrValue(int value){
static int (*fun)(int) = NULL;
if(fun==NULL) fun = (int (*)(int)) R_FindSymbol("GetIndexForAttrValue", "ergm", NULL);
return fun(value);
}
void ChangeStats(unsigned int ntoggles, Vertex *toggletail, Vertex *togglehead, Network *nwp, Model *m){
static void (*fun)(unsigned int,Vertex *,Vertex *,Network *,Model *) = NULL;
if(fun==NULL) fun = (void (*)(unsigned int,Vertex *,Vertex *,Network *,Model *)) R_FindSymbol("ChangeStats", "ergm", NULL);
fun(ntoggles,toggletail,togglehead,nwp,m);
}
