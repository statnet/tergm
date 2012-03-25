#include "edgelist.h"
#include "R_ext/Rdynload.h"
unsigned int dEdgeListSearch(Vertex tail, Vertex head, double *el){
static unsigned int (*fun)(Vertex,Vertex,double *) = NULL;
if(fun==NULL) fun = (unsigned int (*)(Vertex,Vertex,double *)) R_FindSymbol("dEdgeListSearch", "ergm", NULL);
return fun(tail,head,el);
}
unsigned int iEdgeListSearch(Vertex tail, Vertex head, int *el){
static unsigned int (*fun)(Vertex,Vertex,int *) = NULL;
if(fun==NULL) fun = (unsigned int (*)(Vertex,Vertex,int *)) R_FindSymbol("iEdgeListSearch", "ergm", NULL);
return fun(tail,head,el);
}
