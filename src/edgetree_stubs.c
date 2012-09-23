#include "edgetree.h"
#include "R_ext/Rdynload.h"
Network NetworkInitialize(Vertex *tails, Vertex *heads, Edge nedges,Vertex nnodes, int directed_flag, Vertex bipartite,int lasttoggle_flag, int time, int *lasttoggle){
static Network (*fun)(Vertex *,Vertex *,Edge,Vertex,int,Vertex,int,int,int *) = NULL;
if(fun==NULL) fun = (Network (*)(Vertex *,Vertex *,Edge,Vertex,int,Vertex,int,int,int *)) R_FindSymbol("NetworkInitialize", "ergm", NULL);
return fun(tails,heads,nedges,nnodes,directed_flag,bipartite,lasttoggle_flag,time,lasttoggle);
}
void NetworkDestroy(Network *nwp){
static void (*fun)(Network *) = NULL;
if(fun==NULL) fun = (void (*)(Network *)) R_FindSymbol("NetworkDestroy", "ergm", NULL);
fun(nwp);
}
Network NetworkInitializeD(double *tails, double *heads, Edge nedges,Vertex nnodes, int directed_flag, Vertex bipartite,int lasttoggle_flag, int time, int *lasttoggle){
static Network (*fun)(double *,double *,Edge,Vertex,int,Vertex,int,int,int *) = NULL;
if(fun==NULL) fun = (Network (*)(double *,double *,Edge,Vertex,int,Vertex,int,int,int *)) R_FindSymbol("NetworkInitializeD", "ergm", NULL);
return fun(tails,heads,nedges,nnodes,directed_flag,bipartite,lasttoggle_flag,time,lasttoggle);
}
Network * NetworkCopy(Network *dest, Network *src){
static Network * (*fun)(Network *,Network *) = NULL;
if(fun==NULL) fun = (Network * (*)(Network *,Network *)) R_FindSymbol("NetworkCopy", "ergm", NULL);
return fun(dest,src);
}
Edge EdgetreeSearch(Vertex a, Vertex b, TreeNode *edges){
static Edge (*fun)(Vertex,Vertex,TreeNode *) = NULL;
if(fun==NULL) fun = (Edge (*)(Vertex,Vertex,TreeNode *)) R_FindSymbol("EdgetreeSearch", "ergm", NULL);
return fun(a,b,edges);
}
Edge EdgetreeSuccessor(TreeNode *edges, Edge x){
static Edge (*fun)(TreeNode *,Edge) = NULL;
if(fun==NULL) fun = (Edge (*)(TreeNode *,Edge)) R_FindSymbol("EdgetreeSuccessor", "ergm", NULL);
return fun(edges,x);
}
Edge EdgetreePredecessor(TreeNode *edges, Edge x){
static Edge (*fun)(TreeNode *,Edge) = NULL;
if(fun==NULL) fun = (Edge (*)(TreeNode *,Edge)) R_FindSymbol("EdgetreePredecessor", "ergm", NULL);
return fun(edges,x);
}
Edge EdgetreeMinimum(TreeNode *edges, Edge x){
static Edge (*fun)(TreeNode *,Edge) = NULL;
if(fun==NULL) fun = (Edge (*)(TreeNode *,Edge)) R_FindSymbol("EdgetreeMinimum", "ergm", NULL);
return fun(edges,x);
}
Edge EdgetreeMaximum(TreeNode *edges, Edge x){
static Edge (*fun)(TreeNode *,Edge) = NULL;
if(fun==NULL) fun = (Edge (*)(TreeNode *,Edge)) R_FindSymbol("EdgetreeMaximum", "ergm", NULL);
return fun(edges,x);
}
int ToggleEdge(Vertex tail, Vertex head, Network *nwp){
static int (*fun)(Vertex,Vertex,Network *) = NULL;
if(fun==NULL) fun = (int (*)(Vertex,Vertex,Network *)) R_FindSymbol("ToggleEdge", "ergm", NULL);
return fun(tail,head,nwp);
}
int ToggleEdgeWithTimestamp(Vertex tail, Vertex head, Network *nwp){
static int (*fun)(Vertex,Vertex,Network *) = NULL;
if(fun==NULL) fun = (int (*)(Vertex,Vertex,Network *)) R_FindSymbol("ToggleEdgeWithTimestamp", "ergm", NULL);
return fun(tail,head,nwp);
}
int AddEdgeToTrees(Vertex tail, Vertex head, Network *nwp){
static int (*fun)(Vertex,Vertex,Network *) = NULL;
if(fun==NULL) fun = (int (*)(Vertex,Vertex,Network *)) R_FindSymbol("AddEdgeToTrees", "ergm", NULL);
return fun(tail,head,nwp);
}
void AddHalfedgeToTree(Vertex a, Vertex b, TreeNode *edges, Edge next_edge){
static void (*fun)(Vertex,Vertex,TreeNode *,Edge) = NULL;
if(fun==NULL) fun = (void (*)(Vertex,Vertex,TreeNode *,Edge)) R_FindSymbol("AddHalfedgeToTree", "ergm", NULL);
fun(a,b,edges,next_edge);
}
void UpdateNextedge(TreeNode *edges, Edge *nextedge, Network *nwp){
static void (*fun)(TreeNode *,Edge *,Network *) = NULL;
if(fun==NULL) fun = (void (*)(TreeNode *,Edge *,Network *)) R_FindSymbol("UpdateNextedge", "ergm", NULL);
fun(edges,nextedge,nwp);
}
int DeleteEdgeFromTrees(Vertex tail, Vertex head, Network *nwp){
static int (*fun)(Vertex,Vertex,Network *) = NULL;
if(fun==NULL) fun = (int (*)(Vertex,Vertex,Network *)) R_FindSymbol("DeleteEdgeFromTrees", "ergm", NULL);
return fun(tail,head,nwp);
}
int DeleteHalfedgeFromTree(Vertex a, Vertex b, TreeNode *edges,Edge *next_edge){
static int (*fun)(Vertex,Vertex,TreeNode *,Edge *) = NULL;
if(fun==NULL) fun = (int (*)(Vertex,Vertex,TreeNode *,Edge *)) R_FindSymbol("DeleteHalfedgeFromTree", "ergm", NULL);
return fun(a,b,edges,next_edge);
}
int ElapsedTime(Vertex tail, Vertex head, Network *nwp){
static int (*fun)(Vertex,Vertex,Network *) = NULL;
if(fun==NULL) fun = (int (*)(Vertex,Vertex,Network *)) R_FindSymbol("ElapsedTime", "ergm", NULL);
return fun(tail,head,nwp);
}
void TouchEdge(Vertex tail, Vertex head, Network *nwp){
static void (*fun)(Vertex,Vertex,Network *) = NULL;
if(fun==NULL) fun = (void (*)(Vertex,Vertex,Network *)) R_FindSymbol("TouchEdge", "ergm", NULL);
fun(tail,head,nwp);
}
int FindithEdge(Vertex *tail, Vertex *head, Edge i, Network *nwp){
static int (*fun)(Vertex *,Vertex *,Edge,Network *) = NULL;
if(fun==NULL) fun = (int (*)(Vertex *,Vertex *,Edge,Network *)) R_FindSymbol("FindithEdge", "ergm", NULL);
return fun(tail,head,i,nwp);
}
int GetRandEdge(Vertex *tail, Vertex *head, Network *nwp){
static int (*fun)(Vertex *,Vertex *,Network *) = NULL;
if(fun==NULL) fun = (int (*)(Vertex *,Vertex *,Network *)) R_FindSymbol("GetRandEdge", "ergm", NULL);
return fun(tail,head,nwp);
}
int FindithNondge(Vertex *tail, Vertex *head, Edge i, Network *nwp){
static int (*fun)(Vertex *,Vertex *,Edge,Network *) = NULL;
if(fun==NULL) fun = (int (*)(Vertex *,Vertex *,Edge,Network *)) R_FindSymbol("FindithNondge", "ergm", NULL);
return fun(tail,head,i,nwp);
}
int GetRandNonedge(Vertex *tail, Vertex *head, Network *nwp){
static int (*fun)(Vertex *,Vertex *,Network *) = NULL;
if(fun==NULL) fun = (int (*)(Vertex *,Vertex *,Network *)) R_FindSymbol("GetRandNonedge", "ergm", NULL);
return fun(tail,head,nwp);
}
void printedge(Edge e, TreeNode *edges){
static void (*fun)(Edge,TreeNode *) = NULL;
if(fun==NULL) fun = (void (*)(Edge,TreeNode *)) R_FindSymbol("printedge", "ergm", NULL);
fun(e,edges);
}
void InOrderTreeWalk(TreeNode *edges, Edge x){
static void (*fun)(TreeNode *,Edge) = NULL;
if(fun==NULL) fun = (void (*)(TreeNode *,Edge)) R_FindSymbol("InOrderTreeWalk", "ergm", NULL);
fun(edges,x);
}
void NetworkEdgeList(Network *nwp){
static void (*fun)(Network *) = NULL;
if(fun==NULL) fun = (void (*)(Network *)) R_FindSymbol("NetworkEdgeList", "ergm", NULL);
fun(nwp);
}
void ShuffleEdges(Vertex *tails, Vertex *heads, Edge nedges){
static void (*fun)(Vertex *,Vertex *,Edge) = NULL;
if(fun==NULL) fun = (void (*)(Vertex *,Vertex *,Edge)) R_FindSymbol("ShuffleEdges", "ergm", NULL);
fun(tails,heads,nedges);
}
Edge DesignMissing(Vertex a, Vertex b, Network *mnwp){
static Edge (*fun)(Vertex,Vertex,Network *) = NULL;
if(fun==NULL) fun = (Edge (*)(Vertex,Vertex,Network *)) R_FindSymbol("DesignMissing", "ergm", NULL);
return fun(a,b,mnwp);
}
Edge EdgeTree2EdgeList(Vertex *tails, Vertex *heads, Network *nwp, Edge nmax){
static Edge (*fun)(Vertex *,Vertex *,Network *,Edge) = NULL;
if(fun==NULL) fun = (Edge (*)(Vertex *,Vertex *,Network *,Edge)) R_FindSymbol("EdgeTree2EdgeList", "ergm", NULL);
return fun(tails,heads,nwp,nmax);
}
