#include "changestat.h"
#include "R_ext/Rdynload.h"
#include "edgetree.h"
#include "edgelist.h"

double my_choose(double n, int r){
static double (*fun)(double,int) = NULL;
if(fun==NULL) fun = (double (*)(double,int)) R_FindSymbol("my_choose", "ergm", NULL);
return fun(n,r);
}
