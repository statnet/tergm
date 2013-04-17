#ifndef godfather_H 
#define godfather_H

#include <R.h>
#include "model.h"
#include "MCMCDyn.h"

/* Function prototypes */
void godfather_wrapper(int *tails, int *heads, int *time, int *lasttoggle, int *n_edges,
		       int *n_nodes, int *directed_flag, int *bip, 
		       int *nterms, char **funnames, char **sonames, double *inputs,
		       int *total_toggles, int *toggletimes, 
		       int *toggletails, int *toggleheads,
		       int *start_time, int *end_time,
		       double *changestats, 
		       int *maxedges,
		       int *newnetworktails, 
		       int *newnetworkheads, 
		       int *fVerbose, 
		       int *status);
      
#endif

