/*  File src/changestats_duration.h in package tergm, part of the Statnet suite
 *  of packages for network analysis, https://statnet.org .
 *
 *  This software is distributed under the GPL-3 license.  It is free,
 *  open source, and has the attribution requirements (GPL Section 7) at
 *  https://statnet.org/attribution
 *
 *  Copyright 2008-2019 Statnet Commons
 */
#ifndef CHANGESTATS_DURATION_H
#define CHANGESTATS_DURATION_H

#include "changestats_lasttoggle.h"
#include "ergm_storage.h"
#include "tergm_model.h"

X_CHANGESTAT_FN(x_edges_ageinterval_mon);
D_CHANGESTAT_FN(d_edges_ageinterval_mon);
S_CHANGESTAT_FN(s_edges_ageinterval_mon);

X_CHANGESTAT_FN(x_edge_ages_mon);
D_CHANGESTAT_FN(d_edge_ages_mon);
S_CHANGESTAT_FN(s_edge_ages_mon);

X_CHANGESTAT_FN(x_edgecov_ages_mon);
D_CHANGESTAT_FN(d_edgecov_ages_mon);
S_CHANGESTAT_FN(s_edgecov_ages_mon);

X_CHANGESTAT_FN(x_mean_age_mon);
D_CHANGESTAT_FN(d_mean_age_mon);
S_CHANGESTAT_FN(s_mean_age_mon);

X_CHANGESTAT_FN(x_edgecov_mean_age_mon);
D_CHANGESTAT_FN(d_edgecov_mean_age_mon);
S_CHANGESTAT_FN(s_edgecov_mean_age_mon);

X_CHANGESTAT_FN(x_degree_mean_age_mon);
D_CHANGESTAT_FN(d_degree_mean_age_mon);
S_CHANGESTAT_FN(s_degree_mean_age_mon);

X_CHANGESTAT_FN(x_degree_by_attr_mean_age_mon);
D_CHANGESTAT_FN(d_degree_by_attr_mean_age_mon);
S_CHANGESTAT_FN(s_degree_by_attr_mean_age_mon);

X_CHANGESTAT_FN(x_degrange_mean_age_mon);
D_CHANGESTAT_FN(d_degrange_mean_age_mon);
S_CHANGESTAT_FN(s_degrange_mean_age_mon);

X_CHANGESTAT_FN(x_degrange_by_attr_mean_age_mon);
D_CHANGESTAT_FN(d_degrange_by_attr_mean_age_mon);
S_CHANGESTAT_FN(s_degrange_by_attr_mean_age_mon);

#endif








