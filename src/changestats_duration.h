/*  File src/changestats_duration.h in package tergm, part of the Statnet suite
 *  of packages for network analysis, https://statnet.org .
 *
 *  This software is distributed under the GPL-3 license.  It is free,
 *  open source, and has the attribution requirements (GPL Section 7) at
 *  https://statnet.org/attribution
 *
 *  Copyright 2008-2020 Statnet Commons
 */
#ifndef CHANGESTATS_DURATION_H
#define CHANGESTATS_DURATION_H

#include "changestats_lasttoggle.h"
#include "ergm_storage.h"
#include "tergm_model.h"
#include "ergm_Rutil.h"

X_CHANGESTAT_FN(x_edges_ageinterval_mon);
C_CHANGESTAT_FN(c_edges_ageinterval_mon);
S_CHANGESTAT_FN(s_edges_ageinterval_mon);

X_CHANGESTAT_FN(x_edge_ages_mon);
C_CHANGESTAT_FN(c_edge_ages_mon);
S_CHANGESTAT_FN(s_edge_ages_mon);

X_CHANGESTAT_FN(x_edgecov_ages_mon);
C_CHANGESTAT_FN(c_edgecov_ages_mon);
S_CHANGESTAT_FN(s_edgecov_ages_mon);

I_CHANGESTAT_FN(i_mean_age_mon);
X_CHANGESTAT_FN(x_mean_age_mon);
C_CHANGESTAT_FN(c_mean_age_mon);
U_CHANGESTAT_FN(u_mean_age_mon);
S_CHANGESTAT_FN(s_mean_age_mon);

I_CHANGESTAT_FN(i_nodefactor_mean_age);
X_CHANGESTAT_FN(x_nodefactor_mean_age);
C_CHANGESTAT_FN(c_nodefactor_mean_age);
U_CHANGESTAT_FN(u_nodefactor_mean_age);
F_CHANGESTAT_FN(f_nodefactor_mean_age);
S_CHANGESTAT_FN(s_nodefactor_mean_age);

I_CHANGESTAT_FN(i_nodemix_mean_age);
X_CHANGESTAT_FN(x_nodemix_mean_age);
C_CHANGESTAT_FN(c_nodemix_mean_age);
U_CHANGESTAT_FN(u_nodemix_mean_age);
F_CHANGESTAT_FN(f_nodemix_mean_age);
S_CHANGESTAT_FN(s_nodemix_mean_age);

I_CHANGESTAT_FN(i_edgecov_mean_age_mon);
X_CHANGESTAT_FN(x_edgecov_mean_age_mon);
C_CHANGESTAT_FN(c_edgecov_mean_age_mon);
U_CHANGESTAT_FN(u_edgecov_mean_age_mon);
S_CHANGESTAT_FN(s_edgecov_mean_age_mon);

I_CHANGESTAT_FN(i_degree_mean_age_mon);
X_CHANGESTAT_FN(x_degree_mean_age_mon);
C_CHANGESTAT_FN(c_degree_mean_age_mon);
U_CHANGESTAT_FN(u_degree_mean_age_mon);
S_CHANGESTAT_FN(s_degree_mean_age_mon);
F_CHANGESTAT_FN(f_degree_mean_age_mon);

I_CHANGESTAT_FN(i_degree_by_attr_mean_age_mon);
X_CHANGESTAT_FN(x_degree_by_attr_mean_age_mon);
C_CHANGESTAT_FN(c_degree_by_attr_mean_age_mon);
U_CHANGESTAT_FN(u_degree_by_attr_mean_age_mon);
S_CHANGESTAT_FN(s_degree_by_attr_mean_age_mon);
F_CHANGESTAT_FN(f_degree_by_attr_mean_age_mon);

I_CHANGESTAT_FN(i_degrange_mean_age_mon);
X_CHANGESTAT_FN(x_degrange_mean_age_mon);
C_CHANGESTAT_FN(c_degrange_mean_age_mon);
U_CHANGESTAT_FN(u_degrange_mean_age_mon);
S_CHANGESTAT_FN(s_degrange_mean_age_mon);
F_CHANGESTAT_FN(f_degrange_mean_age_mon);

I_CHANGESTAT_FN(i_degrange_by_attr_mean_age_mon);
X_CHANGESTAT_FN(x_degrange_by_attr_mean_age_mon);
C_CHANGESTAT_FN(c_degrange_by_attr_mean_age_mon);
U_CHANGESTAT_FN(u_degrange_by_attr_mean_age_mon);
S_CHANGESTAT_FN(s_degrange_by_attr_mean_age_mon);
F_CHANGESTAT_FN(f_degrange_by_attr_mean_age_mon);

#endif








