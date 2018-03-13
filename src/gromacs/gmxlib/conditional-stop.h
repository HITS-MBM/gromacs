/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright (c) 1991-2000, University of Groningen, The Netherlands.
 * Copyright (c) 2001-2004, The GROMACS development team.
 * Copyright (c) 2017,2018, by the GROMACS development team, led by
 * Mark Abraham, David van der Spoel, Berk Hess, and Erik Lindahl,
 * and including many others, as listed in the AUTHORS file in the
 * top-level source directory and at http://www.gromacs.org.
 *
 * GROMACS is free software; you can redistribute it and/or
 * modify it under the terms of the GNU Lesser General Public License
 * as published by the Free Software Foundation; either version 2.1
 * of the License, or (at your option) any later version.
 *
 * GROMACS is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 * Lesser General Public License for more details.
 *
 * You should have received a copy of the GNU Lesser General Public
 * License along with GROMACS; if not, see
 * http://www.gnu.org/licenses, or write to the Free Software Foundation,
 * Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301  USA.
 *
 * If you want to redistribute modifications to GROMACS, please
 * consider that scientific software is very special. Version
 * control is crucial - bugs must be traceable. We will be happy to
 * consider code for inclusion in the official distribution, but
 * derived work must not be called official GROMACS. Details are found
 * in the README & COPYING files - if they are missing, get the
 * official version at http://www.gromacs.org.
 *
 * To help us fund GROMACS development, we humbly ask that you cite
 * the research papers on the package. Check out http://www.gromacs.org.
 */
#ifndef GMX_CONDITIONAL_STOP_H
#define GMX_CONDITIONAL_STOP_H

#include <cstdio>

#include "gromacs/math/vectypes.h"
#include "gromacs/utility/basedefinitions.h"

struct t_commrec;
struct t_state;
struct t_condstop;
struct t_stopcond;

/*! \brief Compute center of coordinates.
 *
 * This function computes center of given coordinates.
 *
 * \param[in]     x              The coordinates.
 * \param[in]     n              The number of coordinates.
 * \param[out]    coc            The center of coordinates.
 */
void compute_coc(const rvec* x, const int n, rvec &coc,
                 const char* padding = "      ");

/*! \brief Compute the distance and check whether the condition is satisfied.
 *
 * From the condition and two sets of atom coordinates, compute the distance of
 * the correct type and check whether the specified condition is satisfied.
 *
 * \param[in]     cond             The conditional stop struct.
 * \param[in]     x1               The coordinates of the atoms in the first group.
 * \param[in]     nat1             The number of coordinates of the atoms in the first group.
 * \param[in]     x2               The coordinates of the atoms in the second group.
 * \param[in]     nat2             The number of coordinates of the atoms in the second group.
 * \param[out]    distance_squared The computed distance.
 * \returns whether the condition was satisfied.
 */
gmx_bool compute_condition_distance(const t_stopcond *cond,
                                    const rvec *x1, const int nat1,
                                    const rvec *x2, const int nat2,
                                    real &distance_squared,
                                    const char* padding = "      ");

/**/
void
store_condgrp_stateidx(t_condstop* condstop);

/**/
void
bcast_condgrp_statex(t_commrec* cr, t_condstop* condstop, t_state* state_global);

/*! \brief Store the current atom coordinates as origin.
 *
 * This function stores the atom coordinates of the atoms contained in the
 * conditional groups in the stop conditions that have bOrigin enabled. This
 * should normally be ran once at the beginning of the simulation.
 *
 * \param[in]     fplog          The general output file, normally md.log.
 * \param[in]     cr             The struct for communication info.
 * \param[in,out] condstop       The conditional stop struct.
 * \param[in]     state_global   The global state (only atom coordinates are used).
 * \param[in]     state_local    The local state (only atom coordinates are used).
 */
void store_coord_as_origin(FILE                 *fplog,
                           t_commrec            *cr,
                           t_condstop           *condstop,
                           const t_state        *state_global,
                           const t_state        *state_local);

/*! \brief Check if a stop condition is satisfied.
 *
 * This function loops through the stop conditions defined by the user and
 * checks whether any of them is satisfied. If MPI is used, it does so by
 * distributing load on the ranks according to the domain decomposition.
 * This function should be called after store_atom_origin_coordinates so the
 * conditions dependent on the initial state contain the atom coordinates.
 *
 * \param[in]     fplog          The general output file, normally md.log.
 * \param[in]     cr             The struct for communication info.
 * \param[in]     condstop       The conditional stop struct.
 * \param[in]     state_global   The global state (only atom coordinates are used).
 * \param[in]     state_local    The local state (only atom coordinates are used).
 * \returns whether a stop condition was satisfied.
 */
gmx_bool check_stop_conditions(FILE                 *fplog,
                               t_commrec            *cr,
                               const t_condstop     *condstop,
                               const t_state        *state_global,
                               const t_state        *state_local);

#endif
