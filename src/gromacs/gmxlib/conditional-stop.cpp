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
/* This file is completely threadsafe - keep it that way! */
#include "gmxpre.h"

#include "conditional-stop.h"

#include <cmath>
#include <algorithm>
#include <set>
#include <vector>

#include "gromacs/domdec/domdec_struct.h"
#include "gromacs/domdec/ga2la.h"
#include "gromacs/math/vec.h"
#include "gromacs/mdtypes/inputrec.h"
#include "gromacs/mdtypes/state.h"
#include "gromacs/utility/fatalerror.h"
#include "gromacs/utility/smalloc.h"

//! Condition state (unsatisfied, satisfied, or not checked)
enum {
    ecscondUNSAT = FALSE, ecscondSAT = TRUE, ecscondNOCHECK, ecscondNR
};

static const real CONDSTOP_DIST_TOL = 0.0001;

void
compute_coc(const rvec* x, const int n, rvec &coc,
            const char* padding /* = "      " */)
{
    if (debug)
    {
        fprintf(debug, "%sComputing COC of %d atoms\n", padding, n);
        for (int i = 0; i < n; i++)
        {
            fprintf(debug, "%s   atom %d located at %f %f %f\n",
                    padding, i, x[i][0], x[i][1], x[i][2]);
        }
    }
    for (int d = 0; d < DIM; d++)
    {
        coc[d] = 0.0;
        for (int i = 0; i < n; i++)
        {
            coc[d] += x[i][d];
        }
        coc[d] /= n;
    }
    if (debug)
    {
        fprintf(debug, "%sComputed COC at %f %f %f\n", padding, coc[0], coc[1], coc[2]);
    }
}

void
store_condgrp_stateidx(t_condstop* condstop)
{
    std::set<int> indices;
    for (int g = 0; g < condstop->ngrp; g++)
    {
        for (int i = 0; i < condstop->condgrp[g].nat; i++)
        {
            indices.insert(condstop->condgrp[g].ind[i]);
        }
    }

    condstop->nxidx = indices.size();
    snew(condstop->xidx, condstop->nxidx);
    std::copy(indices.begin(), indices.end(), condstop->xidx);
}

void
bcast_condgrp_statex(t_commrec* cr, t_condstop* condstop, t_state* state_global)
{
#if GMX_MPI
    rvec *x_reduced;
    snew(x_reduced, condstop->nxidx);

    if (MASTER(cr))
    {
        for (int i = 0; i < condstop->nxidx; i++)
        {
            int j = condstop->xidx[i];
            for (int d = 0; d < DIM; d++)
            {
                x_reduced[i][d] = state_global->x[j][d];
            }
        }
    }

    MPI_Bcast(x_reduced, condstop->nxidx * sizeof(rvec), MPI_BYTE, MASTERRANK(cr), cr->mpi_comm_mysim);

    if (!MASTER(cr))
    {
        for (int i = 0; i < condstop->nxidx; i++)
        {
            int j = condstop->xidx[i];
            for (int d = 0; d < DIM; d++)
            {
                state_global->x[j][d] = x_reduced[i][d];
            }
        }
    }
#endif
}

gmx_bool
compute_condition_distance(const t_stopcond *cond, const rvec *x1, const int nat1,
                           const rvec *x2, const int nat2, real &distance_squared,
                           const char* padding /* = "      " */)
{
    rvec                        coc1, coc2;
    std::vector<real>::iterator result;
    /* compute distance */
    if (cond->eType == ecsCOCDIST)
    {
        compute_coc(x1, nat1, coc1);
        compute_coc(x2, nat2, coc2);
        distance_squared = distance2(coc1, coc2);
        if (debug)
        {
            fprintf(debug, "%sComputed distance %f between COCs\n%s   group %d COC located at %f %f %f\n%s   group %d COC located at %f %f %f\n",
                    padding, std::sqrt(distance_squared), padding, cond->grp[0], coc1[0], coc1[1], coc1[2],
                    padding, cond->grp[1], coc2[0], coc2[1], coc2[2]);
        }
    }
    else
    {
        std::vector<real> distances_squared;
        for (int j1 = 0; j1 < nat1; j1++)
        {
            for (int j2 = 0; j2 < nat2; j2++)
            {
                real dist2 = distance2(x1[j1], x2[j2]);
                if (debug)
                {
                    fprintf(debug, "%sComputed distance %f between atoms\n%s   group %d atom id %d located at %f %f %f\n%s   group %d atom id %d located at %f %f %f\n",
                            padding, std::sqrt(dist2), padding, cond->grp[0], j1, x1[j1][0], x1[j1][1], x1[j1][2],
                            padding, cond->grp[1], j2, x2[j2][0], x2[j2][1], x2[j2][2]);
                }
                distances_squared.push_back(dist2);
            }
        }
        switch (cond->eType)
        {
            case ecsMINPAIRDIST:
                result = std::min_element(std::begin(distances_squared), std::end(distances_squared));
                break;
            case ecsMAXPAIRDIST:
                result = std::max_element(std::begin(distances_squared), std::end(distances_squared));
                break;
            default:
                gmx_fatal(FARGS, "Unsupported conditional stop type.");
        }
        if (result != std::end (distances_squared))
        {
            distance_squared = *result;
        }
    }

    gmx_bool bSatisfied;
    /* compare the distance with the threshold */
    switch (cond->eDistCr)
    {
        case ecsdistSM:
            bSatisfied = distance_squared < cond->distance * cond->distance;
            break;
        case ecsdistGR:
            bSatisfied = distance_squared > cond->distance * cond->distance;
            break;
        default:
            gmx_fatal(FARGS, "Unsupported conditional stop distance criterion.");
    }

    return bSatisfied;
}

#include "gromacs/domdec/domdec.h"

static gmx_bool
get_atom_coordinates(t_commrec* cr, const t_condstop* condstop,
                     const t_state *state_global, const t_state *state_local,
                     int grpidx, rvec **x, int &nx, const char* padding = "      ")
{
    nx = condstop->condgrp[grpidx].nat;
    snew(*x, nx);
    gmx_bool bLocalGroup = FALSE;
    if (debug)
    {
        fprintf(debug, "%sGetting atom coordinates for group %d with %d atoms\n", padding, grpidx, nx);
    }
    for (int i = 0; i < nx; i++)
    {
        int atom_id_global = condstop->condgrp[grpidx].ind[i];
        int atom_id_local, atom_cell;
        if (DOMAINDECOMP(cr) && ga2la_get(cr->dd->ga2la, atom_id_global,
                                          &atom_id_local, &atom_cell))
        {
            for (int j = 0; j < DIM; j++)
            {
                (*x)[i][j] = state_global->x[atom_id_global][j];
                /* FIXME: this is what should be here instead but is broken until we handle PBC
                 * will give warning: unused parameter ‘state_local’ when compiling
                 * (*x)[i][j] = state_local->x[atom_id_local][j]; */
            }
            bLocalGroup = TRUE;
            if (debug)
            {
                fprintf(debug, "%s   atom %d local id %d global id %d located at %f %f %f\n",
                        padding, i, atom_id_local, atom_id_global, (*x)[i][0], (*x)[i][1], (*x)[i][2]);
            }
        }
        else
        {
            for (int j = 0; j < DIM; j++)
            {
                (*x)[i][j] = state_global->x[atom_id_global][j];
            }
            if (debug)
            {
                fprintf(debug, "%s   atom %d global id %d located at %f %f %f\n",
                        padding, i, atom_id_global, (*x)[i][0], (*x)[i][1], (*x)[i][2]);
            }
        }
    }
    return bLocalGroup;
}

void
store_coord_as_origin(FILE *fplog, t_commrec *cr, t_condstop * condstop,
                      const t_state *state_global, const t_state *state_local)
{
    if (MASTER(cr) && fplog)
    {
        fprintf(fplog, "Conditional stop storing origin coordinates (if any)\n");
    }
    for (int c = 0; c < condstop->ncond; c++)
    {
        if (condstop->cond[c].bOrigin)
        {
            /* if origin = yes only the second group is used */
            get_atom_coordinates(cr, condstop, state_global, state_local,
                                 condstop->cond[c].grp[1], &condstop->cond[c].x,
                                 condstop->cond[c].nx);
            compute_coc(condstop->cond[c].x, condstop->cond[c].nx,
                        condstop->cond[c].coc);
            if (MASTER(cr) && fplog)
            {
                fprintf(fplog, "   Stored origin for stop condition %d group %d (type %s, distance criterion %s, distance threshold %f)\n",
                        c, condstop->cond[c].grp[1], ECONDSTOPTYPE(condstop->cond[c].eType), ECONDSTOPDISTCR(condstop->cond[c].eDistCr), condstop->cond[c].distance);
            }
        }
    }
    if (MASTER(cr) && fplog)
    {
        fprintf(fplog, "\n");
    }
}

gmx_bool
check_stop_conditions(FILE *fplog, t_commrec *cr, const t_condstop *condstop,
                      const t_state *state_global, const t_state *state_local)
{
    gmx_bool bStop      = FALSE;
    int      nnodes     = 1;
    int     *cond, *cond_local;
    real    *cond_dist, *cond_dist_local;

    snew(cond_local, condstop->ncond);
    snew(cond_dist_local, condstop->ncond);

    if (DOMAINDECOMP(cr))
    {
        nnodes = cr->nnodes;
        snew(cond, nnodes * condstop->ncond);
        snew(cond_dist, nnodes * condstop->ncond);
    }
    else
    {
        /* convenient alias to avoid separate code paths for DD vs non-DD */
        cond      = cond_local;
        cond_dist = cond_dist_local;
    }

    if (debug)
    {
        fprintf(debug, "Conditional stop started checking conditions\n");
    }
    for (int c = 0; c < condstop->ncond; c++)
    {
        int                          grp1idx        = condstop->cond[c].grp[0];
        int                          grp2idx        = condstop->cond[c].grp[1];
        gmx_bool                     bLocalGroup1   = FALSE, bLocalGroup2 = FALSE;
        rvec                        *x1, *x2;
        int                          nat1, nat2;
        real                         distance_squared = 0.0;

        if (debug)
        {
            if (condstop->cond[c].bOrigin)
            {
                fprintf(debug, "   Checking stop condition %d type %s containing group %d and its origin\n", c, ECONDSTOPTYPE(condstop->cond[c].eType), condstop->cond[c].grp[1]);
            }
            else
            {
                fprintf(debug, "   Checking stop condition %d type %s containing groups %d and %d\n", c, ECONDSTOPTYPE(condstop->cond[c].eType), condstop->cond[c].grp[0], condstop->cond[c].grp[1]);
            }
        }

        /* get atom coordinates */
        if (condstop->cond[c].bOrigin)
        {
            x1           = condstop->cond[c].x;
            nat1         = condstop->cond[c].nx;
            bLocalGroup1 = TRUE;
        }
        else
        {
            bLocalGroup1 = get_atom_coordinates(cr, condstop, state_global,
                                                state_local, grp1idx, &x1, nat1);
        }
        bLocalGroup2 = get_atom_coordinates(cr, condstop, state_global,
                                            state_local, grp2idx, &x2, nat2);

        /* break if no atoms from either group are in the local domain */
        if (DOMAINDECOMP(cr) && !bLocalGroup1 && !bLocalGroup2)
        {
            cond_local[c] = ecscondNOCHECK;
            break;
        }

        /* compute the distance and compare with the threshold */
        cond_local[c]      = compute_condition_distance(&condstop->cond[c], x1, nat1, x2, nat2, distance_squared);
        cond_dist_local[c] = sqrt(distance_squared);

        if (debug)
        {
            if (cond_local[c])
            {
                fprintf(debug, "      Satisfied ");
            }
            else
            {
                fprintf(debug, "      Unsatisfied ");
            }
            fprintf(debug, "with distance %f which is ", cond_dist_local[c]);
            if (!cond_local[c])
            {
                fprintf(debug, "not ");
            }
            fprintf(debug, "%s than distance threshold %f\n", ECONDSTOPDISTCR(condstop->cond[c].eDistCr), condstop->cond[c].distance);
        }
    }

    /* gather the conditions */
    if (DOMAINDECOMP(cr))
    {
#if GMX_MPI
        MPI_Gather(cond_local, condstop->ncond, MPI_INT, cond, condstop->ncond, MPI_INT, MASTERRANK(cr), cr->mpi_comm_mysim);
        MPI_Gather(cond_dist_local, condstop->ncond * sizeof(real), MPI_BYTE, cond_dist, condstop->ncond * sizeof(real), MPI_BYTE, MASTERRANK(cr), cr->mpi_comm_mysim);
#endif
    }

    if (MASTER(cr))
    {
        int nsatcond = 0;

        if (fplog)
        {
            fprintf(fplog, "Conditional stop finished checking conditions\n");
        }
        for (int c = 0; c < condstop->ncond; c++)
        {
            int nsat = 0, nunsat = 0, nnocheck = 0;
            for (int n = 0; n < nnodes; n++)
            {
                switch (cond[n * condstop->ncond + c])
                {
                    case ecscondNOCHECK:
                        nnocheck++;
                        break;
                    case ecscondSAT:
                        nsat++;
                        break;
                    case ecscondUNSAT:
                        nunsat++;
                        break;
                    default:
                        gmx_fatal(FARGS, "Unsupported conditional stop condition state.");
                }
            }


            if (nsat == 0 && nunsat == 0 && nnocheck > 0)
            {
                gmx_fatal(FARGS, "Conditional stop condition %d state not checked. Please report this bug.", c);
            }
            else if (nsat > 0 && nunsat > 0)
            {
                gmx_fatal(FARGS, "Inconsistent conditional stop condition %d states. Please report this bug.", c);
            }
            else
            {
                std::vector<real> distances;
                for (int n = 0; n < nnodes; n++)
                {
                    switch (cond[n * condstop->ncond + c])
                    {
                        case ecscondNOCHECK:
                            break;
                        case ecscondSAT:
                        case ecscondUNSAT:
                            distances.push_back(cond_dist[n * condstop->ncond + c]);
                            break;
                        default:
                            gmx_fatal(FARGS, "Unsupported conditional stop condition state.");
                    }
                }

                for (auto &dist : distances)
                {
                    if (std::abs(dist - distances[0]) > CONDSTOP_DIST_TOL)
                    {
                        gmx_fatal(FARGS, "Conditional stop computed condition %d distances %f and %f that differ above numerical tolerance. Please report this bug.", c, distances[0], dist);
                    }
                }

                if (nsat > 0)
                {
                    bStop = TRUE;
                    nsatcond++;
                }

                if (fplog)
                {
                    if (nsat > 0)
                    {
                        fprintf(fplog, "   Satisfied ");
                    }
                    else
                    {
                        fprintf(fplog, "   Unsatisfied ");
                    }
                    fprintf(fplog, "stop condition %d containing groups %d and %d with %s distance %f ",
                            c, condstop->cond[c].grp[0], condstop->cond[c].grp[1],
                            ECONDSTOPTYPE(condstop->cond[c].eType), distances[0]);
                    if (nsat == 0)
                    {
                        fprintf(fplog, "not ");
                    }
                    fprintf(fplog, "%s than %f\n", ECONDSTOPDISTCR(condstop->cond[c].eDistCr),
                            condstop->cond[c].distance);
                }
            }
        }
        if (fplog)
        {
            fprintf(fplog, "Conditional stop found %d satisfied and %d unsatisfied stop conditions (%d total)\n\n", nsatcond, condstop->ncond - nsatcond, condstop->ncond);
        }
    }

    /* broadcast whether to stop */
    if (DOMAINDECOMP(cr))
    {
#if GMX_MPI
        MPI_Bcast(&bStop, 1, MPI_INT, MASTERRANK(cr), cr->mpi_comm_mysim);
#endif
    }

    return bStop;
}
