/*
 * This file is part of the GROMACS molecular simulation package.
 *
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
/*! \internal \file
 * \brief
 * Tests for conditional stop.
 *
 * \author Vedran Miletić <vedran@miletic.net>
 * \ingroup module_gmxlib
 */
#include "gmxpre.h"

#include <cmath>

#include <algorithm>

#include <gtest/gtest.h>

#include "gromacs/math/vec.h"
#include "gromacs/gmxlib/conditional-stop.h"
#include "gromacs/mdtypes/commrec.h"
#include "gromacs/mdtypes/inputrec.h"
#include "gromacs/mdtypes/state.h"
#include "gromacs/utility/smalloc.h"

#include "testutils/refdata.h"
#include "testutils/testasserts.h"
#include "testutils/testfilemanager.h"

namespace gmx
{
namespace
{

using gmx::test::defaultRealTolerance;

const int                   natmax = 4;

class ConditionalStopTest : public ::testing::Test
{
    protected:
        ConditionalStopTest() {}

        void test(rvec x[], int ind[][natmax], int nat[], int dist_criterion[], real dist_threshold[])
        {
            rvec       xgrpatom[2][natmax];
            t_condstop condstop;
            condstop.ngrp  = 2;
            condstop.ncond = 1;

            snew(condstop.condgrp, condstop.ngrp);
            for (int g = 0; g < 2; g++)
            {
                t_condgrp *cgrp = &condstop.condgrp[g];
                cgrp->nat       = nat[g];
                snew(cgrp->ind, cgrp->nat);
                for (int i = 0; i < cgrp->nat; i++)
                {
                    cgrp->ind[i] = ind[g][i];
                    for (int d = 0; d < DIM; d++)
                    {
                        xgrpatom[g][i][d] = x[ind[g][i]][d];
                    }
                }
            }

            gmx_bool bCond   = FALSE;
            real     dist_sq = 1.0;

            snew(condstop.cond, condstop.ncond);
            t_stopcond *cond = &condstop.cond[0];
            cond->ngrp = 2;

            // trivial test of origin = yes
            cond->grp[0]   = 0;
            cond->grp[1]   = 0;
            cond->bOrigin  = TRUE;
            cond->x        = xgrpatom[0];
            cond->nx       = condstop.condgrp[0].nat;
            cond->eType    = ecsCOCDIST;
            cond->eDistCr  = ecsdistSM;
            cond->distance = 0.1;
            bCond          = compute_condition_distance(cond, cond->x, cond->nx,
                                                        xgrpatom[0], condstop.condgrp[0].nat, dist_sq);
            EXPECT_TRUE(bCond);
            cond->bOrigin = FALSE;
            cond->x       = nullptr;
            cond->nx      = 0;

            cond->grp[0]   = 0;
            cond->grp[1]   = 1;
            cond->eType    = ecsMINPAIRDIST;
            cond->eDistCr  = dist_criterion[0];
            cond->distance = dist_threshold[0];
            bCond          = compute_condition_distance(cond, xgrpatom[0], condstop.condgrp[0].nat,
                                                        xgrpatom[1], condstop.condgrp[1].nat, dist_sq);
            EXPECT_FALSE(bCond);

            cond->eDistCr  = dist_criterion[1];
            cond->distance = dist_threshold[1];
            bCond          = compute_condition_distance(cond, xgrpatom[0], condstop.condgrp[0].nat,
                                                        xgrpatom[1], condstop.condgrp[1].nat, dist_sq);
            EXPECT_TRUE(bCond);

            cond->eType    = ecsMAXPAIRDIST;
            cond->eDistCr  = dist_criterion[2];
            cond->distance = dist_threshold[2];
            bCond          = compute_condition_distance(cond, xgrpatom[0], condstop.condgrp[0].nat,
                                                        xgrpatom[1], condstop.condgrp[1].nat, dist_sq);
            EXPECT_FALSE(bCond);

            cond->eDistCr  = dist_criterion[3];
            cond->distance = dist_threshold[3];
            bCond          = compute_condition_distance(cond, xgrpatom[0], condstop.condgrp[0].nat,
                                                        xgrpatom[1], condstop.condgrp[1].nat, dist_sq);
            EXPECT_TRUE(bCond);

            cond->eType    = ecsCOCDIST;
            cond->eDistCr  = dist_criterion[4];
            cond->distance = dist_threshold[4];
            bCond          = compute_condition_distance(cond, xgrpatom[0], condstop.condgrp[0].nat,
                                                        xgrpatom[1], condstop.condgrp[1].nat, dist_sq);
            EXPECT_FALSE(bCond);

            cond->eDistCr  = dist_criterion[5];
            cond->distance = dist_threshold[5];
            bCond          = compute_condition_distance(cond, xgrpatom[0], condstop.condgrp[0].nat,
                                                        xgrpatom[1], condstop.condgrp[1].nat, dist_sq);
            EXPECT_TRUE(bCond);

            for (int g = 0; g < condstop.ngrp; g++)
            {
                t_condgrp *cgrp = &condstop.condgrp[g];
                sfree(cgrp->ind);
            }
            sfree(condstop.condgrp);
            for (int c = 0; c < condstop.ncond; c++)
            {
                t_stopcond *cond = &condstop.cond[c];
                sfree(cond->x);
            }
            sfree(condstop.cond);
        }
};

TEST_F (ConditionalStopTest, 2Groups1Condition)
{
    rvec x[] = { { 1.0, 2.0, 3.0 },
                 { 2.0, 4.0, 6.0 },
                 { 3.0, 6.0, 9.0 },
                 { 15.0, 2.0, 10.0 },
                 { 12.0, 4.0, 7.0 },
                 { 18.0, 3.0, 14.0 },
                 { 27.0, 5.0, 8.0 } };
    int  ind[][natmax] = { { 0, 1, 2, -1 },
                           { 3, 4, 5, 6 } };
    int  nat[]            = { 3, 4 };
    int  dist_criterion[] = { ecsdistSM, ecsdistSM, ecsdistGR, ecsdistGR, ecsdistSM, ecsdistSM };
    // squared distances are 89.0, 710.0, 270.3125
    real dist_threshold[] = { 9.0, 9.5, 26.7, 26.5, 16.3, 16.6 };

    test(x, ind, nat, dist_criterion, dist_threshold);
}

}

}
