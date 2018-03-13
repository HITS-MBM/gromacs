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
#include "gmxpre.h"

#include <assert.h>
#include <stdlib.h>
#include <string.h>

#include "gromacs/fileio/readinp.h"
#include "gromacs/gmxpreprocess/readir.h"
#include "gromacs/mdtypes/inputrec.h"
#include "gromacs/utility/cstringutil.h"
#include "gromacs/utility/fatalerror.h"
#include "gromacs/utility/smalloc.h"

char **read_condstopparams(int *ninp_p, t_inpfile **inp_p, t_condstop *condstop,
                           warninp_t wi)
{
    int           ninp, nscan;
    t_inpfile    *inp;
    const char   *tmp;
    char        **grpbuf;
    char          buf[STRLEN], groups[STRLEN], wbuf[STRLEN];
    t_stopcond   *cond;

    ninp   = *ninp_p;
    inp    = *inp_p;

    /* read conditional stop parameters */
    CTYPE("Number of steps between between conditional stop checks");
    ITYPE("conditional-stop-nsteps",      condstop->nsteps, 1000);
    CTYPE("Number of conditional stop groups");
    ITYPE("conditional-stop-ngroups",     condstop->ngrp, 2);

    if (condstop->ngrp < 1)
    {
        gmx_fatal(FARGS, "conditional-stop-ngroups should be >= 1");
    }

    condstop->ngrp += 1;

    snew(condstop->condgrp, condstop->ngrp);

    /* Conditional group options */
    CTYPE("Group parameters");

    /* Read the conditional groups */
    snew(grpbuf, condstop->ngrp);

    /* Group 0 is the absolute reference, we don't read anything for 0 */
    for (int groupNum = 1; groupNum < condstop->ngrp; groupNum++)
    {
        snew(grpbuf[groupNum], STRLEN);
        sprintf(buf, "conditional-stop-group%d-name", groupNum);
        STYPE(buf,              grpbuf[groupNum], "");
    }

    CTYPE("Number of conditional stop conditions");
    ITYPE("conditional-stop-nconds",     condstop->ncond, 1);

    if (condstop->ncond < 1)
    {
        gmx_fatal(FARGS, "conditional-stop-nconds should be >= 1");
    }

    snew(condstop->cond, condstop->ncond);

    /* Stop condition options */
    CTYPE("Stop condition parameters");

    /* Read the stop conditions */
    for (int condNum = 1; condNum < condstop->ncond + 1; condNum++)
    {
        cond = &condstop->cond[condNum - 1];
        sprintf(buf, "conditional-stop-cond%d-ngroups", condNum);
        ITYPE(buf,              cond->ngrp, 2);
        sprintf(buf, "conditional-stop-cond%d-groups", condNum);
        STYPE(buf,              groups, "");

        nscan = sscanf(groups, "%d %d",
                       &cond->grp[0], &cond->grp[1]);
        if (nscan != cond->ngrp)
        {
            sprintf(wbuf, "%s should contain %d conditional group indices",
                    buf, cond->ngrp);
            set_warning_line(wi, NULL, -1);
            warning_error(wi, wbuf);
        }
        for (int g = 0; g < cond->ngrp; g++)
        {
            if (cond->grp[g] < 0 || cond->grp[g] >= condstop->ngrp)
            {
                /* Quit with a fatal error to avoid invalid memory access */
                gmx_fatal(FARGS, "%s contains an invalid pull group %d, you should have %d <= group <= %d",
                          buf, cond->grp[g], 0, condstop->ngrp - 1);
            }
        }

        sprintf(buf, "conditional-stop-cond%d-type", condNum);
        EETYPE(buf,             cond->eType, ecs_names);
        sprintf(buf, "conditional-stop-cond%d-distance-criterion", condNum);
        EETYPE(buf,             cond->eDistCr, ecsdist_names);
        sprintf(buf, "conditional-stop-cond%d-distance", condNum);
        RTYPE(buf,              cond->distance, 0.0);
        sprintf(buf, "conditional-stop-cond%d-origin", condNum);
        EETYPE(buf,             cond->bOrigin, yesno_names);

        if (cond->bOrigin)
        {
            if (cond->grp[0] != 0)
            {
                gmx_fatal(FARGS, "When using the condstop-cond1-origin option, the first group in the condition must be the reference group (group id 0).");
            }
        }
    }

    *ninp_p   = ninp;
    *inp_p    = inp;

    return grpbuf;
}

void make_condstop_groups(t_condstop *condstop, char **cgnames,
                          const t_blocka *grps, char **gnames)
{
    int           ig = -1;
    t_condgrp    *cgrp;

    /* Absolute reference group (might not be used) is special */
    cgrp          = &condstop->condgrp[0];
    cgrp->nat     = 0;

    for (int g = 1; g < condstop->ngrp; g++)
    {
        cgrp = &condstop->condgrp[g];

        if (strcmp(cgnames[g], "") == 0)
        {
            gmx_fatal(FARGS, "Pull option conditional-stop-group%d required by grompp has not been set.", g);
        }

        ig        = search_string(cgnames[g], grps->nr, gnames);
        cgrp->nat = grps->index[ig+1] - grps->index[ig];

        fprintf(stderr, "Conditional stop group %d '%s' has %d atoms\n",
                g, cgnames[g], cgrp->nat);

        if (cgrp->nat == 0)
        {
            gmx_fatal(FARGS, "Conditional stop group %d '%s' is empty", g, cgnames[g]);
        }

        snew(cgrp->ind, cgrp->nat);
        for (int i = 0; i < cgrp->nat; i++)
        {
            cgrp->ind[i] = grps->a[grps->index[ig]+i];
        }
    }
}
