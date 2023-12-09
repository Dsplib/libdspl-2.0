/*
* Copyright (c) 2015-2022 Sergey Bakhurin
* Digital Signal Processing Library [http://dsplib.org]
*
* This file is part of libdspl-2.0.
*
* is free software: you can redistribute it and/or modify
* it under the terms of the GNU Lesser    General Public License as published by
* the Free Software Foundation, either version 3 of the License, or
* (at your option) any later version.
*
* DSPL is distributed in the hope that it will be useful,
* but WITHOUT ANY WARRANTY; without even the implied warranty of
* MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.    See the
* GNU General Public License for more details.
*
* You should have received a copy of the GNU Lesser General Public License
* along with Foobar.    If not, see <http://www.gnu.org/licenses/>.
*/


#include <stdlib.h>
#include <math.h>
#include <time.h>

#include "dspl.h"
#include "randomgen.h"



#ifdef DOXYGEN_ENGLISH

#endif
#ifdef DOXYGEN_RUSSIAN

#endif
int randu_mrg32k3a (double* u, int n, random_t* prnd)
{

    if(!u || !prnd)
        return ERROR_PTR;
    if(n < 1)
        return ERROR_SIZE;

    long z;
    double xn, yn, *x, *y;
    int k;

    x = prnd->mrg32k3a_x;
    y = prnd->mrg32k3a_y;
    for(k = 0; k < n; k++)
    {
        /* Component x[n] */
        xn = MRG32K3A_A12 * x[1] - MRG32K3A_A13 * x[2];

        z = (long)(xn / MRG32K3A_M1);
        xn -= (double)z * MRG32K3A_M1;
        if (xn < 0.0)
            xn += MRG32K3A_M1;

        x[2] = x[1];
        x[1] = x[0];
        x[0] = xn;

        /* Component y[n] */
        yn = MRG32K3A_A21 * y[0] - MRG32K3A_A23 * y[2];
        z = (long)(yn / MRG32K3A_M2);
        yn -= (double)z * MRG32K3A_M2;
        if (yn < 0.0)
             yn += MRG32K3A_M2;

        y[2] = y[1];
        y[1] = y[0];
        y[0] = yn;

        /* Combination */
        u[k] = (xn <= yn) ? ((xn - yn + MRG32K3A_M1) * MRG32K3A_NORM):
                            (xn - yn) * MRG32K3A_NORM;
    }
    return RES_OK;
}




