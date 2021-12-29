/*
* Copyright (c) 2015-2019 Sergey Bakhurin
* Digital Signal Processing Library [http://dsplib.org]
*
* This file is part of DSPL.
*
* is free software: you can redistribute it and/or modify
* it under the terms of the GNU General Public License as published by
* the Free Software Foundation, either version 3 of the License, or
* (at your option) any later version.
*
* DSPL is distributed in the hope that it will be useful,
* but WITHOUT ANY WARRANTY; without even the implied warranty of
* MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.    See the
* GNU General Public License for more details.
*
* You should have received a copy of the GNU General Public License
* along with Foobar.    If not, see <http://www.gnu.org/licenses/>.
*/


#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "dspl.h"



#ifdef DOXYGEN_ENGLISH

#endif
#ifdef DOXYGEN_RUSSIAN

#endif
int DSPL_API ellip_rat(double* w, int n, int ord, double k, double* u)
{
    double t, xi, w2, xi2, k2;
    int i, m, r, L;

    if(!u || !w)
        return ERROR_PTR;
    if(n<1)
        return ERROR_SIZE;
    if(k < 0.0 || k>= 1.0)
        return ERROR_ELLIP_MODULE;

    r = ord%2;
    L = (ord-r)/2;

    if(r)
        memcpy(u, w, n*sizeof(double));
    else
    {
        for(m = 0; m < n; m++)
        {
            u[m] = 1.0;
        }
    }

    k2 = k*k;
    for(i = 0; i < L; i++)
    {
        t = (double)(2*i+1) / (double)ord;
        ellip_cd(&t, 1, k, &xi);
        xi2 = xi*xi;
        for(m = 0; m < n; m++)
        {
            w2 = w[m]*w[m];
            u[m] *= (w2 - xi2) / (1.0 - w2 * k2 * xi2);
            u[m] *= (1.0 - k2*xi2) / (1.0 - xi2);
        }
    }
    return RES_OK;
}


