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
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <float.h>
#include "dspl.h"

#include "blas.h"


#ifdef DOXYGEN_ENGLISH

#endif
#ifdef DOXYGEN_RUSSIAN

#endif
int DSPL_API matrix_transpose_hermite(complex_t* a, int n, int m, complex_t* b)
{
    int p, q, i, j, aind, bind;
    
    if(!a || !b)
        return ERROR_PTR;
    if(n < 1 || m < 1)
        return ERROR_MATRIX_SIZE;
    
    for(p = 0; p < n - DSPL_MATRIX_BLOCK; p+=DSPL_MATRIX_BLOCK)
    {
        for(q = 0; q < m - DSPL_MATRIX_BLOCK; q+=DSPL_MATRIX_BLOCK)
        {
            for(i = 0; i < DSPL_MATRIX_BLOCK; i++)
            {
                for(j = 0; j < DSPL_MATRIX_BLOCK; j++)
                {
                    aind = (q+j) * n + p + i;
                    bind = (p+i) * m + q + j;
                    RE(b[bind]) =    RE(a[aind]);
                    IM(b[bind]) = -IM(a[aind]);
                }
            }
        }
    }
    for(i = p; i < n; i++)
    {
        for(j = 0; j < m; j++)
        {
            RE(b[i*m + j]) = RE(a[j*n+i]);
            IM(b[i*m + j]) = -IM(a[j*n+i]);
        }
    }

    for(i = 0; i < p; i++)
    {
        for(j = q; j < m; j++)
        {
            RE(b[i*m + j]) =    RE(a[j*n+i]);
            IM(b[i*m + j]) = -IM(a[j*n+i]);
        }
    }
    
    return RES_OK;
}

