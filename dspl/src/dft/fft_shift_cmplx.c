/*
* Copyright (c) 2015-2019 Sergey Bakhurin
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
#include <stdio.h>
#include <string.h>
#include <float.h>

#include "dspl.h"
#include "dspl_internal.h"


#ifdef DOXYGEN_ENGLISH

#endif
#ifdef DOXYGEN_RUSSIAN

#endif
int DSPL_API fft_shift_cmplx(complex_t* x, int n, complex_t* y)
{
    int n2, r;
    int k;
    complex_t tmp;
    complex_t *buf;

    if(!x || !y)
        return ERROR_PTR;

    if(n<1)
        return ERROR_SIZE;

    r = n%2;
    if(!r)
    {
        n2 = n>>1;
        for(k = 0; k < n2; k++)
        {
            RE(tmp) = RE(x[k]);
            IM(tmp) = IM(x[k]);

            RE(y[k]) = RE(x[k+n2]);
            IM(y[k]) = IM(x[k+n2]);

            RE(y[k+n2]) = RE(tmp);
            IM(y[k+n2]) = IM(tmp);
        }
    }
    else
    {
        n2 = (n+1) >> 1;
        buf = (complex_t*) malloc(n2*sizeof(complex_t));
        memcpy(buf, x, n2*sizeof(complex_t));
        memcpy(y, x+n2, (n2-1)*sizeof(complex_t));
        memcpy(y+n2-1, buf, n2*sizeof(complex_t));
        free(buf);
    }
    return RES_OK;
}

