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
#include <string.h>
#include <math.h>
#include "dspl.h"






#ifdef DOXYGEN_ENGLISH

#endif
#ifdef DOXYGEN_RUSSIAN

#endif
int DSPL_API freqs_cmplx(double* b, double* a, int ord,
                         complex_t* s, int n, complex_t *h)
{
    complex_t *bc = NULL;
    complex_t *ac = NULL;
    complex_t num, den;
    double mag;
    int k;
    int res;

    if(!b || !a || !s || !h)
        return ERROR_PTR;
    if(ord<0)
        return ERROR_FILTER_ORD;
    if(n<1)
        return ERROR_SIZE;


    bc = (complex_t*) malloc((ord+1) * sizeof(complex_t));
    res = re2cmplx(b, ord+1, bc);

    if( res!=RES_OK )
        goto exit_label;

    ac = (complex_t*) malloc((ord+1) * sizeof(complex_t));
    res = re2cmplx(a, ord+1, ac);
    if( res!=RES_OK )
        goto exit_label;

    for(k = 0; k < n; k++)
    {
        res = polyval_cmplx(bc, ord, s+k, 1, &num);
        if(res != RES_OK)
            goto exit_label;
        res = polyval_cmplx(ac, ord, s+k, 1, &den);
        if(res != RES_OK)
            goto exit_label;
        mag = ABSSQR(den);
        if(mag == 0.0)
        {
            res = ERROR_DIV_ZERO;
            goto exit_label;
        }
        mag = 1.0 / mag;
        RE(h[k]) = CMCONJRE(num, den) * mag;
        IM(h[k]) = CMCONJIM(num, den) * mag;

    }
    res = RES_OK;
    exit_label:
    if(bc)
        free(bc);
    if(ac)
        free(ac);
    return res;
}
