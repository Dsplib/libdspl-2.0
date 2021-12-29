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
int DSPL_API freqs2time(double* b, double* a, int ord, double fs,
                        int n, fft_t* pfft, double *t, double *h)
{
    double *w = NULL;
    complex_t *hs = NULL;
    complex_t *ht = NULL;
    int err, k;

    if(!b || !a || !t || !h)
        return ERROR_PTR;
    if(ord<1)
        return ERROR_FILTER_ORD;
    if(n<1)
        return ERROR_SIZE;

    w    = (double*)malloc(n*sizeof(double));
    hs = (complex_t*)malloc(n*sizeof(complex_t));


    err = linspace(-fs*0.5, fs*0.5, n, DSPL_PERIODIC, w);
    if(err != RES_OK)
        goto exit_label;

    err = freqs(b, a, ord, w, n, hs);
    if(err != RES_OK)
        goto exit_label;

    err = fft_shift_cmplx(hs, n, hs);
    if(err != RES_OK)
        goto exit_label;

    ht = (complex_t*)malloc(n*sizeof(complex_t));

    err = ifft_cmplx(hs, n, pfft, ht);
    if(err != RES_OK)
    {
        err = idft_cmplx(hs, n, ht);
        if(err != RES_OK)
            goto exit_label;
    }

    for(k = 0; k < n; k++)
    {
        t[k] = (double)k/fs;
        h[k] = RE(ht[k]) * fs;
    }

exit_label:
    if(w)
        free(w);
    if(hs)
        free(hs);
    if(ht)
        free(ht);
    return err;
}
