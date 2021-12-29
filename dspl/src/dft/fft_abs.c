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
int DSPL_API fft_abs(double* x, int n, fft_t* pfft, 
                     double fs, int flag,
                     double* mag, double* freq)
{
    int k, err = RES_OK;
    complex_t *X = NULL;
    if(!x || !pfft)
        return ERROR_PTR;

    if(n<1)
        return ERROR_SIZE;

    if(mag)
    {    
        X = (complex_t*)malloc(n*sizeof(complex_t));
        err = fft(x, n, pfft, X);
        if(err!=RES_OK)
            goto error_proc;
        

        for(k = 0; k < n; k++)
            mag[k] = ABS(X[k]);
        if(flag & DSPL_FLAG_FFT_SHIFT)
        {
            err = fft_shift(mag, n, mag);
            if(err!=RES_OK)
                goto error_proc;
        }
    }
    
    if(freq)
    {
        if(flag & DSPL_FLAG_FFT_SHIFT)
            if(n%2)
                err = linspace(-fs*0.5 + fs*0.5/(double)n, 
                                fs*0.5 - fs*0.5/(double)n, 
                                n, DSPL_SYMMETRIC, freq);
            else
                err = linspace(-fs*0.5, fs*0.5, n, DSPL_PERIODIC, freq);
        else
            err = linspace(0, fs, n, DSPL_PERIODIC, freq);
    } 

error_proc:
    if(X)
        free(X);

    return err;
}

