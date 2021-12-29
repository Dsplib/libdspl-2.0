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
int DSPL_API fft_mag_cmplx(complex_t* x, int n, fft_t* pfft, 
                           double fs, int flag,
                           double* mag, double* freq)
{
    int k, err;
     fft_t *cfft = NULL;
    
    if(pfft)
        cfft = pfft;
    else
    {
        cfft = (fft_t*) malloc(sizeof(fft_t));
        memset(cfft, 0, sizeof(fft_t));
    }
    
    err = fft_abs_cmplx(x, n, cfft, fs, flag, mag, freq);
    if(err != RES_OK)
        goto error_proc;
      
    if(mag)
    {
        if(flag & DSPL_FLAG_LOGMAG)
            for(k = 0; k < n; k++)
                mag[k] = 20.0 * log10(mag[k] + DBL_EPSILON);
        else
            for(k = 0; k < n; k++)
                mag[k] *= mag[k];
    }
error_proc:
    if(cfft && cfft != pfft)
        free(cfft);
    return err;
}

