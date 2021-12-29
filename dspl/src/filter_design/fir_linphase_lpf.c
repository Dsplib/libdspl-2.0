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

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "dspl.h"
#include "dspl_internal.h"



/******************************************************************************
 * Linear phase lowpass filter
 ******************************************************************************/
int fir_linphase_lpf(int ord, double wp, int win_type,
                     double win_param, double* h)
{
    int n, err = RES_OK;
    double *w = NULL;


    w = (double*)malloc((ord+1)*sizeof(double));

    err = linspace(-(double)ord*0.5, (double)ord*0.5, ord+1, DSPL_SYMMETRIC, w);

    if(err!=RES_OK)
        goto error_proc;

    err = sinc(w, ord+1, M_PI*wp, h);

    if(err!=RES_OK)
        goto error_proc;

    err = window(w, ord+1, win_type | DSPL_SYMMETRIC, win_param);

    if(err!=RES_OK)
        goto error_proc;

    for(n = 0; n < ord+1; n++)
        h[n] *= w[n] * wp;

error_proc:
    if(w)
        free(w);
    return err;
}


