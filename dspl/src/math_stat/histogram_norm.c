/*
* Copyright (c) 2015-2022 Sergey Bakhurin
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
#include "dspl.h"





#ifdef DOXYGEN_ENGLISH

#endif
#ifdef DOXYGEN_RUSSIAN

#endif
int DSPL_API histogram_norm(double* y, int n, int nh, double* x, double* w)
{
    double *pedges = NULL;
    int k, res;
    
    if(!y || !x || !w)
        return ERROR_PTR;

    if(n<1 || nh<1)
        return ERROR_SIZE;
    
    pedges = (double*)malloc((nh+1)*sizeof(double));
    
    res = histogram(y, n, nh, pedges, w);
    if(res != RES_OK)
        goto exit_label;

    for(k = 1; k < nh+1; k++)
    {
        x[k-1] = 0.5*(pedges[k] + pedges[k-1]);
        w[k-1] /= ((double)n * (pedges[k] - pedges[k-1]));
    }
    
    res = RES_OK;    
exit_label: 

    if(pedges)
        free(pedges);
    return res;

}


