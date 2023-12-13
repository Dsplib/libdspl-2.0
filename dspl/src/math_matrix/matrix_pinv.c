/*
* Copyright (c) 2015-2024 Sergey Bakhurin
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
int DSPL_API matrix_pinv(double* a, int n, int m, double* tol, 
                         double* inv, int* info)
{
    int err, mn, i, j;
    double eps;


    double* u  = NULL;
    double* vt = NULL;
    double* v  = NULL;
    double* ut = NULL;
    double* s  = NULL;
    
    
    mn = (m > n) ? n : m;
    
    u = (double*) malloc(n*n*sizeof(double));
    if(!u)
    {
        err = ERROR_MALLOC;
        goto exit_label;
    }
    
    vt = (double*) malloc(m*m*sizeof(double));
    if(!vt)
    {
        err = ERROR_MALLOC;
        goto exit_label;
    }
    
    s = (double*) malloc(mn*sizeof(double));
    if(!s)
    {
        err = ERROR_MALLOC;
        goto exit_label;
    }
    
    err = matrix_svd(a, n, m, u, s, vt, info);
    if(err != RES_OK)
        goto exit_label;
    
    
      
    if(tol)
        eps = *tol;
    else
    {
        double smax;
        double mx = (n > m) ? (double)n : (double)m;
        err = minmax(s, mn, NULL, &smax);
        eps = DBL_EPSILON * mx * smax;
    }
    
    for(i = 0; i < mn; i++)
        if(s[i] > eps)
            s[i] = 1.0 / s[i];
        else
            s[i] = 0.0;

    v = (double*) malloc(m*m*sizeof(double));
    if(!v)
    {
        err = ERROR_MALLOC;
        goto exit_label;
    }
    err = matrix_transpose(vt, m, m, v);
    if(err != RES_OK)
        goto exit_label;
    
    ut = (double*) malloc(n*n*sizeof(double));
    if(!ut)
    {
        err = ERROR_MALLOC;
        goto exit_label;
    }
    err = matrix_transpose(u, n, n, ut);
    if(err != RES_OK)
        goto exit_label;
      
    for(i = 0; i < mn; i++)
        for(j = 0; j < m; j++)
            v[j + i*m] *= s[i];
    
    if(mn < m)
        memset(v+ mn*m, 0, (m-mn)*sizeof(double));
    
    err = matrix_mul(v, m, n, ut, n, n, inv);

exit_label:
    if(u)
        free(u);
    if(vt)
        free(vt);
    if(s)
        free(s);
    if(v)
        free(v);
    if(ut)
        free(ut);
    
    return err;
}


