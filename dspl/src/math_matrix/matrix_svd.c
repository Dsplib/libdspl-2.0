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
int DSPL_API matrix_svd(double* a, int n, int m, 
                        double* u, double* s, double* vt, int* info)
{
    char jobz = 'A';
    double* work = NULL;
    int* iwork = NULL;
    int lwork, mn, mx, err;
    int pi;
    
    if(!a || !u || !s || !vt)
        return ERROR_PTR;
    if(n < 1 || m < 1)
        return ERROR_SIZE;

    if(n > m)
    {
        mn = m;
        mx = n;
    }
    else
    { 
        mn = n; 
        mx = m;
    }

    err = RES_OK;
    
    lwork = 4 * mn * mn + 6 * mn + mx;
    work = (double*) malloc(lwork*sizeof(double));
    iwork = (int*) malloc(8*mn*sizeof(int));
    dgesdd_(&jobz, &n, &m, a, &n, s, u, &n, vt, &m, work, &lwork, iwork, &pi);
    
    if(info)
        *info = pi;
      
    if(pi)
        err = ERROR_LAPACK;
    if(work)
        free(work);
    if(iwork)
        free(iwork);
    return err;
}
