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
int DSPL_API matrix_print_cmplx(complex_t* a, int n, int m, 
                                const char* name, const char* format)
{
    int p,q;

    if(!a)
        return ERROR_PTR;
    if(n < 1 || m < 1)
        return ERROR_MATRIX_SIZE;

    if(!a)
        return ERROR_PTR;
    if(n < 1 || m < 1)
        return ERROR_SIZE;
        
    printf("\n%s = [ %% size [%d x %d] type: complex", name, n, m);
    
    for(p = 0; p < n; p++)
    {
        printf("\n");
        for(q = 0; q < m; q++)
        {
            printf(format, RE(a[q*n + p]), IM(a[q*n + p]));
            if(q == m-1)
                printf(";");
            else
                printf(", ");
        }
    }
    printf("];\n");

    return RES_OK;
}

