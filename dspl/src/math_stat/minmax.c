/*
* Copyright (c) 2015-2024 Sergey Bakhurin
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
int DSPL_API minmax(double* x, int n, double* xmin, double* xmax)
{
    int k; 
    double min, max;
    
    if(!x)
        return ERROR_PTR;

    if(n<1)
        return ERROR_SIZE;

    min = max = x[0];
    for(k = 1; k < n; k++)
    {
        min = x[k] < min ? x[k] : min;
        max = x[k] > max ? x[k] : max;
    }
    
    if(xmin)
        *xmin = min;
    if(xmax)
        *xmax = max;
    
    return RES_OK;
}


