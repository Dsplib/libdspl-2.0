/*
* Copyright (c) 2015-2020 Sergey Bakhurin
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



int DSPL_API find_nearest(double* x, int n, double val, int *idx, double* dist)
{
    double mind, dv;
    int iv, i;
    
    if(!x)
        return ERROR_PTR;
    if(n < 1)
        return ERROR_SIZE;
    
    mind = fabs(x[0] - val);
    iv = 0;
    for(i = 1; i < n; i++)
    {
        dv = fabs(x[i] - val);
        if( dv < mind)
        {
            mind =  dv;
            iv = i;
        }
    }
    
    if(idx)
        *idx = iv;
    if(dist)
        *dist = fabs(x[iv] - val);
    
    return RES_OK;
    
}