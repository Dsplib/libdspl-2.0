/*
* Copyright (c) 2015-2018 Sergey Bakhurin
* Digital Signal Processing Library [http://dsplib.org]
*
* This file is part of libdspl-2.0.
*  
* is free software: you can redistribute it and/or modify
* it under the terms of the GNU Lesser  General Public License as published by
* the Free Software Foundation, either version 3 of the License, or
* (at your option) any later version.
*
* DSPL is distributed in the hope that it will be useful,
* but WITHOUT ANY WARRANTY; without even the implied warranty of
* MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
* GNU General Public License for more details.
*
* You should have received a copy of the GNU Lesser General Public License
* along with Foobar.  If not, see <http://www.gnu.org/licenses/>.
*/


#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "dspl.h"



void dscal_(int*, double*, double*, int*);




/**************************************************************************************************
DSCAL scales a vector by a constant
***************************************************************************************************/
int DSPL_API blas_dscal(int n, double a, double* x, int incx)
{
    if(!x)
        return ERROR_PTR;
    if(n < 1)
        return ERROR_SIZE;
    if(incx < 1)
        return ERROR_INC_SIZE;

    dscal_(&n, &a, x, &incx);
    
    return RES_OK;
}









