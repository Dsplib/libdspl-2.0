/*
* Copyright (c) 2015-2019 Sergey Bakhurin
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
* along with Foobar. If not, see <http://www.gnu.org/licenses/>.
*/


#include <stdio.h>
#include <stdlib.h>
#include "dspl.h"


#ifdef DOXYGEN_ENGLISH

#endif
#ifdef DOXYGEN_RUSSIAN

#endif
int DSPL_API trapint(double* x, double* y, int n, double* sum)
{
    int k;

    if(!x || !y)
        return ERROR_PTR;
    if(n<2)
        return ERROR_SIZE;
    *sum = 0.0;

    for(k = 1; k < n; k++)
        *sum += 0.5 * (x[k] - x[k-1]) * (y[k] + y[k-1]);

    return RES_OK;
}



#ifdef DOXYGEN_ENGLISH

#endif
#ifdef DOXYGEN_RUSSIAN

#endif
int DSPL_API trapint_cmplx(double* x, complex_t* y, int n, complex_t* sum)
{
    int k;
    double dx;
    if(!x || !y)
        return ERROR_PTR;
    if(n<2)
        return ERROR_SIZE;
    RE(*sum) = IM(*sum) = 0.0;

    for(k = 1; k < n; k++)
    {
        dx = 0.5 * (x[k] - x[k-1]);
        RE(*sum) +=    dx * (RE(y[k]) + RE(y[k-1]));
        IM(*sum) +=    dx * (IM(y[k]) + IM(y[k-1]));
    }
    return RES_OK;
}

