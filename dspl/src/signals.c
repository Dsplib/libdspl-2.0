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
#include <string.h>
#include <math.h>
#include "dspl.h"


#ifdef DOXYGEN_ENGLISH

#endif
#ifdef DOXYGEN_RUSSIAN

#endif
int    DSPL_API signal_pimp(double* t, size_t n, double amp,
                            double tau, double dt, double period, double* y)
{
    size_t k;
    double ll, lr, p2, tp;

    if(!t || !y)
        return ERROR_PTR;
    if(n < 1)
        return ERROR_SIZE;
    if(tau < 0.0 || period < 0.0)
        return ERROR_NEGATIVE;


    ll = -0.5 * tau;
    lr =    0.5 * tau;
    p2 = period*0.5;
    for(k = 0; k < n; k++)
    {
        tp = dmod(t[k] - dt + p2, period) - p2;
        y[k] = (tp < ll || tp > lr) ? 0.0 : amp;
    }
    return RES_OK;
}




#ifdef DOXYGEN_ENGLISH

#endif
#ifdef DOXYGEN_RUSSIAN

#endif
int    DSPL_API signal_saw(double* t, size_t n, double amp,
                           double dt, double period, double* y)
{
    size_t k;
    double p2, tp;

    if(!t || !y)
        return ERROR_PTR;
    if(n < 1)
        return ERROR_SIZE;
    if(period < 0.0)
        return ERROR_NEGATIVE;

    p2 = period*0.5;
    for(k = 0; k < n; k++)
    {
        tp = dmod(t[k] - dt + p2, period) - p2;
        y[k] = amp * tp;
    }
    return RES_OK;
}

