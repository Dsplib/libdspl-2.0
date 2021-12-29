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
int DSPL_API fourier_series_dec_cmplx(double* t, complex_t* s, int nt,
            double period, int nw, double* w, complex_t* y)
{
    int k, m;
    double dw = M_2PI / period;
    complex_t e[2];

    if(!t || !s || !w || !y)
        return ERROR_PTR;
    if(nt<1 || nw < 1)
        return ERROR_SIZE;
    if(period <= 0.0)
        return ERROR_NEGATIVE;

    memset(y, 0 , nw*sizeof(complex_t));

    for(k = 0; k < nw; k++)
    {
        w[k] = (k - nw/2) * dw;
        RE(e[1]) =    RE(s[0]) * cos(w[k] * t[0]) +
        IM(s[0]) * sin(w[k] * t[0]);
        IM(e[1]) = -RE(s[0]) * sin(w[k] * t[0]) +
        IM(s[0]) * cos(w[k] * t[0]);
        for(m = 1; m < nt; m++)
        {
            RE(e[0]) = RE(e[1]);
            IM(e[0]) = IM(e[1]);
            RE(e[1]) =     RE(s[m]) * cos(w[k] * t[m]) +
            IM(s[m]) * sin(w[k] * t[m]);
            IM(e[1]) = -RE(s[m]) * sin(w[k] * t[m]) +
            IM(s[m]) * cos(w[k] * t[m]);
            RE(y[k]) += 0.5 * (RE(e[0]) + RE(e[1]))*(t[m] - t[m-1]);
            IM(y[k]) += 0.5 * (IM(e[0]) + IM(e[1]))*(t[m] - t[m-1]);
        }
        RE(y[k]) /= period;
        IM(y[k]) /= period;
    }

    if(!(nw%2))
        RE(y[0]) = RE(y[1]) = 0.0;

    return RES_OK;
}
