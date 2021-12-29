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

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "dspl.h"



#ifdef DOXYGEN_ENGLISH

#endif
#ifdef DOXYGEN_RUSSIAN

#endif
double DSPL_API filter_ws1(int ord, double rp, double rs, int type)
{
    double es2, ep2, gs2, x, ws;

    if(ord<1 || rp < 0.0 || rs < 0.0)
        return -1.0;

    es2 = pow(10.0, rs*0.1) - 1.0;
    ep2 = pow(10.0, rp*0.1) - 1.0;
    gs2 = 1.0 / (1.0 + es2);

    x = (1.0 - gs2) / (gs2 * ep2);

    switch( type & DSPL_FILTER_APPROX_MASK)
    {
        case DSPL_FILTER_BUTTER:
            ws = pow(x, 0.5 / (double)ord);
            break;
        case DSPL_FILTER_CHEBY1:
        case DSPL_FILTER_CHEBY2:
            x = sqrt(x) + sqrt(x - 1.0);
            x = log(x) / (double)ord;
            ws    = 0.5 * (exp(-x) + exp(x));
            break;
        case DSPL_FILTER_ELLIP:
        {
            double k, k1;
            complex_t y, z;
            int res;
            k = sqrt(ep2 / es2);
            res = ellip_modulareq(rp, rs, ord, &k1);
            if(res != RES_OK)
            {
                ws = -1.0;
                break;
            }
            RE(z) = sqrt(x);
            IM(z) = 0.0;

            res = ellip_acd_cmplx(&z, 1, k, &y);
            if(res != RES_OK)
            {
                ws = -1.0;
                break;
            }
            RE(y) /= (double)ord;
            IM(y) /= (double)ord;
            res = ellip_cd_cmplx(&y, 1, k1, &z);
            if(res != RES_OK)
            {
                ws = -1.0;
                break;
            }
            ws = RE(z);
            break;
        }
        default:
            ws    = -1.0;
            break;
    }
    return ws;
}

