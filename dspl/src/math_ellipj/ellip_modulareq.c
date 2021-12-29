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
* along with Foobar.    If not, see <http://www.gnu.org/licenses/>.
*/


#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "dspl.h"


#ifdef DOXYGEN_ENGLISH

#endif
#ifdef DOXYGEN_RUSSIAN

#endif
int DSPL_API ellip_modulareq(double rp, double rs, int ord, double *k)
{
    double ep, es, ke, kp, t, sn = 0.0;
    int i, L, r;

    if(rp < 0 || rp == 0)
        return ERROR_FILTER_RP;
    if(rs < 0 || rs == 0)
        return ERROR_FILTER_RS;
    if(ord < 1)
        return ERROR_FILTER_ORD;
    if(!k)
        return ERROR_PTR;


    ep = sqrt(pow(10.0, rp*0.1)-1.0);
    es = sqrt(pow(10.0, rs*0.1)-1.0);

    ke = ep/es;

    ke = sqrt(1.0 - ke*ke);

    r = ord % 2;
    L = (ord-r)/2;

    kp = 1.0;
    for(i = 0; i < L; i++)
    {
        t = (double)(2*i+1) / (double)ord;
        ellip_sn(&t, 1, ke, &sn);
        sn*=sn;
        kp *= sn*sn;
    }

    kp *= pow(ke, (double)ord);
    *k = sqrt(1.0 - kp*kp);

    return RES_OK;

}

