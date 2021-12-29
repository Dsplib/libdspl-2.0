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
int DSPL_API poly_z2a_cmplx(complex_t* z, int nz, int ord, complex_t* a)
{
    int k, ind, res;
    complex_t x[2];

    if(!z || !a)
        return ERROR_PTR;
    if(nz < 0)
        return ERROR_SIZE;
    if(nz > ord || ord < 1)
        return ERROR_POLY_ORD;

    RE(x[1]) = 1.0;
    IM(x[1]) = 0.0;

    memset(a, 0, (ord+1) * sizeof(complex_t));

    RE(a[0]) = 1.0;
    ind = 1;
    for(k = 0; k < nz; k++)
    {
        RE(x[0]) = -RE(z[k]);
        IM(x[0]) = -IM(z[k]);
        res = conv_cmplx(a, ind, x, 2, a);
        if(res!=RES_OK)
            return res;
        ind++;
    }

    return RES_OK;
}


