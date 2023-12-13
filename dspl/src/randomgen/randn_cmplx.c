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


#include <stdlib.h>
#include <math.h>
#include <time.h>

#include "dspl.h"
#include "randomgen.h"



#ifdef DOXYGEN_ENGLISH

#endif
#ifdef DOXYGEN_RUSSIAN

#endif
int DSPL_API randn_cmplx(complex_t* x, int n, complex_t* mu, 
                         double sigma, random_t* prnd)
{
    int err, i;
    
    err = randn((double*)x, 2*n, 0.0, sigma / M_SQRT2, prnd);
    if(err!= RES_OK)
        return err;
    if(mu)
    {
        for(i = 0; i < n; i++)
        {
            RE(x[i]) += RE(mu[0]);
            IM(x[i]) += IM(mu[0]);
        }
    }
    return err;
}

