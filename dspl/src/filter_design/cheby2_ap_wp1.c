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
int DSPL_API cheby2_ap_wp1(double rp, double rs, int ord, double* b, double* a)
{
    int err;
    double es, gp, alpha, beta, y, wp;

    if(rp <= 0)
        return    ERROR_FILTER_RP;

    err = cheby2_ap(rs, ord, b, a);
    if(err!=RES_OK)
        goto exit_label;

    es = sqrt(pow(10.0, rs*0.1) - 1.0);
    gp = pow(10.0, -rp*0.05);
    alpha = gp * es / sqrt(1.0 - gp*gp);
    beta = alpha + sqrt(alpha * alpha - 1.0);
    y = log(beta)/ (double)ord;
    wp = 2.0 / (exp(y) + exp(-y));
    
    err = low2low(b, a, ord, wp, 1.0, b, a);

exit_label:
    return err;
}

