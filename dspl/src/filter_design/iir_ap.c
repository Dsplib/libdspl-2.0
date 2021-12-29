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
#include "dspl_internal.h"



/******************************************************************************
Analog prototype for IIR
*******************************************************************************/
int iir_ap(double rp, double rs, int ord, int type, double* b, double* a)
{
    int err;
    switch(type & DSPL_FILTER_APPROX_MASK)
    {
        case DSPL_FILTER_BUTTER:
            err = butter_ap(rp, ord, b, a);
            break;
        case DSPL_FILTER_CHEBY1:
            err = cheby1_ap(rp, ord, b, a);
            break;
        case DSPL_FILTER_CHEBY2:
            err = cheby2_ap_wp1(rp, rs, ord, b, a);
            break;
        case DSPL_FILTER_ELLIP:
            err = ellip_ap(rp, rs, ord, b, a);
            break;
        default:
            err = ERROR_FILTER_APPROX;
    }
    return err;
}




