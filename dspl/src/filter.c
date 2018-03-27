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
#include "dspl.h"





/**************************************************************************************************
Zeros and poles to filter coefficients recalc
***************************************************************************************************/
int DSPL_API filter_zp2ab(complex_t *z, int nz, complex_t *p, int np, int ord, double* b, double* a)
{
    complex_t *acc = NULL;
    int res;

    if(!z || !p || !b || !a)
        return ERROR_PTR;
    if(nz < 0 || np < 0)
        return ERROR_SIZE;
    if(nz < ord || np < ord)
        return ERROR_POLY_ORD;

    acc = (complex_t*) malloc((ord+1) * sizeof(complex_t));
    res = poly_z2a_cmplx(z, nz, ord, acc);
    if(res != RES_OK)
        goto exit_label;    
    res = cmplx2re(acc, ord+1, b, NULL);
    if(res != RES_OK)
        goto exit_label;

    res = poly_z2a_cmplx(p, np, ord, acc);
    if(res != RES_OK)
        goto exit_label;    
    res = cmplx2re(acc, ord+1, a, NULL);
    if(res != RES_OK)
        goto exit_label;

exit_label:
    if(acc)
        free(acc);
    return res;
    
}






