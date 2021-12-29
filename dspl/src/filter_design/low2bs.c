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
int DSPL_API low2bs(double* b, double* a, int ord,
                    double w0, double wsl, double wsh,
                    double* beta, double* alpha)
{

    double den[3] = {0.0, 0.0, 1.0};
    double num[3] = {0.0, 0.0, 0.0};

    if(!b || !a || !beta || !alpha)
        return ERROR_PTR;
    if(ord < 1)
        return ERROR_FILTER_ORD;
    if(w0 <= 0.0 || wsl <= 0.0 || wsh <= 0.0 || wsh <= wsl)
        return ERROR_FILTER_FT;

    den[0] = (wsh * wsl) / (w0 * w0);
    num[1] = (wsh - wsl) / w0;

    return ratcompos(b, a, ord, num, den, 2, beta, alpha);
}

