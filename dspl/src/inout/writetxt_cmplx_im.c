/*
* Copyright (c) 2015-2024 Sergey Bakhurin
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
#include "dspl.h"






#ifdef DOXYGEN_ENGLISH

#endif
#ifdef DOXYGEN_RUSSIAN

#endif
int DSPL_API writetxt_cmplx_im(double* x, complex_t *y, int n, char* fn)
{
    int k;
    FILE* pFile = NULL;

    if(!x)
        return ERROR_PTR;
    if(n < 1)
        return ERROR_SIZE;
    if(!fn)
        return ERROR_FNAME;

    pFile = fopen(fn, "w+");
    if(pFile == NULL)
        return ERROR_FOPEN;

    if(x)
        for(k = 0; k < n; k++)
            fprintf(pFile, "%+.12E\t%+.12E\n", x[k], IM(y[k]));
    else
        for(k = 0; k < n; k++)
            fprintf(pFile, "%+.12E\n", IM(y[k]));

    fclose(pFile);
    return RES_OK;
}

