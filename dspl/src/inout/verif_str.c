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
void DSPL_API verif_str(double* yout, int nout, 
                        char* str_msg, char* outfn, char* logfn)
{
    char str[VERIF_STR_BUF] = {0};
    char msg[VERIF_STR_BUF] = {0};
    double *y = NULL;
    double derr = 0.0;
    int n, m, verr, type;

    sprintf(str, "%s", str_msg);
    while(strlen(str) < VERIF_STR_LEN)
        str[strlen(str)] = VERIF_CHAR_POINT;

    readbin(outfn, (void**)(&y), &n, &m, &type);

    if(nout != n*m)
    {
        sprintf(msg, "FAILED (out size [%d] neq [%d])", n, nout);
        strcat(str, msg);
        addlog(str, logfn);
        printf("%s\n", str);
        return;
    }
    
    if(type!=DAT_DOUBLE)
    {
        sprintf(msg, "FAILED (type is complex)");
        strcat(str, msg);
        addlog(str, logfn);
        printf("%s\n", str);
        return;
    }

    verr = verif(yout, y, nout, VERIF_LEVEL_DOUBLE, &derr);
    if(verr == DSPL_VERIF_SUCCESS)
        sprintf(msg, "ok (err = %12.4E)", derr);
    else
        sprintf(msg, "FAILED (err = %12.4E)", derr);
    strcat(str, msg);
    addlog(str, logfn);
    printf("%s\n", str);
}

