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
int DSPL_API verif_data_gen(int len, int type, char* fn)
{
    double    *pd = NULL;
    complex_t *pc = NULL;
    random_t rnd = {0};
    int err;
    if(len < 1)
        return ERROR_SIZE;
    if(!fn)
        return ERROR_FNAME;
    
    err = random_init(&rnd, RAND_TYPE_MRG32K3A, NULL);
    if(err != RES_OK)
        goto exit_label;
    switch(type & DAT_MASK)
    {
        case DAT_DOUBLE:
            pd = (double*)malloc(len*sizeof(double));
            if(!pd)
            {
                err = ERROR_MALLOC;
                goto exit_label;
            }
            err  = randn(pd, len, 1.0, 10.0, &rnd);
            if(err != RES_OK)
                goto exit_label;
            
            err = writebin(pd, len, 1, type, fn);
            if(err != RES_OK)
                goto exit_label;
            break;
        case DAT_COMPLEX:
            pc = (complex_t*)malloc(len*sizeof(complex_t));
            if(!pc)
            {
                err = ERROR_MALLOC;
                goto exit_label;
            }
            err  = randn((double*)pc, 2*len, 1.0, 10.0, &rnd);
            if(err != RES_OK)
                goto exit_label;
            
            err = writebin(pc, len, 1, type, fn);
            if(err != RES_OK)
                goto exit_label;
            break;
        default:
            err = ERROR_DAT_TYPE;
    }
    
exit_label:
    if(pd)
        free(pd);
    if(pc)
        free(pc);
    return err;
}


