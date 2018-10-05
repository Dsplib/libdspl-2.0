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


#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "dspl.h"


/**************************************************************************************************
Concntenate arrays
***************************************************************************************************/
int DSPL_API concat(void* a, size_t na, void *b, size_t nb, void* c)
{
        if(!a || !b || !c || c == b)
                return ERROR_PTR;
        if(na < 1 || nb < 1)
                return ERROR_SIZE;
        
    if(c != a)
        memcpy(c, a, na);

    memcpy((char*)c+na, b, nb);
        return RES_OK;	
}






/**************************************************************************************************
Flip real array in place 
***************************************************************************************************/
int DSPL_API flipip(double* x, int n)
{
    int k;
    double tmp;
    if(!x)
        return ERROR_PTR;
    if(n<1) 
        return ERROR_SIZE;
    
    for(k = 0; k < n/2; k++)
    {
        tmp = x[k];
        x[k] = x[n-1-k];
        x[n-1-k] = tmp;
    }
    return RES_OK;
}



/**************************************************************************************************
Flip complex array in place 
***************************************************************************************************/
int DSPL_API flipip_cmplx(complex_t* x, int n)
{
    int k;
    complex_t tmp;
    if(!x)
        return ERROR_PTR;
    if(n<1) 
        return ERROR_SIZE;
    
    for(k = 0; k < n/2; k++)
    {
        RE(tmp) = RE(x[k]);
        RE(x[k]) = RE(x[n-1-k]);
        RE(x[n-1-k]) = RE(tmp);

        IM(tmp) = IM(x[k]);
        IM(x[k]) = IM(x[n-1-k]);
        IM(x[n-1-k]) = IM(tmp);
    }
    return RES_OK;
}





