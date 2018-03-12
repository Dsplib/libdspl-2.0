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
#include <math.h>
#include "dspl.h"



/**************************************************************************************************
Chebyshev polynomials of the first kind
***************************************************************************************************/
int DSPL_API cheby_poly1(double* x, int n, int ord, double* y)
{
    int k, m;
    double t[2];

    if(!x || !y)
        return ERROR_PTR;
    if(n < 1)
        return ERROR_SIZE;
    if(ord<0)
        return ERROR_POLY_ORD;
    if(ord==0)
    {
        for(k = 0; k < n; k++)
        {
            y[k] = 1.0;             
        }
        return RES_OK;
    }

    if(ord==1)
    {
        memcpy(y, x, n*sizeof(double));
        return RES_OK;
    }

    for(k = 0; k < n; k++)
    {
        m = 2;
        t[1]  = x[k];
        t[0]  = 1.0;  
        while(m <= ord)
        {   
            y[k] = 2.0 * x[k] *t[1] - t[0];
            t[0] = t[1];
            t[1] = y[k]; 
            m++;           
        }         
    }
    return RES_OK;
}




/**************************************************************************************************
Chebyshev polynomials of the second kind
***************************************************************************************************/
int DSPL_API cheby_poly2(double* x, int n, int ord, double* y)
{
    int k, m;
    double t[2];

    if(!x || !y)
        return ERROR_PTR;
    if(n < 1)
        return ERROR_SIZE;
    if(ord<0)
        return ERROR_POLY_ORD;
    if(ord==0)
    {
        for(k = 0; k < n; k++)
        {
            y[k] = 1.0;             
        }
        return RES_OK;
    }

    if(ord==1)
    {
        for(k = 0; k < n; k++)
        {
            y[k] = 2.0*x[n];             
        };
        return RES_OK;
    }

    for(k = 0; k < n; k++)
    {
        m = 2;
        t[1]  = 2.0*x[n];
        t[0]  = 1.0;  
        while(m <= ord)
        {   
            y[k] = 2.0 * x[k] *t[1] - t[0];
            t[0] = t[1];
            t[1] = y[k];
            m++;            
        }         
    }
    return RES_OK;



}
