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
Goertzel algorithm for real vector
**************************************************************************************************/
int DSPL_API goertzel(double *x, int n, int *ind, int k, complex_t *y)
{
	
    int m, p;
    double wR, wI;
    double alpha;
    double v[3];

    if(!x || !y || !ind)
		return ERROR_PTR;

    if(n < 1 || k < 1)
		return ERROR_SIZE;

	for(p = 0; p < k; p++)
	{
		wR = cos(M_2PI * (double)ind[p] / (double)n);
		wI = sin(M_2PI * (double)ind[p] / (double)n);
		
		alpha = 2.0 * wR;
		v[0] = v[1] = v[2] = 0.0;
		
		for(m = 0; m < n; m++)
		{
			v[2] = v[1];
			v[1] = v[0];
			v[0] = x[m]+alpha*v[1] - v[2];
		}
		RE(y[p]) = wR * v[0] - v[1];
		IM(y[p]) = wI * v[0];
	}	
    
	return RES_OK;		
}





/**************************************************************************************************
Goertzel algorithm for complex vector
**************************************************************************************************/
int DSPL_API goertzel_cmplx(complex_t *x, int n, int *ind, int k, complex_t *y)
{
	
    int m, p;
    complex_t w;
    double alpha;
    complex_t v[3];
	
    if(!x || !y || !ind)
		return ERROR_PTR;

    if(n < 1 || k < 1)
		return ERROR_SIZE;

	for(p = 0; p < k; p++)
	{
		RE(w) = cos(M_2PI * (double)ind[p] / (double)n);
		IM(w) = sin(M_2PI * (double)ind[p] / (double)n);
		
		alpha = 2.0 * RE(w);
		memset(v, 0, 3*sizeof(complex_t)); 
		
		for(m = 0; m < n; m++)
		{
			RE(v[2]) = RE(v[1]);
			RE(v[1]) = RE(v[0]);
			RE(v[0]) = RE(x[m]) + alpha * RE(v[1]) - RE(v[2]);

            IM(v[2]) = IM(v[1]);
			IM(v[1]) = IM(v[0]);
			IM(v[0]) = IM(x[m]) + alpha * IM(v[1]) - IM(v[2]);
		}

        RE(y[p]) = CMRE(w, v[0]) - RE(v[1]);
		IM(y[p]) = CMIM(w, v[0]) - IM(v[1]);
	}	
    
	return RES_OK;		
}
