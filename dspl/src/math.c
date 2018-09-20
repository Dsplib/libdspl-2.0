/*
* Copyright (c) 2015-2017 Sergey Bakhurin
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
* MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
* GNU General Public License for more details.
*
* You should have received a copy of the GNU General Public License
* along with Foobar.  If not, see <http://www.gnu.org/licenses/>.
*/


#include <stdio.h>
#include <stdlib.h>
#include "dspl.h"


/*******************************************************************************
module operator for double
*******************************************************************************/
double DSPL_API dmod (double x, double y)
{
	if(y == 0.0)
		return x;
	return  x - floor(x/y) * y;
}





/*******************************************************************************
sinc(x) = sin(pi*x)/(pi*x)
*******************************************************************************/
int DSPL_API sinc(double* x, int n, double a, double* y)
{
	int k;

	if(!x || !y)
		return ERROR_PTR;
	if(n<1)
		return ERROR_SIZE;

	for(k = 0; k < n; k++)
		y[k] = (x[k]==0.0) ? 1.0 : sin(a*x[k])/(a*x[k]);

	return RES_OK;


}

