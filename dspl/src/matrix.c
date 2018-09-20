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


void transpose(double* a, int n, int m, double* b);
void transpose_cmplx(complex_t* a, int n, int m, complex_t* b);
void transpose_hermite(complex_t* a, int n, int m, complex_t* b);




/*******************************************************************************
Real matrx transpose
*******************************************************************************/
int DSPL_API matrix_create(matrix_t* a, int n, int m, int type)
{
	if(!a)
		return ERROR_PTR;
	if(n < 1 || m < 1)
		return ERROR_MATRIX_SIZE;

	if(a->dat)
	{
		a->dat = (type & DAT_MASK) ?
			(void*) realloc(a->dat, n*m*sizeof(complex_t)):
			(void*) realloc(a->dat, n*m*sizeof(double));
	}
	else
	{
		a->dat = (type & DAT_MASK) ?
			(void*) malloc(n*m*sizeof(complex_t)):
			(void*) malloc(n*m*sizeof(double));
	}

	a->n = n;
	a->m = m;
	a->type = type;

	return RES_OK;
}











/*******************************************************************************
Real matrx transpose
*******************************************************************************/
void transpose(double* a, int n, int m, double* b)
{
	int p, q, i, j, aind, bind;

	for(p = 0; p < n - DSPL_MATRIX_BLOCK; p+=DSPL_MATRIX_BLOCK)
	{
		for(q = 0; q < m - DSPL_MATRIX_BLOCK; q+=DSPL_MATRIX_BLOCK)
		{
			for(i = 0; i < DSPL_MATRIX_BLOCK; i++)
			{
				for(j = 0; j < DSPL_MATRIX_BLOCK; j++)
				{
					aind = (q+j) * n + p + i;
					bind = (p+i) * m + q + j;
					b[bind] = a[aind];
				}
			}
		}
	}
	for(i = p; i < n; i++)
		for(j = 0; j < m; j++)
			b[i*m + j] = a[j*n+i];

	for(i = 0; i < p; i++)
		for(j = q; j < m; j++)
			b[i*m + j] = a[j*n+i];

}




/*******************************************************************************
Complex matrx transpose
*******************************************************************************/
void transpose_cmplx(complex_t* a, int n, int m, complex_t* b)
{
	int p, q, i, j, aind, bind;

	for(p = 0; p < n - DSPL_MATRIX_BLOCK; p+=DSPL_MATRIX_BLOCK)
	{
		for(q = 0; q < m - DSPL_MATRIX_BLOCK; q+=DSPL_MATRIX_BLOCK)
		{
			for(i = 0; i < DSPL_MATRIX_BLOCK; i++)
			{
				for(j = 0; j < DSPL_MATRIX_BLOCK; j++)
				{
					aind = (q+j) * n + p + i;
					bind = (p+i) * m + q + j;
					RE(b[bind]) = RE(a[aind]);
					IM(b[bind]) = IM(a[aind]);
				}
			}
		}
	}
	for(i = p; i < n; i++)
	{
		for(j = 0; j < m; j++)
		{
			RE(b[i*m + j]) = RE(a[j*n+i]);
			IM(b[i*m + j]) = IM(a[j*n+i]);
		}
	}

	for(i = 0; i < p; i++)
	{
		for(j = q; j < m; j++)
		{
			RE(b[i*m + j]) = RE(a[j*n+i]);
			IM(b[i*m + j]) = IM(a[j*n+i]);
		}
	}
}




/*******************************************************************************
Hermite matrx transpose
*******************************************************************************/
void transpose_hermite(complex_t* a, int n, int m, complex_t* b)
{
	int p, q, i, j, aind, bind;

	for(p = 0; p < n - DSPL_MATRIX_BLOCK; p+=DSPL_MATRIX_BLOCK)
	{
		for(q = 0; q < m - DSPL_MATRIX_BLOCK; q+=DSPL_MATRIX_BLOCK)
		{
			for(i = 0; i < DSPL_MATRIX_BLOCK; i++)
			{
				for(j = 0; j < DSPL_MATRIX_BLOCK; j++)
				{
					aind = (q+j) * n + p + i;
					bind = (p+i) * m + q + j;
					RE(b[bind]) =  RE(a[aind]);
					IM(b[bind]) = -IM(a[aind]);
				}
			}
		}
	}
	for(i = p; i < n; i++)
	{
		for(j = 0; j < m; j++)
		{
			RE(b[i*m + j]) = RE(a[j*n+i]);
			IM(b[i*m + j]) = IM(a[j*n+i]);
		}
	}

	for(i = 0; i < p; i++)
	{
		for(j = q; j < m; j++)
		{
			RE(b[i*m + j]) =  RE(a[j*n+i]);
			IM(b[i*m + j]) = -IM(a[j*n+i]);
		}
	}
}

