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


#include <stdlib.h>
#include <string.h>
#include "dspl.h"


/**************************************************************************************************
Composition of rational functions
***************************************************************************************************/
int rcompos(double *b, double *a, int n, double *c, double *d, int p, double *beta, double *alpha)
{	
	int k2, i, k,  pn, pd, ln, ld, k2s, nk2s;
	double *num = NULL, *den = NULL, *ndn = NULL, *ndd = NULL;
	int res;
	
	if (!a || !b || !c || !d || !beta || !alpha)
		return ERROR_PTR;
		
	if(n < 0 || p < 0)
		return ERROR_SIZE;


	k2   = (n*p)+1; 
	k2s  = k2*sizeof(double);	    /* alpha and beta size			*/
	nk2s = (n+1)*k2*sizeof(double); /* num, den, ndn and ndd size	*/

	num = (double*)malloc(nk2s);
	den = (double*)malloc(nk2s);
	ndn = (double*)malloc(nk2s);
	ndd = (double*)malloc(nk2s);

	memset(num,   0, nk2s);
	memset(den,   0, nk2s);
	memset(ndn,   0, nk2s);
	memset(ndd,   0, nk2s);

	                                
	num[0] = den[0] = 1.0;
	pn = 0;
	ln = 1;
	for(i = 1; i < n+1; i++)
	{
		res = conv(num+pn, ln, c, p+1, num+pn+k2);
		if(res!=RES_OK)
			goto exit_label;
		res = conv(den+pn, ln, d, p+1, den+pn+k2);
		if(res!=RES_OK)
			goto exit_label;
		pn += k2;
		ln += p;
	}

	pn = 0; 
	pd = n*k2;
	ln = 1;
	ld = k2;

	for (i = 0; i < n+1; i++)
	{
		res = conv(num + pn, ln, den + pd, ld, ndn + i*k2);
		if(res!=RES_OK)
			goto exit_label;
		ln += p;
		ld -= p;
		pn += k2;
		pd -= k2;
	}

	for (i = 0; i < n+1; i++)
	{
		for (k = 0; k < k2; k++)
		{
			ndd[i*k2 + k]  = ndn[i*k2 + k] * a[i];
			ndn[i*k2 + k] *= b[i];
		}
	}

	memset(alpha, 0, k2s);
	memset(beta,  0, k2s);

	for (k = 0; k < k2; k++)
	{
		for (i = 0; i < n+1; i++)
		{
			beta[k]  += ndn[i*k2 + k];
			alpha[k] += ndd[i*k2 + k];
		}
	}
	
	res = RES_OK;
exit_label:
	if(num)
		free(num);
	if(den)
		free(den);
	if(ndn)
		free(ndn);
	if(ndd)
		free(ndd);

	return res;
}
