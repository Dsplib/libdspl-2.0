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
#include "fftw3.h"



/**************************************************************************************************
Create FFT object
***************************************************************************************************/
int DSPL_API fft_create( fft_t *pfft, int n)
{
    if(!pfft)
        return ERROR_PTR;
    if(n<1)
        return ERROR_SIZE;

    printf("1\n");
    if(pfft->in)
    {
        if(pfft->size != n)
            pfft->in = (complex_t*) realloc(pfft->in, n*sizeof(complex_t));
    }    
    else
        pfft->in = (complex_t*) malloc(n*sizeof(complex_t));
 
    if(pfft->out)
    {
        if(pfft->size != n)
            pfft->out = (complex_t*) realloc(pfft->out, n*sizeof(complex_t));
    }    
    else
        pfft->out = (complex_t*) malloc(n*sizeof(complex_t));
    
    pfft->size = n;
    printf("2\n");

    if(pfft->pfftw)
        fftw_destroy_plan(pfft->pfftw);

    printf("3\n");
    pfft->pfftw = (void*)fftw_plan_dft_1d(n, pfft->in, pfft->out, FFTW_FORWARD, FFTW_ESTIMATE);  
    printf("4\n");

    if(!pfft->pfftw)
        return ERROR_FFT_CREATE;

    return RES_OK;
}



/**************************************************************************************************
Destroy FFT object
***************************************************************************************************/
void DSPL_API fft_destroy(fft_t *pfft)
{
    if(pfft)
        return;

    if(pfft->pfftw)
        fftw_destroy_plan(pfft->pfftw);

    if(pfft->in)
        free(pfft->in);
    
    if(pfft->out)
        free(pfft->out);

    pfft->size = 0;
}



/**************************************************************************************************
Destroy FFT object
***************************************************************************************************/
int DSPL_API fft_cmplx(complex_t *x, int n, fft_t* pfft, complex_t* y)
{
    if(!x || !y || !pfft)
        return ERROR_PTR;

    if(n<1)
        return ERROR_SIZE;

    if(n != pfft->size)
       fft_create(pfft, n);

    memcpy(pfft->in, x, n*sizeof(complex_t));
    
    fftw_execute(pfft->pfftw);

    memcpy(y, pfft->out, n*sizeof(complex_t));
    return RES_OK;

}




/**************************************************************************************************
FFT shifting
***************************************************************************************************/
int DSPL_API fft_shift(double* x, int n, double* y)
{
	int n2, r;
	int k;
	double tmp;
	double *buf;
	
	if(!x || !y)
		return ERROR_PTR;

	if(n<1)
		return ERROR_SIZE;
		
	r = n%2;
	if(!r)
	{
		n2 = n>>1;
		for(k = 0; k < n2; k++)
		{
			tmp = x[k];
			y[k] = x[k+n2];
			y[k+n2] = tmp;
		}			
	}
	else
	{
		n2 = (n-1) >> 1;
		buf = (double*) malloc(n2*sizeof(double));
		memcpy(buf, x, n2*sizeof(double));
		memcpy(y, x+n2, (n2+1)*sizeof(double));
		memcpy(y+n2+1, buf, n2*sizeof(double));
		free(buf);
	}	
	return RES_OK;
}








