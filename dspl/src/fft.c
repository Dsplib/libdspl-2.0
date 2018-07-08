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

int fft_bit_reverse(complex_t* x, complex_t* y, int n, int p2);

int fft_dit(fft_t *pfft, int n, complex_t* y);

void fft_dit_krn(	complex_t *x0, complex_t *x1, complex_t *w, int n, 
			complex_t *y0, complex_t *y1);

int fft_p2(int n);





/*******************************************************************************
COMPLEX vector IFFT
*******************************************************************************/
int DSPL_API ifft_cmplx(complex_t *x, int n, fft_t* pfft, complex_t* y)
{
 	int err, k;
 	double norm;
 	
 	if(!x || !pfft || !y)
		return ERROR_PTR;
 	if(n<1)
		return ERROR_SIZE;
 	
 	
 	err = fft_create(pfft, n);
 	if(err != RES_OK)
		return err;
 	
 	memcpy(pfft->t1, x, n*sizeof(complex_t));
 	for(k = 0; k < n; k++)
		IM(pfft->t1[k]) = -IM(pfft->t1[k]);
 	
 	err = fft_dit(pfft, n, y);
 	
 	if(err!=RES_OK)
		return err;
 	
 	norm = 1.0 / (double)n;
 	for(k = 0; k < n; k++)
 	{
 		RE(y[k]) =  RE(y[k])*norm;
 		IM(y[k]) = -IM(y[k])*norm;
 	}
 	return RES_OK; 
}







/*******************************************************************************
Real vector FFT
*******************************************************************************/
int DSPL_API fft(double *x, int n, fft_t* pfft, complex_t* y)
{
	int err;
	
	if(!x || !pfft || !y)
		return ERROR_PTR;
	if(n<1)
		return ERROR_SIZE;
	
	
	err = fft_create(pfft, n);
	if(err != RES_OK)
		return err;
	
	re2cmplx(x, n, pfft->t1);
	
	return fft_dit(pfft, n, y);
}




/*******************************************************************************
COMPLEX vector FFT
*******************************************************************************/
int DSPL_API fft_cmplx(complex_t *x, int n, fft_t* pfft, complex_t* y)
{
    int err;

    if(!x || !pfft || !y)
        return ERROR_PTR;
    if(n<1)
        return ERROR_SIZE;


    err = fft_create(pfft, n);
    if(err != RES_OK)
        return err;
    
    memcpy(pfft->t1, x, n*sizeof(complex_t));

    return fft_dit(pfft, n, y);
}



/*******************************************************************************
FFT bit reverse
*******************************************************************************/
int fft_bit_reverse(complex_t* x, complex_t* y, int n, int p2)
{
	static unsigned char rb_table[256] = 
	{
		0x00, 0x80, 0x40, 0xC0, 0x20, 0xA0, 0x60, 0xE0, 
		0x10, 0x90, 0x50, 0xD0, 0x30, 0xB0, 0x70, 0xF0, 
		0x08, 0x88, 0x48, 0xC8, 0x28, 0xA8, 0x68, 0xE8,
		0x18, 0x98, 0x58, 0xD8, 0x38, 0xB8, 0x78, 0xF8, 
		0x04, 0x84, 0x44, 0xC4, 0x24, 0xA4, 0x64, 0xE4, 
		0x14, 0x94, 0x54, 0xD4, 0x34, 0xB4, 0x74, 0xF4, 
		0x0C, 0x8C, 0x4C, 0xCC, 0x2C, 0xAC, 0x6C, 0xEC, 
		0x1C, 0x9C, 0x5C, 0xDC, 0x3C, 0xBC, 0x7C, 0xFC, 
		0x02, 0x82, 0x42, 0xC2, 0x22, 0xA2, 0x62, 0xE2, 
		0x12, 0x92, 0x52, 0xD2, 0x32, 0xB2, 0x72, 0xF2, 
		0x0A, 0x8A, 0x4A, 0xCA, 0x2A, 0xAA, 0x6A, 0xEA, 
		0x1A, 0x9A, 0x5A, 0xDA, 0x3A, 0xBA, 0x7A, 0xFA,
		0x06, 0x86, 0x46, 0xC6, 0x26, 0xA6, 0x66, 0xE6, 
		0x16, 0x96, 0x56, 0xD6, 0x36, 0xB6, 0x76, 0xF6, 
		0x0E, 0x8E, 0x4E, 0xCE, 0x2E, 0xAE, 0x6E, 0xEE, 
		0x1E, 0x9E, 0x5E, 0xDE, 0x3E, 0xBE, 0x7E, 0xFE,
		0x01, 0x81, 0x41, 0xC1, 0x21, 0xA1, 0x61, 0xE1, 
		0x11, 0x91, 0x51, 0xD1, 0x31, 0xB1, 0x71, 0xF1,
		0x09, 0x89, 0x49, 0xC9, 0x29, 0xA9, 0x69, 0xE9, 
		0x19, 0x99, 0x59, 0xD9, 0x39, 0xB9, 0x79, 0xF9, 
		0x05, 0x85, 0x45, 0xC5, 0x25, 0xA5, 0x65, 0xE5,
		0x15, 0x95, 0x55, 0xD5, 0x35, 0xB5, 0x75, 0xF5,
		0x0D, 0x8D, 0x4D, 0xCD, 0x2D, 0xAD, 0x6D, 0xED, 
		0x1D, 0x9D, 0x5D, 0xDD, 0x3D, 0xBD, 0x7D, 0xFD,
		0x03, 0x83, 0x43, 0xC3, 0x23, 0xA3, 0x63, 0xE3, 
		0x13, 0x93, 0x53, 0xD3, 0x33, 0xB3, 0x73, 0xF3, 
		0x0B, 0x8B, 0x4B, 0xCB, 0x2B, 0xAB, 0x6B, 0xEB, 
		0x1B, 0x9B, 0x5B, 0xDB, 0x3B, 0xBB, 0x7B, 0xFB,
		0x07, 0x87, 0x47, 0xC7, 0x27, 0xA7, 0x67, 0xE7, 
		0x17, 0x97, 0x57, 0xD7, 0x37, 0xB7, 0x77, 0xF7, 
		0x0F, 0x8F, 0x4F, 0xCF, 0x2F, 0xAF, 0x6F, 0xEF, 
		0x1F, 0x9F, 0x5F, 0xDF, 0x3F, 0xBF, 0x7F, 0xFF
	};

	if(!x || !y)
		return ERROR_PTR;
	if(n<1 || p2 < 1)
		return ERROR_SIZE;
	
	unsigned int v, c;
	
	for(v = 0; v < n; v++)
	{
		c = (unsigned int)
			((rb_table[ v        & 0xff] << 24) |
			(rb_table[(v >>  8) & 0xff] << 16) | 
			(rb_table[(v >> 16) & 0xff] <<  8) | 
			(rb_table[(v >> 24) & 0xff])) >> (32 - p2);
			
		RE(y[c]) = RE(x[v]);
		IM(y[c]) = IM(x[v]);
	} 
	return RES_OK;
}





/*******************************************************************************
FFT create
*******************************************************************************/
int DSPL_API fft_create(fft_t *pfft, int n)
{
    int p2, k,r,m,addr,s;
    double phi;

    p2 = fft_p2(n);
    if(p2 < 1)
        return ERROR_FFT_SIZE;
    if(n < pfft->n+1)
        return RES_OK;

    pfft->n = n;
    pfft->p2 = p2;

    pfft->w = pfft->w   ? (complex_t*) realloc(pfft->w,  n*sizeof(complex_t)):
                          (complex_t*) malloc(           n*sizeof(complex_t));


    pfft->t0 = pfft->t0 ? (complex_t*) realloc(pfft->t0, n*sizeof(complex_t)):
                          (complex_t*) malloc(           n*sizeof(complex_t));


    pfft->t1 = pfft->t1 ? (complex_t*) realloc(pfft->t1, n*sizeof(complex_t)):
                          (complex_t*) malloc(           n*sizeof(complex_t));

    m = 0;
    addr = 0;

    for(k = 0; k < p2; k++)
    {
        s = 1<<m;
        for( r = 0; r < s; r++)
        {
            phi = -M_2PI *(double)r / (double)(2*s);
            RE(pfft->w[addr+r]) = cos(phi);
            IM(pfft->w[addr+r]) = sin(phi);
        }
        addr+=s;
        m++;
    }
    return RES_OK;
}



/*******************************************************************************
FFT decimation in time
*******************************************************************************/
int fft_dit(fft_t *pfft, int n, complex_t* y)
{
    int k,s,m,waddr, dm,p2;
    complex_t *t, *t0, *t1, *w;
    int err;

    p2 = fft_p2(n);
    if(p2<0)
        return ERROR_FFT_SIZE;

    t0 = pfft->t0;
    t1 = pfft->t1;
    w  = pfft->w;

    s = n>>1;
    m = 1;
    waddr = 0;
 
    err = fft_bit_reverse(t1, t0, n, p2);
    if(err!= RES_OK)
        return err;
   
    while(s)
    {
        dm = m<<1;
        if(s > 1)
        {
            for(k = 0; k < n; k+=dm)
            {
                fft_dit_krn(t0+k, t0+k+m, w+waddr, m, t1+k, t1+k+m);
            }
            t = t1;
            t1 = t0;
            t0 = t;
            waddr+=m;
            m <<= 1;
        }
        else
        {
            fft_dit_krn(t0, t0+m, w+waddr, m, y, y+m);
        }
        s >>= 1;     
    }
    

    return RES_OK;
}





/*******************************************************************************
FFT decimation in time kernel
*******************************************************************************/
void fft_dit_krn(complex_t *x0, complex_t *x1, complex_t *w, int n, 
				 complex_t *y0, complex_t *y1)
{
    int k;
    complex_t mul;
    for(k = 0; k < n; k++)
    {
        RE(mul) = CMRE(x1[k], w[k]);
        IM(mul) = CMIM(x1[k], w[k]);

        RE(y0[k]) = RE(x0[k]) + RE(mul);
        IM(y0[k]) = IM(x0[k]) + IM(mul);
        
        RE(y1[k]) = RE(x0[k]) - RE(mul);
        IM(y1[k]) = IM(x0[k]) - IM(mul);
    } 
    
}





/*******************************************************************************
FFT free
*******************************************************************************/
void DSPL_API fft_free(fft_t *pfft)
{
    if(!pfft)
        return;
    if(pfft->w)
        free(pfft->w);
    if(pfft->t0)
        free(pfft->t0);
    if(pfft->t1)
        free(pfft->t1);
}





/*******************************************************************************
FFT power of 2
*******************************************************************************/
int fft_p2(int n)
{
    int p = (n>0) ? n : 0;
    int p2 = 0;
    while((p = p>>1))
        p2++;
    if((1<<p2)!=n)
        return -1;
    return p2;
}







/*******************************************************************************
FFT shifting
*******************************************************************************/
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



/**************************************************************************************************
FFT shifting for complex vector
***************************************************************************************************/
int DSPL_API fft_shift_cmplx(complex_t* x, int n, complex_t* y)
{
	int n2, r;
	int k;
	complex_t tmp;
	complex_t *buf;
	
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
			RE(tmp) = RE(x[k]);
			IM(tmp) = IM(x[k]);
                        
                        RE(y[k]) = RE(x[k+n2]);
			IM(y[k]) = IM(x[k+n2]);
                        
                        RE(y[k+n2]) = RE(tmp);
                        IM(y[k+n2]) = IM(tmp);
		}			
	}
	else
	{
		n2 = (n-1) >> 1;
		buf = (complex_t*) malloc(n2*sizeof(complex_t));
		memcpy(buf, x, n2*sizeof(complex_t));
		memcpy(y, x+n2, (n2+1)*sizeof(complex_t));
		memcpy(y+n2+1, buf, n2*sizeof(complex_t));
		free(buf);
	}	
	return RES_OK;
}








