/*
* Copyright (c) 2015-2019 Sergey Bakhurin
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
#include <stdio.h>
#include <string.h>
#include "dspl.h"
#include "dspl_internal.h"

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

   err = fft_krn(pfft->t1, pfft->t0, pfft, n, 0);

   if(err!=RES_OK)
    return err;

   norm = 1.0 / (double)n;
   for(k = 0; k < n; k++)
   {
     RE(y[k]) =  RE(pfft->t0[k])*norm;
     IM(y[k]) = -IM(pfft->t0[k])*norm;
   }
   return RES_OK;
}


/*******************************************************************************
Real vector FFT
*******************************************************************************/
int DSPL_API fft(double* x, int n, fft_t* pfft, complex_t* y)
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

  return fft_krn(pfft->t1, y, pfft, n, 0);
}




/*******************************************************************************
COMPLEX vector FFT
*******************************************************************************/
int DSPL_API fft_cmplx(complex_t* x, int n, fft_t* pfft, complex_t* y)
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

    return fft_krn(pfft->t1, y, pfft, n, 0);
}



/*******************************************************************************
composite FFT kernel
*******************************************************************************/
int fft_krn(complex_t* t0, complex_t* t1, fft_t* p, int n, int addr)
{
  int n1, n2, k, m, i;
  complex_t *pw = p->w+addr;
  complex_t tmp;
  
  n1 = 1;
  if(n%16== 0) { n1 = 16;   goto label_size;  }
  if(n%7 == 0) { n1 =  7;   goto label_size;  }
  if(n%8 == 0) { n1 =  8;   goto label_size;  }
  if(n%5 == 0) { n1 =  5;   goto label_size;  }
  if(n%4 == 0) { n1 =  4;   goto label_size;  }
  if(n%3 == 0) { n1 =  3;   goto label_size;  }
  if(n%2 == 0) { n1 =  2;   goto label_size;  }

label_size:
  if(n1 == 1)
  {
    for(k = 0; k < n; k++)
    {
      RE(t1[k]) = IM(t1[k]) = 0.0;
      for(m = 0; m < n; m++)
      {
        i = (k*m) % n;
        RE(tmp) = CMRE(t0[m], pw[i]);
        IM(tmp) = CMIM(t0[m], pw[i]);
        RE(t1[k]) += RE(tmp);
        IM(t1[k]) += IM(tmp);
      }
    }
  }
  else
  {
    n2 = n / n1;
    
    if(n2>1)
    {
      memcpy(t1, t0, n*sizeof(complex_t));
      matrix_transpose_cmplx(t1, n2, n1, t0);
    }
   
    if(n1 == 16)
      for(k = 0; k < n2; k++)
        dft16(t0+16*k, t1+16*k);
        
    if(n1 == 7)
      for(k = 0; k < n2; k++)
        dft7(t0+7*k, t1+7*k);
      
    if(n1 == 8)
      for(k = 0; k < n2; k++)
        dft8(t0+8*k, t1+8*k); 
        
    if(n1 == 5)
      for(k = 0; k < n2; k++)
        dft5(t0+5*k, t1+5*k);
   
    if(n1 == 4)
      for(k = 0; k < n2; k++)
        dft4(t0+4*k, t1+4*k);
    
    if(n1 == 3)
      for(k = 0; k < n2; k++)
        dft3(t0+3*k, t1+3*k);
    
    if(n1 == 2)
      for(k = 0; k < n2; k++)
        dft2(t0+2*k, t1+2*k);

    if(n2 > 1)
    {

      for(k =0; k < n; k++)
      {
        RE(t0[k]) = CMRE(t1[k], pw[k]);
        IM(t0[k]) = CMIM(t1[k], pw[k]);
      }

      matrix_transpose_cmplx(t0, n1, n2, t1);
      
      for(k = 0; k < n1; k++)
      {
        fft_krn(t1+k*n2, t0+k*n2, p, n2, addr+n);
      }
      matrix_transpose_cmplx(t0, n2, n1, t1);
    }
  }
  return RES_OK;

}




/*******************************************************************************
FFT create for composite N
*******************************************************************************/
int DSPL_API fft_create(fft_t* pfft, int n)
{

  int n1, n2, addr, s, k, m, nw, err;
  double phi;
  s = n;
  nw = addr = 0;

  if(pfft->n == n)
    return RES_OK;

  while(s > 1)
  {
    n2 = 1;
    if(s%16== 0) { n2 = 16; goto label_size; }
    if(s%7 == 0) { n2 =  7; goto label_size; }
    if(s%8 == 0) { n2 =  8; goto label_size; }
    if(s%5 == 0) { n2 =  5; goto label_size; }
    if(s%4 == 0) { n2 =  4; goto label_size; }
    if(s%3 == 0) { n2 =  3; goto label_size; }
    if(s%2 == 0) { n2 =  2; goto label_size; }


label_size:
    if(n2 == 1)
    {
      if(s > FFT_COMPOSITE_MAX)
      {
        err = ERROR_FFT_SIZE;
        goto error_proc;
      }
      
      nw += s;
      pfft->w = pfft->w ? (complex_t*) realloc(pfft->w,  nw*sizeof(complex_t)):
                          (complex_t*) malloc(           nw*sizeof(complex_t));
      for(k = 0; k < s; k++)
      {
        phi = - M_2PI * (double)k / (double)s;
        RE(pfft->w[addr]) = cos(phi);
        IM(pfft->w[addr]) = sin(phi);
        addr++;
      }
	  s = 1;
    }
    else
    {
      n1 = s / n2;
      nw += s;
      pfft->w = pfft->w ? (complex_t*) realloc(pfft->w,  nw*sizeof(complex_t)):
                          (complex_t*) malloc(           nw*sizeof(complex_t));

      for(k = 0; k < n1; k++)
      {
        for(m = 0; m < n2; m++)
        {
          phi = - M_2PI * (double)(k*m) / (double)s;
          RE(pfft->w[addr]) = cos(phi);
          IM(pfft->w[addr]) = sin(phi);
          addr++;
        }
      }
    }
    s /= n2;    
  }

  pfft->t0 = pfft->t0 ? (complex_t*) realloc(pfft->t0, n*sizeof(complex_t)):
                        (complex_t*) malloc(           n*sizeof(complex_t));

  pfft->t1 = pfft->t1 ? (complex_t*) realloc(pfft->t1, n*sizeof(complex_t)):
                        (complex_t*) malloc(           n*sizeof(complex_t));
                        
  pfft->n = n;

  return RES_OK;
error_proc:
  if(pfft->t0) free(pfft->t0);
  if(pfft->t1) free(pfft->t1);
  if(pfft->w)  free(pfft->w);
  pfft->n = 0;
  return err;
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
  memset(pfft, 0, sizeof(fft_t));
}



/*******************************************************************************
FFT magnitude for the real signal
*******************************************************************************/
int DSPL_API fft_mag(double* x, int n, fft_t* pfft, 
                     double fs, int flag,
                     double* mag, double* freq)
{
  int k, err = RES_OK;
  complex_t *X = NULL;
  
  if(!x || !pfft)
    return ERROR_PTR;
  
  if(n<1)
    return ERROR_SIZE;
  
  if(mag)
  {  
    X = (complex_t*)malloc(n*sizeof(complex_t));
    err = fft(x, n, pfft, X);
    if(err!=RES_OK)
      goto error_proc;
    
    if(flag & DSPL_FLAG_LOGMAG)
      for(k = 0; k < n; k++)
        mag[k] = 10.0*log10(ABSSQR(X[k]));
    else
      for(k = 0; k < n; k++)
        mag[k] = ABS(X[k]);
    if(flag & DSPL_FLAG_FFT_SHIFT)
    {
      err = fft_shift(mag, n, mag);
      if(err!=RES_OK)
        goto error_proc;
    }
  }
  
  if(freq)
  {
    if(flag & DSPL_FLAG_FFT_SHIFT)
      if(n%2)
        err = linspace(-fs*0.5 + fs*0.5/(double)n, 
                        fs*0.5 - fs*0.5/(double)n, 
                        n, DSPL_SYMMETRIC, freq);
      else
        err = linspace(-fs*0.5, fs*0.5, n, DSPL_PERIODIC, freq);
    else
      err = linspace(0, fs, n, DSPL_PERIODIC, freq);
  } 
       
error_proc:  
  if(X)
    free(X);
    
  return err;
}







/*******************************************************************************
FFT magnitude for the complex signal
*******************************************************************************/
int DSPL_API fft_mag_cmplx(complex_t* x, int n, fft_t* pfft, 
                           double fs, int flag,
                           double* mag, double* freq)
{
  int k, err = RES_OK;
  complex_t *X = NULL;
  
  if(!x || !pfft)
    return ERROR_PTR;
  
  if(n<1)
    return ERROR_SIZE;
  
  if(mag)
  {  
    X = (complex_t*)malloc(n*sizeof(complex_t));
    err = fft_cmplx(x, n, pfft, X);
    if(err!=RES_OK)
      goto error_proc;
    
    if(flag & DSPL_FLAG_LOGMAG)
      for(k = 0; k < n; k++)
        mag[k] = 10.0*log10(ABSSQR(X[k]));
    else
      for(k = 0; k < n; k++)
        mag[k] = ABS(X[k]);
    if(flag & DSPL_FLAG_FFT_SHIFT)
    {
      err = fft_shift(mag, n, mag);
      if(err!=RES_OK)
        goto error_proc;
    }
  }
  
  if(freq)
  {
    if(flag & DSPL_FLAG_FFT_SHIFT)
      if(n%2)
        err = linspace(-fs*0.5 + fs*0.5/(double)n, 
                        fs*0.5 - fs*0.5/(double)n, 
                        n, DSPL_SYMMETRIC, freq);
      else
        err = linspace(-fs*0.5, fs*0.5, n, DSPL_PERIODIC, freq);
    else
      err = linspace(0, fs, n, DSPL_PERIODIC, freq);
  }     
error_proc:  
  if(X)
    free(X);
    
  return err;
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
    n2 = (n+1) >> 1;
    buf = (double*) malloc(n2*sizeof(double));
    memcpy(buf, x, n2*sizeof(double));
    memcpy(y, x+n2, (n2-1)*sizeof(double));
    memcpy(y+n2-1, buf, n2*sizeof(double));
    free(buf);
  }
  return RES_OK;
}



/*******************************************************************************
FFT shifting for complex vector
*******************************************************************************/
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
    n2 = (n+1) >> 1;
    buf = (complex_t*) malloc(n2*sizeof(complex_t));
    memcpy(buf, x, n2*sizeof(complex_t));
    memcpy(y, x+n2, (n2-1)*sizeof(complex_t));
    memcpy(y+n2-1, buf, n2*sizeof(complex_t));
    free(buf);
  }
  return RES_OK;
}

