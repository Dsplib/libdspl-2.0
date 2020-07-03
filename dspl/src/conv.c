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
#include <string.h>
#include "dspl.h"



/*******************************************************************************
Real vectors linear convolution.
--------------------------------------------------------------------------------
Documented: RU, EN
*******************************************************************************/
int DSPL_API conv(double* a, int na, double* b, int nb, double* c)
{
  int k;
  int n;

  double *t;
  size_t bufsize;

  if(!a || !b || !c)
    return ERROR_PTR;
  if(na < 1 || nb < 1)
    return ERROR_SIZE;


  bufsize = (na + nb - 1) * sizeof(double);

  if((a != c) && (b != c))
    t = c;
  else
    t = (double*)malloc(bufsize);

  memset(t, 0, bufsize);

  for(k = 0; k < na; k++)
    for(n = 0; n < nb; n++)
      t[k+n] += a[k]*b[n];

  if(t!=c)
  {
    memcpy(c, t, bufsize);
    free(t);
  }
  return RES_OK;
}






/******************************************************************************
 Complex vectors linear convolution.
--------------------------------------------------------------------------------
Documented: RU, EN
*******************************************************************************/
int DSPL_API conv_cmplx(complex_t* a, int na, complex_t* b,
                        int nb, complex_t* c)
{
  int k;
  int n;

  complex_t *t;
  size_t bufsize;

  if(!a || !b || !c)
    return ERROR_PTR;
  if(na < 1 || nb < 1)
    return ERROR_SIZE;

  bufsize = (na + nb - 1) * sizeof(complex_t);

  if((a != c) && (b != c))
    t = c;
  else
    t = (complex_t*)malloc(bufsize);

  memset(t, 0, bufsize);

  for(k = 0; k < na; k++)
  {
    for(n = 0; n < nb; n++)
          {
      RE(t[k+n]) += CMRE(a[k], b[n]);
      IM(t[k+n]) += CMIM(a[k], b[n]);
    }
  }

  if(t!=c)
  {
    memcpy(c, t, bufsize);
    free(t);
  }

  return RES_OK;
}





/******************************************************************************
Real vectors fast linear convolution by using fast Fourier transform
--------------------------------------------------------------------------------
Documented: RU, EN
*******************************************************************************/
int DSPL_API conv_fft(double* a, int na,   double* b, int nb,
                      fft_t* pfft,  int nfft, double* c)
{
  complex_t *pa = NULL, *pb = NULL, *pc = NULL;
  int err;
  
  if(!a || !b || !c || !pfft)
    return ERROR_PTR;
  if(na<1 || nb < 1)
    return ERROR_SIZE;
  if(nfft<2)
    return ERROR_FFT_SIZE;
  
  pa = (complex_t*) malloc(na*sizeof(complex_t));
  pb = (complex_t*) malloc(nb*sizeof(complex_t));
  pc = (complex_t*) malloc((na+nb-1)*sizeof(complex_t));
  
  re2cmplx(a, na, pa);
  re2cmplx(b, nb, pb);
  
  err = conv_fft_cmplx(pa, na, pb, nb, pfft, nfft, pc);
  if(err != RES_OK)
    goto exit_label;
  
  err = cmplx2re(pc, na+nb-1, c, NULL);
  
exit_label:
  if(pa) free(pa);
  if(pb) free(pb);
  if(pc) free(pc);
  
  return err;
}




/******************************************************************************
Complex vectors fast linear convolution by using fast Fourier
transform algorithms
--------------------------------------------------------------------------------
Documented: RU, EN
*******************************************************************************/
int DSPL_API conv_fft_cmplx(complex_t* a, int na,   complex_t* b, int nb,
                            fft_t* pfft,  int nfft, complex_t* c)
{
  
  int La, Lb, Lc, Nz, n, p0, p1, ind, err;
  complex_t *pa, *pb;
  complex_t *pt, *pA, *pB, *pC;
  
  if(!a || !b || !c)
    return ERROR_PTR;
  if(na < 1 || nb < 1)
    return ERROR_SIZE;
  
  if(na >= nb)
  {
    La = na;
    Lb = nb;
    pa = a; 
    pb = b;
  }
  else
  {
    La = nb;
    pa = b;
    Lb = na;
    pb = a;
  }
    
  Lc = La + Lb - 1;
  Nz = nfft - Lb;

  if(Nz <= 0)
    return ERROR_FFT_SIZE;

  pt = (complex_t*)malloc(nfft*sizeof(complex_t));
  pB = (complex_t*)malloc(nfft*sizeof(complex_t));  
  pA = (complex_t*)malloc(nfft*sizeof(complex_t));  
  pC = (complex_t*)malloc(nfft*sizeof(complex_t));

  memset(pt,    0,  nfft*sizeof(complex_t));
  memcpy(pt+Nz, pb, Lb*sizeof(complex_t));

  err = fft_cmplx(pt, nfft, pfft, pB);
  if(err != RES_OK)
    goto exit_label;

  p0 = -Lb;
  p1 = p0 + nfft;
  ind = 0;
  while(ind < Lc)
  {
    if(p0 >=0)
    {
      if(p1 < La)
        err = fft_cmplx(pa + p0, nfft, pfft, pA);
      else
      {
        memset(pt, 0, nfft*sizeof(complex_t));
        memcpy(pt, pa+p0, (nfft+La-p1)*sizeof(complex_t));
        err = fft_cmplx(pt, nfft, pfft, pA);
      }
    }
    else
    {
      memset(pt, 0, nfft*sizeof(complex_t));
      if(p1 < La)        
        memcpy(pt - p0, pa, (nfft+p0)*sizeof(complex_t));
      else
        memcpy(pt - p0, pa, La * sizeof(complex_t));
      err = fft_cmplx(pt, nfft, pfft, pA);
    }
    
    if(err != RES_OK)
      goto exit_label;

    for(n = 0; n < nfft; n++)
    {
      RE(pC[n]) = CMRE(pA[n], pB[n]);
      IM(pC[n]) = CMIM(pA[n], pB[n]);
    }


    if(ind+nfft < Lc)
      err = ifft_cmplx(pC, nfft, pfft, c+ind);
    else
    {
      err = ifft_cmplx(pC, nfft, pfft, pt);
      memcpy(c+ind, pt, (Lc-ind)*sizeof(complex_t));
    }
    if(err != RES_OK)
      goto exit_label;
    
    p0  += Nz;
    p1  += Nz;
    ind += Nz;
  }
 
exit_label: 
  if(pt) free(pt);
  if(pB) free(pB);
  if(pA) free(pA);
  if(pC) free(pC);
  
  return err;
}



/*******************************************************************************
Real IIR filtration
--------------------------------------------------------------------------------
Documented: RU, EN
*******************************************************************************/
int DSPL_API filter_iir(double* b, double* a, int ord,
                        double* x, int n, double* y)
{
  double* buf = NULL;
  double* an  = NULL;
  double  u;
  int   k;
  int   m;
  int   count;

  if(!b || !x || !y)
    return  ERROR_PTR;

  if(ord < 1 || n < 1)
      return ERROR_SIZE;

  if(a && a[0]==0.0)
    return ERROR_FILTER_A0;

  count = ord + 1;
  buf = (double*) malloc(count*sizeof(double));
  an =  (double*) malloc(count*sizeof(double));

  memset(buf, 0, count*sizeof(double));

  if(!a)
    memset(an, 0, count*sizeof(double));
  else
    for(k = 0; k < count; k++)
      an[k] = a[k] / a[0];

  for(k = 0; k < n; k++)
  {
    for(m = ord; m > 0; m--)
      buf[m] = buf[m-1];
    u = 0.0;
    for(m = ord; m > 0; m--)
      u += buf[m]*an[m];

    buf[0] = x[k] - u;
    y[k] = 0.0;
    for(m = 0; m < count; m++)
      y[k] += buf[m] * b[m];
  }

  if(buf)
    free(buf);
  if(an)
    free(an);
  return RES_OK;
}

