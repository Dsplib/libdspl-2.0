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
Real vectors linear convolution
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
Complex vectors linear convolution
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






/*******************************************************************************
 Complex vectors FFT linear convolution
 ******************************************************************************/
int DSPL_API conv_fft_cmplx(complex_t* a, int na, complex_t* b, int nb,
                            fft_t* pfft, complex_t* c)
{
  complex_t *pa = NULL;
  complex_t *pb = NULL;
  complex_t *pc = NULL;
  complex_t *pA = NULL;
  complex_t *pB = NULL;
  complex_t *pC = NULL;

  int nfft, nfft2, n, npos, err;
  int ma, mb;
  complex_t *ta, *tb;

  if(!a || !b || !c)
    return ERROR_PTR;
  if(na < 1 || nb < 1)
    return ERROR_SIZE;


  if(na > nb)
  {
    ma = na;
    mb = nb;
    ta = a;
    tb = b;
  }
  else
  {
    ma = nb;
    mb = na;
    ta = b;
    tb = a;
  }
  if(ma > 2*mb)
  {
    nfft = 4;
    n = mb-1;
    while(n>>=1)
      nfft <<= 1;
    nfft2 = nfft >> 1;

    pa = (complex_t*)malloc(nfft * sizeof(complex_t));
    pb = (complex_t*)malloc(nfft * sizeof(complex_t));
    pc = (complex_t*)malloc(nfft * sizeof(complex_t));
    pA = (complex_t*)malloc(nfft * sizeof(complex_t));
    pB = (complex_t*)malloc(nfft * sizeof(complex_t));
    pC = (complex_t*)malloc(nfft * sizeof(complex_t));

    npos = -nfft2;
    memset(pa, 0, nfft*sizeof(complex_t));
    memset(pb, 0, nfft*sizeof(complex_t));

    memcpy(pa + nfft2, ta, nfft2 * sizeof(complex_t));
    memcpy(pb,     tb,    mb * sizeof(complex_t));

    err = fft_cmplx(pa, nfft, pfft, pA);
    if(err != RES_OK)
      goto exit_label;

    err = fft_cmplx(pb, nfft, pfft, pB);
    if(err != RES_OK)
      goto exit_label;

    for(n = 0; n < nfft; n++)
    {
      RE(pC[n]) = CMRE(pA[n], pB[n]);
      IM(pC[n]) = CMIM(pA[n], pB[n]);
    }

    err = ifft_cmplx(pC, nfft, pfft, pc);
    if(err != RES_OK)
      goto exit_label;

    memcpy(c, pc+nfft2, nfft2*sizeof(complex_t));

    npos = 0;
    while(npos < ma)
    {
      if(npos+nfft > ma)
      {
        memset(pa, 0, nfft * sizeof(complex_t));
        memcpy(pa, ta+npos, (ma - npos) * sizeof(complex_t));
        err = fft_cmplx(pa, nfft, pfft, pA);


      }
      else
        err = fft_cmplx(ta+npos, nfft, pfft, pA);
      if(err != RES_OK)
        goto exit_label;
      for(n = 0; n < nfft; n++)
      {
        RE(pC[n]) = CMRE(pA[n], pB[n]);
        IM(pC[n]) = CMIM(pA[n], pB[n]);
      }

      err = ifft_cmplx(pC, nfft, pfft, pc);
      if(err != RES_OK)
        goto exit_label;
      if(npos+nfft <= ma+mb-1)
        memcpy(c+npos+nfft2, pc+nfft2,
            nfft2*sizeof(complex_t));
      else
      {
        if(ma+mb-1-npos-nfft2 > 0)
        {
          memcpy(c+npos+nfft2, pc+nfft2,(ma+mb-1-npos-nfft2)*sizeof(complex_t));
        }
      }
      npos+=nfft2;
    }
  }
  else
  {
    nfft = 4;
    n = ma - 1;
    while(n>>=1)
      nfft <<= 1;

    pa = (complex_t*)malloc(nfft * sizeof(complex_t));
    pb = (complex_t*)malloc(nfft * sizeof(complex_t));
    pc = (complex_t*)malloc(nfft * sizeof(complex_t));
    pA = (complex_t*)malloc(nfft * sizeof(complex_t));
    pB = (complex_t*)malloc(nfft * sizeof(complex_t));
    pC = (complex_t*)malloc(nfft * sizeof(complex_t));


    memset(pa, 0, nfft*sizeof(complex_t));
    memset(pb, 0, nfft*sizeof(complex_t));

    memcpy(pa, ta, ma * sizeof(complex_t));
    memcpy(pb, tb, mb * sizeof(complex_t));

    err = fft_cmplx(pa, nfft, pfft, pA);
    if(err != RES_OK)
      goto exit_label;
    err = fft_cmplx(pb, nfft, pfft, pB);
    if(err != RES_OK)
      goto exit_label;
    for(n = 0; n < nfft; n++)
    {
      RE(pC[n]) = CMRE(pA[n], pB[n]);
      IM(pC[n]) = CMIM(pA[n], pB[n]);
    }
    err = ifft_cmplx(pC, nfft, pfft, pc);
    if(err != RES_OK)
      goto exit_label;
    memcpy(c, pc, (ma+mb-1)*sizeof(complex_t));
  }
exit_label:
  if(pa)  free(pa);
  if(pb)  free(pb);
  if(pc)  free(pc);
  if(pA)  free(pA);
  if(pB)  free(pB);
  if(pB)  free(pC);

  return err;
}





/*******************************************************************************
IIR FILTER for real vector
*******************************************************************************/
int DSPL_API filter_iir(double* b, double* a, int ord,
                        double* x, int n, double* y)
{
  double* buf = NULL;
  double* an  = NULL;
  double  u;
  int   k;
  int    m;
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

