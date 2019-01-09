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
#include <math.h>
#include "dspl.h"




/******************************************************************************
Complex frequency response of an analog filter H(s)
*******************************************************************************/
int DSPL_API freqs(double* b, double* a, int ord,
                   double* w, int n, complex_t *h)
{
  complex_t jw;
  complex_t *bc = NULL;
  complex_t *ac = NULL;
  complex_t num, den;
  double mag;
  int k;
  int res;

  if(!b || !a || !w || !h)
    return ERROR_PTR;
  if(ord<0)
    return ERROR_FILTER_ORD;
  if(n<1)
    return ERROR_SIZE;

  RE(jw) = 0.0;

  bc = (complex_t*) malloc((ord+1) * sizeof(complex_t));
  res = re2cmplx(b, ord+1, bc);

  if( res!=RES_OK )
    goto exit_label;

  ac = (complex_t*) malloc((ord+1) * sizeof(complex_t));
  res = re2cmplx(a, ord+1, ac);
  if( res!=RES_OK )
    goto exit_label;

  for(k = 0; k < n; k++)
  {
    IM(jw) = w[k];
    res = polyval_cmplx(bc, ord, &jw, 1, &num);
    if(res != RES_OK)
      goto exit_label;
    res = polyval_cmplx(ac, ord, &jw, 1, &den);
    if(res != RES_OK)
      goto exit_label;
    mag = ABSSQR(den);
    if(mag == 0.0)
    {
      res = ERROR_DIV_ZERO;
      goto exit_label;
    }
    mag = 1.0 / mag;
    RE(h[k]) = CMCONJRE(num, den) * mag;
    IM(h[k]) = CMCONJIM(num, den) * mag;
  }
  res = RES_OK;
exit_label:
  if(bc)
    free(bc);
  if(ac)
    free(ac);
  return res;
}






/******************************************************************************
 * Complex frequency response of an analog filter H(s), s is complex variable
 ******************************************************************************/
int DSPL_API freqs_cmplx(double* b, double* a, int ord,
                         complex_t* s, int n, complex_t *h)
{
  complex_t *bc = NULL;
  complex_t *ac = NULL;
  complex_t num, den;
  double mag;
  int k;
  int res;

  if(!b || !a || !s || !h)
    return ERROR_PTR;
  if(ord<0)
    return ERROR_FILTER_ORD;
  if(n<1)
    return ERROR_SIZE;


  bc = (complex_t*) malloc((ord+1) * sizeof(complex_t));
  res = re2cmplx(b, ord+1, bc);

  if( res!=RES_OK )
    goto exit_label;

  ac = (complex_t*) malloc((ord+1) * sizeof(complex_t));
  res = re2cmplx(a, ord+1, ac);
  if( res!=RES_OK )
    goto exit_label;

  for(k = 0; k < n; k++)
  {
    res = polyval_cmplx(bc, ord, s+k, 1, &num);
    if(res != RES_OK)
      goto exit_label;
    res = polyval_cmplx(ac, ord, s+k, 1, &den);
    if(res != RES_OK)
      goto exit_label;
    mag = ABSSQR(den);
    if(mag == 0.0)
    {
      res = ERROR_DIV_ZERO;
      goto exit_label;
    }
    mag = 1.0 / mag;
    RE(h[k]) = CMCONJRE(num, den) * mag;
    IM(h[k]) = CMCONJIM(num, den) * mag;

  }
  res = RES_OK;
  exit_label:
  if(bc)
    free(bc);
  if(ac)
    free(ac);
  return res;
}







/******************************************************************************
impulse response of an analog filter H(s)
*******************************************************************************/
int DSPL_API freqs2time(double* b, double* a, int ord, double fs,
                        int n, fft_t* pfft, double *t, double *h)
{
  double *w = NULL;
  complex_t *hs = NULL;
  complex_t *ht = NULL;
  int err, k;

  if(!b || !a || !t || !h)
    return ERROR_PTR;
  if(ord<1)
    return ERROR_FILTER_ORD;
  if(n<1)
    return ERROR_SIZE;

  w  = (double*)malloc(n*sizeof(double));
  hs = (complex_t*)malloc(n*sizeof(complex_t));


  err = linspace(-fs*0.5, fs*0.5, n, DSPL_PERIODIC, w);
  if(err != RES_OK)
    goto exit_label;

  err = freqs(b, a, ord, w, n, hs);
  if(err != RES_OK)
    goto exit_label;

  err = fft_shift_cmplx(hs, n, hs);
  if(err != RES_OK)
    goto exit_label;

  ht = (complex_t*)malloc(n*sizeof(complex_t));

  err = ifft_cmplx(hs, n, pfft, ht);
  if(err != RES_OK)
  {
    err = idft_cmplx(hs, n, ht);
    if(err != RES_OK)
      goto exit_label;
  }

  for(k = 0; k < n; k++)
  {
    t[k] = (double)k/fs;
    h[k] = RE(ht[k]) * fs;
  }

exit_label:
  if(w)
    free(w);
  if(hs)
    free(hs);
  if(ht)
    free(ht);
  return err;
}





/******************************************************************************
Magnitude, phase response and group delay of an analog filter H(s)
*******************************************************************************/
int DSPL_API freqs_resp(double* b, double* a, int ord,
                        double* w, int n, int flag,
                        double *h, double* phi, double* tau)
{
  int res, k;

  complex_t *hc = NULL;
  double *phi0 = NULL;
  double *phi1 = NULL;
  double *w0   = NULL;
  double *w1   = NULL;

  if(!b || !a || !w)
    return ERROR_PTR;
  if(ord < 1)
    return ERROR_FILTER_ORD;
  if(n < 1)
    return ERROR_SIZE;


  hc = (complex_t*) malloc (n*sizeof(complex_t));
  res = freqs(b, a, ord, w, n, hc);
  if(res != RES_OK)
    goto exit_label;


  if(h)
  {
    if(flag & DSPL_FLAG_LOG)
    {
      for(k = 0; k < n; k++)
        h[k] = 10.0 * log10(ABSSQR(hc[k]));
    }
    else
    {
      for(k = 0; k < n; k++)
        h[k] = sqrt(ABSSQR(hc[k]));
    }
  }


  if(phi)
  {
    for(k = 0; k < n; k++)
      phi[k] = atan2(IM(hc[k]), RE(hc[k]));

    if(flag & DSPL_FLAG_UNWRAP)
    {
      res = unwrap(phi, n, M_2PI, 0.8);
      if(res != RES_OK)
        goto exit_label;
    }
  }


  if(tau)
  {
    phi0 = (double*) malloc(n*sizeof(double));
    phi1 = (double*) malloc(n*sizeof(double));
    w0   = (double*) malloc(n*sizeof(double));
    w1   = (double*) malloc(n*sizeof(double));

    w0[0] = w[0] - (w[1] - w[0])*0.02;
    w1[0] = w[0] + (w[1] - w[0])*0.02;

    for(k = 1; k < n; k++)
    {
      w0[k] = w[k] - (w[k] - w[k-1])*0.02;
      w1[k] = w[k] + (w[k] - w[k-1])*0.02;
    }
    res = freqs_resp(b, a, ord, w0, n, DSPL_FLAG_UNWRAP, NULL, phi0, NULL);
    if(res != RES_OK)
      goto exit_label;
    res = freqs_resp(b, a, ord, w1, n, DSPL_FLAG_UNWRAP, NULL, phi1, NULL);
    if(res != RES_OK)
      goto exit_label;
    for(k = 0; k < n; k++)
        tau[k] = (phi0[k] - phi1[k])/(w1[k] - w0[k]);
  }


exit_label:
  if(hc)
    free(hc);
  if(phi0)
    free(phi0);
  if(phi1)
    free(phi1);
  if(w0)
    free(w0);
  if(w1)
    free(w1);
  return res;
}






/*******************************************************************************
Complex frequency response of a digital filter H(z)
*******************************************************************************/
int DSPL_API freqz(double* b, double* a, int ord, double* w,
              int n, complex_t *h)
{
  complex_t jw;
  complex_t *bc = NULL;
  complex_t *ac = NULL;
  complex_t num, den;
  double mag;
  int k;
  int res;

  if(!b || !w || !h)
    return ERROR_PTR;
  if(ord<0)
    return ERROR_FILTER_ORD;
  if(n<1)
    return ERROR_SIZE;


  bc = (complex_t*) malloc((ord+1) * sizeof(complex_t));
  res = re2cmplx(b, ord+1, bc);
  if( res!=RES_OK )
    goto exit_label;

  if(a)
  {
    // IIR filter if a != NULL
    ac = (complex_t*) malloc((ord+1) * sizeof(complex_t));
    res = re2cmplx(a, ord+1, ac);
    if( res!=RES_OK )
      goto exit_label;
    for(k = 0; k < n; k++)
    {
      RE(jw) =  cos(w[k]);
      IM(jw) = -sin(w[k]);
      res = polyval_cmplx(bc, ord, &jw, 1, &num);
      if(res != RES_OK)
        goto exit_label;
      res = polyval_cmplx(ac, ord, &jw, 1, &den);
      if(res != RES_OK)
        goto exit_label;
      mag = ABSSQR(den);
      if(mag == 0.0)
      {
        res = ERROR_DIV_ZERO;
        goto exit_label;
      }
      mag = 1.0 / mag;
      RE(h[k]) = CMCONJRE(num, den) * mag;
      IM(h[k]) = CMCONJIM(num, den) * mag;
    }
  }
  else
  {
    // FIR filter if a == NULL
    for(k = 0; k < n; k++)
    {
      RE(jw) =  cos(w[k]);
      IM(jw) = -sin(w[k]);
      res = polyval_cmplx(bc, ord, &jw, 1, h+k);
      if(res != RES_OK)
        goto exit_label;
    }
  }
  res = RES_OK;
exit_label:
  if(bc)
    free(bc);
  if(ac)
    free(ac);
  return res;
}








/*******************************************************************************
Unwrap function
*******************************************************************************/
int DSPL_API unwrap(double* phi, int n, double lev, double mar)
{
  double a[2] = {0.0, 0.0};
  double d;
  double th;
  int k;
  int flag = 1;


  if(!phi)
    return ERROR_PTR;

  if(n<1)
    return ERROR_SIZE;

  if(lev<=0 || mar <=0)
    return ERROR_UNWRAP;

  th = mar*lev;
  while(flag)
  {
    flag = 0;
    a[0] = a[1] = 0.0;
    for(k = 0; k<n-1; k++)
    {
      d = phi[k+1] - phi[k];
      if( d > th)
      {
        a[0] -= lev;
        flag = 1;
      }
      if( d < -th)
      {
        a[0] += lev;
        flag = 1;
      }
      phi[k]+=a[1];
      a[1] = a[0];
    }
    phi[n-1]+=a[1];
  }

    return RES_OK;
}

