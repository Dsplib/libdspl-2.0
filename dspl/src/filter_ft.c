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

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "dspl.h"


/******************************************************************************
 * Recalculate ws frequency to 1 rad/s for HPF and BANDSTOP filters
 ******************************************************************************/
double DSPL_API filter_ws1(int ord, double rp, double rs, int type)
{
  double es2, ep2, gs2, x, ws;
  
  if(ord<1 || rp < 0.0 || rs < 0.0)
    return -1.0;
  
  es2 = pow(10.0, rs*0.1) - 1.0;
  ep2 = pow(10.0, rp*0.1) - 1.0;
  gs2 = 1.0 / (1.0 + es2);
  
  x = (1.0 - gs2) / (gs2 * ep2);
  
  switch( type & DSPL_FILTER_APPROX_MASK)
  {
    case DSPL_FILTER_BUTTER:
      ws = pow(x, 0.5 / (double)ord);
      break;
    case DSPL_FILTER_CHEBY1:
    case DSPL_FILTER_CHEBY2:
      x = sqrt(x) + sqrt(x - 1.0);
      x = log(x) / (double)ord;
      ws  = 0.5 * (exp(-x) + exp(x));
      break;
    case DSPL_FILTER_ELLIP:
    {
      double k, k1;
      complex_t y, z;
      int res;
      k = sqrt(ep2 / es2);
      res = ellip_modulareq(rp, rs, ord, &k1);
      if(res != RES_OK)
      {
        ws = -1.0;
        break;
      }
      RE(z) = sqrt(x);
      IM(z) = 0.0;
      
      res = ellip_acd_cmplx(&z, 1, k, &y);
      if(res != RES_OK)
      {
        ws = -1.0;
        break;
      }
      RE(y) /= (double)ord;
      IM(y) /= (double)ord;
      res = ellip_cd_cmplx(&y, 1, k1, &z);
      if(res != RES_OK)
      {
        ws = -1.0;
        break;
      }
      ws = RE(z);
      break;
    }
    default:
      ws  = -1.0;
      break;
  }
  return ws;
}

/******************************************************************************
 * low 2 bandpass transformation
 ******************************************************************************/
int DSPL_API low2bp(double* b, double* a, int ord,
                    double w0, double wpl, double wph,
                    double* beta, double* alpha)
{

  double num[3] = {0.0, 0.0, 1.0};
  double den[3] = {0.0, 0.0, 0.0};

  if(!b || !a || !beta || !alpha)
    return ERROR_PTR;
  if(ord < 1)
    return ERROR_FILTER_ORD;
  if(w0 <= 0.0 || wpl <= 0.0 || wph <= 0.0 || wph <= wpl)
    return ERROR_FILTER_FT;

  num[0] = (wph * wpl) / (w0 * w0);
  den[1] = (wph - wpl) / w0;

  return ratcompos(b, a, ord, num, den, 2, beta, alpha);
}





/******************************************************************************
 * low 2 bandstop transformation
 ******************************************************************************/
int DSPL_API low2bs(double* b, double* a, int ord,
                    double w0, double wsl, double wsh,
                    double* beta, double* alpha)
{

  double den[3] = {0.0, 0.0, 1.0};
  double num[3] = {0.0, 0.0, 0.0};

  if(!b || !a || !beta || !alpha)
    return ERROR_PTR;
  if(ord < 1)
    return ERROR_FILTER_ORD;
  if(w0 <= 0.0 || wsl <= 0.0 || wsh <= 0.0 || wsh <= wsl)
    return ERROR_FILTER_FT;

  den[0] = (wsh * wsl) / (w0 * w0);
  num[1] = (wsh - wsl) / w0;

  return ratcompos(b, a, ord, num, den, 2, beta, alpha);
}




/******************************************************************************
 * low 2 high transformation
 ******************************************************************************/
int DSPL_API low2high(double* b, double* a, int ord, double w0, double w1,
                      double* beta, double* alpha)
{

  double num[2] = {0.0, 0.0};
  double den[2] = {0.0, 1.0};

  if(!b || !a || !beta || !alpha)
    return ERROR_PTR;
  if(ord < 1)
    return ERROR_FILTER_ORD;
  if(w0 <= 0.0 || w1 <= 0.0)
    return ERROR_FILTER_FT;

  num[0] = w1 / w0;

  return ratcompos(b, a, ord, num, den, 1, beta, alpha);
}







/******************************************************************************
 low 2 low transformation
*******************************************************************************/
int DSPL_API low2low(double* b, double* a, int ord, double w0, double w1,
                     double* beta, double* alpha)
{

  double num[2] = {0.0, 1.0};
  double den[2] = {0.0, 0.0};

  if(!b || !a || !beta || !alpha)
    return ERROR_PTR;
  if(ord < 1)
    return ERROR_FILTER_ORD;
  if(w0 <= 0.0 || w1 <= 0.0)
    return ERROR_FILTER_FT;

  den[0] = w1 / w0;

  return ratcompos(b, a, ord, num, den, 1, beta, alpha);
}




/******************************************************************************
Rational composition
*******************************************************************************/
int DSPL_API ratcompos(double* b, double* a, int n,
                       double* c, double* d, int p,
                       double* beta, double* alpha)
{

  int k2, i, k,  pn, pd, ln, ld, k2s, nk2s;
  double *num = NULL, *den = NULL, *ndn = NULL, *ndd = NULL;
  int res;

  if (!a || !b || !c || !d || !beta || !alpha)
  {
    res = ERROR_PTR;
    goto exit_label;
  }
  if(n < 1 || p < 1)
  {
    res =  ERROR_SIZE;
    goto exit_label;
  }

  k2   = (n*p)+1;
  k2s  = k2*sizeof(double);       /* alpha and beta size    */
  nk2s = (n+1)*k2*sizeof(double); /* num, den, ndn and ndd size  */

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
      ndd[i*k2 + k] = ndn[i*k2 + k] * a[i];
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

