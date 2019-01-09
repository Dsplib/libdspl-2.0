/*
* Copyright (c) 2015-2019 Sergey Bakhurin
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



/******************************************************************************
Acos complex
*******************************************************************************/
int DSPL_API acos_cmplx(complex_t* x, int n, complex_t *y)
{
  int k, res;
  double pi2 = 0.5 * M_PI;

  res = asin_cmplx(x, n, y);
  if(res != RES_OK)
    return res;

  for(k = 0; k < n; k++)
  {
    RE(y[k]) = pi2 - RE(y[k]);
    IM(y[k]) =   - IM(y[k]);
  }
  return RES_OK;
}




/******************************************************************************
Asin complex
*******************************************************************************/
int DSPL_API asin_cmplx(complex_t* x, int n, complex_t *y)
{
  int k;
  complex_t tmp;
  if(!x || !y)
    return ERROR_PTR;
  if(n < 1)
    return ERROR_SIZE;

  for(k = 0; k < n; k++)
  {
    RE(tmp) = 1.0 - CMRE(x[k], x[k]); // 1-x[k]^2
    IM(tmp) =   - CMIM(x[k], x[k]); // 1-x[k]^2
    sqrt_cmplx(&tmp, 1, y+k);     // sqrt(1 - x[k]^2)
    RE(y[k]) -= IM(x[k]);       // j * x[k] + sqrt(1 - x[k]^2)
    IM(y[k]) += RE(x[k]);       // j * x[k] + sqrt(1 - x[k]^2)
    log_cmplx(y+k, 1, &tmp);      // log( j * x[k] + sqrt(1 - x[k]^2) )
    RE(y[k]) =  IM(tmp);        // -j * log( j * x[k] + sqrt(1 - x[k]^2) )
    IM(y[k]) = -RE(tmp);        // -j * log( j * x[k] + sqrt(1 - x[k]^2) )
  }
  return RES_OK;
}



/******************************************************************************
convert double array to a complex array
*******************************************************************************/
int DSPL_API re2cmplx(double* x, int n, complex_t *y)
{
  int k;
  if(!x || !y)
    return ERROR_PTR;
  if(n < 1)
    return ERROR_SIZE;

  for(k = 0; k < n; k++)
  {
    RE(y[k]) = x[k];
    IM(y[k]) = 0.0;
  }
  return RES_OK;
}





/******************************************************************************
convert complex array to a re and im arrays
*******************************************************************************/
int DSPL_API cmplx2re(complex_t* x, int n, double *re, double *im)
{
  int k;
  if(!x)
    return ERROR_PTR;
  if(n < 1)
    return ERROR_SIZE;

  if(re)
  {
    for(k = 0; k < n; k++)
      re[k] = RE(x[k]);
  }
  if(im)
  {
    for(k = 0; k < n; k++)
      im[k] = IM(x[k]);
  }
  return RES_OK;
}






/******************************************************************************
Complex cosine
*******************************************************************************/
int DSPL_API cos_cmplx(complex_t* x, int n, complex_t *y)
{
  int k;
  double ep, em, sx, cx;
  if(!x || !y)
    return ERROR_PTR;
  if(n < 1)
    return ERROR_SIZE;

  for(k = 0; k < n; k++)
  {
    ep = exp( IM(x[k]));
    em = exp(-IM(x[k]));
    sx = 0.5 * sin(RE(x[k]));
    cx = 0.5 * cos(RE(x[k]));
    RE(y[k]) = cx * (em + ep);
    IM(y[k]) = sx * (em - ep);
  }
  return RES_OK;
}




/******************************************************************************
Complex cosine
*******************************************************************************/
int DSPL_API sin_cmplx(complex_t* x, int n, complex_t *y)
{
  int k;
  double ep, em, sx, cx;
  if(!x || !y)
    return ERROR_PTR;
  if(n < 1)
    return ERROR_SIZE;

  for(k = 0; k < n; k++)
  {
    ep = exp( IM(x[k]));
    em = exp(-IM(x[k]));
    sx = 0.5 * sin(RE(x[k]));
    cx = 0.5 * cos(RE(x[k]));
    RE(y[k]) = sx * (em + ep);
    IM(y[k]) = cx * (ep - em);
  }
  return RES_OK;
}




/******************************************************************************
Logarithm complex
*******************************************************************************/
int DSPL_API log_cmplx(complex_t* x, int n, complex_t *y)
{
  int k;
  if(!x || !y)
    return ERROR_PTR;
  if(n < 1)
    return ERROR_SIZE;

  for(k = 0; k < n; k++)
  {
    RE(y[k]) = 0.5 * log(ABSSQR(x[k]));
    IM(y[k]) = atan2(IM(x[k]), RE(x[k]));
  }
  return RES_OK;
}




/******************************************************************************
SQRT complex
*******************************************************************************/
int DSPL_API sqrt_cmplx(complex_t* x, int n, complex_t *y)
{
  int k;
  double r, zr;
  complex_t t;
  if(!x || !y)
    return ERROR_PTR;
  if(n < 1)
    return ERROR_SIZE;

  for(k = 0; k < n; k++)
  {
    r = ABS(x[k]);
    RE(t) = RE(x[k]) + r;
    IM(t) = IM(x[k]);
    zr = 1.0 / ABS(t);
    r = sqrt(r);
    RE(y[k]) = RE(t) * zr * r;
    IM(y[k]) = IM(t) * zr * r;

  }
  return RES_OK;
}

