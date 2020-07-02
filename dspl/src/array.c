/*
* Copyright (c) 2015-2020 Sergey Bakhurin
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
#include "blas.h"


/******************************************************************************
Vector linear transformation 
--------------------------------------------------------------------------------
Documented: RU, EN
*******************************************************************************/
int array_scale_lin(double* x,   int n, 
    double xmin, double xmax, double dx,
    double h,    double* y)
{
  double kx;
  int k;  
  if(!x)
    return ERROR_PTR;
  if(n<1)
    return ERROR_SIZE;
  if(h<0.0)
    return ERROR_NEGATIVE;

  if(xmin >= xmax)
    return ERROR_MIN_MAX;

  kx = h / (xmax - xmin);

  for(k = 0; k < n; k++)
    y[k] = (x[k] - xmin) * kx + dx;

  return RES_OK;
}




/******************************************************************************
Concatenate arrays `a` and `b`
--------------------------------------------------------------------------------
Documented: RU, EN
*******************************************************************************/
int DSPL_API concat(void* a, size_t na, void* b, size_t nb, void* c)
{
  if(!a || !b || !c || c == b)
    return ERROR_PTR;
  if(na < 1 || nb < 1)
    return ERROR_SIZE;

  if(c != a)
    memcpy(c, a, na);

  memcpy((char*)c+na, b, nb);
  return RES_OK;
}





/******************************************************************************
Real vector decimation
--------------------------------------------------------------------------------
Documented: RU, EN
*******************************************************************************/
int DSPL_API decimate(double* x, int n, int d, double* y, int* cnt)
{
  int k = 0, i = 0;
  if(!x || !y)
    return ERROR_PTR;
  if(n < 1)
    return ERROR_SIZE;
  if(d < 1)
    return ERROR_NEGATIVE;

  k = i = 0;
  while(k + d <= n)
  {
    y[i] = x[k];
    k+=d;
    i++;
  }
  if(cnt)
    *cnt = i;

  return RES_OK;
}




/******************************************************************************
Complex vector decimation
--------------------------------------------------------------------------------
Documented: RU, EN
*******************************************************************************/
int DSPL_API decimate_cmplx(complex_t* x, int n, int d, complex_t* y, int* cnt)
{
  int k = 0, i = 0;
  if(!x || !y)
    return ERROR_PTR;
  if(n < 1)
    return ERROR_SIZE;
  if(d < 1)
    return ERROR_NEGATIVE;

  k = i = 0;
  while(k + d < n)
  {
    RE(y[i]) = RE(x[k]);
    IM(y[i]) = IM(x[k]);
    k+=d;
    i++;
  }
  if(cnt)
    *cnt = i;

  return RES_OK;
}





/******************************************************************************
Flip real vector `x` in place
--------------------------------------------------------------------------------
Documented: RU, EN
*******************************************************************************/
int DSPL_API flipip(double* x, int n)
{
  int k;
  double tmp;
  if(!x)
    return ERROR_PTR;
  if(n<1)
    return ERROR_SIZE;

  for(k = 0; k < n/2; k++)
  {
    tmp = x[k];
    x[k] = x[n-1-k];
    x[n-1-k] = tmp;
  }
  return RES_OK;
  
}



/******************************************************************************
Flip complex vector `x` in place
--------------------------------------------------------------------------------
Documented: RU, EN
*******************************************************************************/
int DSPL_API flipip_cmplx(complex_t* x, int n)
{
  int k;
  complex_t tmp;
  if(!x)
    return ERROR_PTR;
  if(n<1)
    return ERROR_SIZE;

  for(k = 0; k < n/2; k++)
  {
    RE(tmp) = RE(x[k]);
    RE(x[k]) = RE(x[n-1-k]);
    RE(x[n-1-k]) = RE(tmp);

    IM(tmp) = IM(x[k]);
    IM(x[k]) = IM(x[n-1-k]);
    IM(x[n-1-k]) = IM(tmp);
  }
  return RES_OK;
}






/*******************************************************************************
Linspace array filling
--------------------------------------------------------------------------------
Documented: RU, EN
*******************************************************************************/
int DSPL_API linspace(double x0, double x1, int n, int type, double* x)
{
  double dx;
  int k;

  if(n < 2)
    return ERROR_SIZE;
  if(!x)
  return ERROR_PTR;

  switch (type)
  {
    case DSPL_SYMMETRIC:
      dx = (x1 - x0)/(double)(n-1);
      x[0] = x0;
      for(k = 1; k < n; k++)
        x[k] = x[k-1] + dx;
      break;
    case DSPL_PERIODIC:
      dx = (x1 - x0)/(double)n;
      x[0] = x0;
      for(k = 1; k < n; k++)
        x[k] = x[k-1] + dx;
      break;
    default:
      return ERROR_SYM_TYPE;
  }
  return RES_OK;
}





/*******************************************************************************
Logspace array filling
--------------------------------------------------------------------------------
Documented: RU, EN
*******************************************************************************/
int DSPL_API logspace(double x0, double x1, int n, int type, double* x)
{
  double mx, a, b;
  int k;

  if(n < 2)
    return ERROR_SIZE;
  if(!x)
    return ERROR_PTR;

  a = pow(10.0, x0);
  b = pow(10.0, x1);

  switch (type)
  {
    case DSPL_SYMMETRIC:
      mx = pow(b/a, 1.0/(double)(n-1));
      x[0] = a;
      for(k = 1; k < n; k++)
        x[k] = x[k-1] * mx;
      break;
    case DSPL_PERIODIC:
      mx = pow(b/a, 1.0/(double)n);
      x[0] = a;
      for(k = 1; k < n; k++)
        x[k] = x[k-1] * mx;
      break;
    default:
      return ERROR_SYM_TYPE;
  }
  return RES_OK;
}


/*******************************************************************************
Ones double array
--------------------------------------------------------------------------------
Documented: RU, EN
*******************************************************************************/
int DSPL_API ones(double* x, int n)
{
  int i;
  if(!x)
    return ERROR_PTR;
  if(n<1)
    return ERROR_SIZE;
  for(i = 0; i < n; i++)
    x[i] = 1.0;
 return RES_OK;
}


/******************************************************************************
Real arrays verification
--------------------------------------------------------------------------------
Documented: RU, EN
*******************************************************************************/
int DSPL_API verif(double* x,  double* y, size_t n, double eps, double* err)
{
  double d, maxd; 
  size_t k; 
  int res;
  if(!x || !y)
    return ERROR_PTR;
  if(n < 1)
    return ERROR_SIZE;
  if(eps <= 0.0 )
    return ERROR_NEGATIVE;
    
  maxd = -100.0;
  
  for(k = 0; k < n; k++)
  {
    d = fabs(x[k] - y[k]);
    if(fabs(x[k]) > 0.0)
    {
      d = d / fabs(x[k]);
      if(d > maxd)
        maxd = d;
    }
  }
  if(err) 
    *err = maxd;
    
  if(maxd > eps)
    res = DSPL_VERIF_FAILED;
  else
    res = DSPL_VERIF_SUCCESS;
 
  return res;
}



/******************************************************************************
Complex arrays verification
--------------------------------------------------------------------------------
Documented: RU, EN
*******************************************************************************/
int DSPL_API verif_cmplx(complex_t* x,  complex_t* y, size_t n, 
         double eps, double* err)
{
  
  complex_t d;
  double mx, md, maxd;
  size_t k;
  int res;
  if(!x || !y)
    return ERROR_PTR;
  if(n < 1)
    return ERROR_SIZE;
  if(eps <= 0.0 )
    return ERROR_NEGATIVE;
    
  maxd = -100.0;
  
  for(k = 0; k < n; k++)
  {
    RE(d) = RE(x[k]) - RE(y[k]);
    IM(d) = IM(x[k]) - IM(y[k]);
    md = ABS(d);
    mx = ABS(x[k]);
    if(mx > 0.0)
    {
      md = md / mx;
      if(md > maxd)
        maxd = md;
    }
  }
  if(err)
    *err = maxd;
    
  if(maxd > eps)
    res = DSPL_VERIF_FAILED;
  else
    res = DSPL_VERIF_SUCCESS;
 
  return res;
}
