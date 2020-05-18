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
Polynomial zeros to coefficients
*******************************************************************************/
int DSPL_API poly_z2a_cmplx(complex_t* z, int nz, int ord, complex_t* a)
{
  int k, ind, res;
  complex_t x[2];

  if(!z || !a)
    return ERROR_PTR;
  if(nz < 0)
    return ERROR_SIZE;
  if(nz > ord || ord < 1)
    return ERROR_POLY_ORD;

  RE(x[1]) = 1.0;
  IM(x[1]) = 0.0;

  memset(a, 0, (ord+1) * sizeof(complex_t));

  RE(a[0]) = 1.0;
  ind = 1;
  for(k = 0; k < nz; k++)
  {
    RE(x[0]) = -RE(z[k]);
    IM(x[0]) = -IM(z[k]);
    res = conv_cmplx(a, ind, x, 2, a);
    if(res!=RES_OK)
      return res;
    ind++;
  }

  return RES_OK;
}





/******************************************************************************
Real polynomial roots calculation
*******************************************************************************/
int DSPL_API polyroots(double* a, int ord, complex_t* r, int* info)
{
  complex_t *t = NULL;
  int m;
  int err;
  
  if(!a || !r)
    return ERROR_PTR;
  if(ord<0)
    return ERROR_POLY_ORD;
  if(a[ord] == 0.0)
    return ERROR_POLY_AN;
  
  t = (complex_t*)malloc(ord * ord * sizeof(complex_t));
  if(!t)
    return ERROR_MALLOC;
  
  for(m = 0; m < ord-1; m++)
  {
    RE(t[m * (ord+1) + 1]) = 1.0;
    RE(t[m + ord * (ord - 1)]) = -a[m] / a[ord];
  }
  RE(t[ord * ord - 1]) = -a[ord-1] / a[ord];

  err = matrix_eig_cmplx(t, ord, r, info);
  
  return err;
}




/******************************************************************************
Real polynomial evaluation
*******************************************************************************/
int DSPL_API polyval(double* a, int ord, double* x, int n, double* y)
{
  int k, m;

  if(!a || !x || !y)
    return ERROR_PTR;
  if(ord<0)
    return ERROR_POLY_ORD;
  if(n<1)
    return ERROR_SIZE;

  for(k = 0; k < n; k++)
  {
    y[k] = a[ord];
    for(m = ord-1; m>-1; m--)
      y[k] = y[k]*x[k] + a[m];
  }
  return RES_OK;
}





/*******************************************************************************
Complex polynomial evaluation
*******************************************************************************/
int DSPL_API polyval_cmplx(complex_t* a, int ord,
                           complex_t* x, int n, complex_t* y)
{
  int k, m;
  complex_t t;

  if(!a || !x || !y)
    return ERROR_PTR;
  if(ord<0)
    return ERROR_POLY_ORD;
  if(n<1)
    return ERROR_SIZE;

  for(k = 0; k < n; k++)
  {
    RE(y[k]) = RE(a[ord]);
    IM(y[k]) = IM(a[ord]);
    for(m = ord-1; m>-1; m--)
    {
      RE(t) = CMRE(y[k], x[k]);
      IM(t) = CMIM(y[k], x[k]);
      RE(y[k]) = RE(t) + RE(a[m]);
      IM(y[k]) = IM(t) + IM(a[m]);
    }
  }
  return RES_OK;
}

