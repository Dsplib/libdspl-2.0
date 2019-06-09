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
\fn int concat(void* a, size_t na, void* b, size_t nb, void* c) 
\brief Concatenate arrays `a` and `b`

Let's arrays `a` and `b` are vectors:<BR>
`a = [a(0), a(1), ... a(na-1)]`, <BR>
`b = [b(0), b(1), ... b(nb-1)]`, <BR>
concatenation of these arrays will be array `c` size `na+nb`:<BR>
`c = [a(0), a(1), ... a(na-1), b(0), b(1), ... b(nb-1)]`.


\param[in]  a   Pointer to the first array `a`.<BR>
                Array `a` size is `na` bytes.<BR><BR>

\param[in]  na  Array `a` size (bytes).<BR><BR>

\param[in]  b   Pointer to the second array `b`.<BR>
                Array `b` size is `nb` bytes.<BR><BR>

\param[in]  nb  Array `a` size (bytes).<BR><BR>

\param[out] c   Pointer to the concatenation result array `c`.<BR>
                Array `c` size is `na + nb` bytes.<BR>
                Memory must be allocated.<BR><BR>

\return
`RES_OK` if function returns successfully.<BR>
 Else \ref ERROR_CODE_GROUP "code error".
 
Function uses pointer type `void*` and can be useful for an arrays 
concatenation with different types.<BR>
For example two `double` arrays concatenation:
\code
double a[3] = {1.0, 2.0, 3.0};
double b[2] = {4.0, 5.0};
double c[5];

concat((void*)a, 3*sizeof(double), (void*)b, 2*sizeof(double), (void*)c);
\endcode 
Vector `c` keeps follow data:
\code
c = [1.0, 2.0, 3.0, 4.0, 5.0]
\endcode 

\author Sergey Bakhurin www.dsplib.org
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
decimate real vector
*******************************************************************************/
int DSPL_API decimate(double* x, int n, int dec, double* y, int* cnt)
{
  int k = 0, i = 0;
  if(!x || !y)
    return ERROR_PTR;
  if(n < 1)
    return ERROR_SIZE;
  if(dec < 1)
    return ERROR_NEGATIVE;

  k = i = 0;
  while(k + dec < n)
  {
    y[i] = x[k];
    k+=dec;
    i++;
  }
  if(cnt)
    *cnt = i;

  return RES_OK;
}




/******************************************************************************
decimate complex vector
*******************************************************************************/
int DSPL_API decimate_cmplx(complex_t* x, int n, int dec,
                            complex_t* y, int* cnt)
{
  int k = 0, i = 0;
  if(!x || !y)
    return ERROR_PTR;
  if(n < 1)
    return ERROR_SIZE;
  if(dec < 1)
    return ERROR_NEGATIVE;

  k = i = 0;
  while(k + dec < n)
  {
    RE(y[i]) = RE(x[k]);
    IM(y[i]) = IM(x[k]);
    k+=dec;
    i++;
  }
  if(cnt)
    *cnt = i;

  return RES_OK;
}

/******************************************************************************
Find max(|a|)
*******************************************************************************/
int DSPL_API find_max_abs(double* a, int n, double* m, int* ind)
{
  int k, i;
  double t;
  if(!a)
    return ERROR_PTR;
  if(n < 1)
    return ERROR_SIZE;
  t = fabs(a[0]);
  i = 0;
  for(k = 1; k < n; k++)
  {
    if(fabs(a[k]) > t)
    {
      t = a[k];
      i = k;
    }
  }
  if(m)
    *m = t;
  if(ind)
    *ind = i;
  return RES_OK;
}


/******************************************************************************
Flip real array in place
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
Flip complex array in place
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




/******************************************************************************
Verif double
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
Verif double
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
