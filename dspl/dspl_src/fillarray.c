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
#include <math.h>
#include "dspl.h"



/*******************************************************************************
Linspace array filling
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

