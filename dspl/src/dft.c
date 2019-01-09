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



/******************************************************************************
Real vector DFT
*******************************************************************************/
int DSPL_API dft(double* x, int n, complex_t *y)
{
  int k;
  int m;
  double divn;
  double phi;


  if(!x || !y)
    return ERROR_PTR;

  if(n<1)
    return ERROR_SIZE;

  divn = 1.0 / (double)n;

  for(k = 0; k < n; k++)
  {
    RE(y[k]) = IM(y[k]) = 0.0;
    for(m = 0; m < n; m++)
    {
      phi  = -M_2PI * divn * (double)k * (double)m;
      RE(y[k]) += x[m] * cos(phi);
      IM(y[k]) += x[m] * sin(phi);
    }
  }
  return RES_OK;
}



/******************************************************************************
Complex vector DFT
*******************************************************************************/
int DSPL_API dft_cmplx(complex_t* x, int n, complex_t *y)
{
  int k;
  int m;
  double divn;
  double phi;
  complex_t e;

  if(!x || !y)
    return ERROR_PTR;

  if(n<1)
    return ERROR_SIZE;

  divn = 1.0 / (double)n;

  for(k = 0; k < n; k++)
  {
    RE(y[k]) = IM(y[k]) = 0.0;
    for(m = 0; m < n; m++)
    {
      phi  = -M_2PI * divn * (double)k * (double)m;
      RE(e) = cos(phi);
      IM(e) = sin(phi);
      RE(y[k]) += CMRE(x[m], e);
      IM(y[k]) += CMIM(x[m], e);
    }
  }
  return RES_OK;
}





/******************************************************************************
Complex vector DFT
*******************************************************************************/
int DSPL_API idft_cmplx(complex_t* x, int n, complex_t *y)
{
  int k;
  int m;
  double divn;
  double phi;
  complex_t e;

  if(!x || !y)
    return ERROR_PTR;

  if(n<1)
    return ERROR_SIZE;

  divn = 1.0 / (double)n;

  for(k = 0; k < n; k++)
  {
    RE(y[k]) = IM(y[k]) = 0.0;
    for(m = 0; m < n; m++)
    {
      phi  =  M_2PI * divn * (double)k * (double)m;
      RE(e) = cos(phi);
      IM(e) = sin(phi);
      RE(y[k]) += CMRE(x[m], e);
      IM(y[k]) += CMIM(x[m], e);
    }
    RE(y[k]) /= (double)n;
    IM(y[k]) /= (double)n;
  }
  return RES_OK;
}

