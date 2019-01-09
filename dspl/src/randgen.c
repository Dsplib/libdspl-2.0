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
#include <time.h>

#include "dspl.h"
#include "dspl_internal.h"




/******************************************************************************
Uniform random numbers generator
*******************************************************************************/
int DSPL_API randu(double* x, int n)
{
   int k,m;
   unsigned int x1[4], x2[4], y;

   if(!x)
    return  ERROR_PTR;
   if(n<1)
    return ERROR_SIZE;

   x1[1] = rand();
   x2[1] = rand();
   x1[2] = rand();
   x2[2] = rand();
   x1[3] = rand();
   x2[3] = rand();
   for(k = 0; k<n; k++)
   {
     x1[0] = (63308 * x1[2] - 183326*x1[3]) % DSPL_RAND_MOD_X1;
     x2[0] = (86098 * x2[1] - 539608*x2[3]) % DSPL_RAND_MOD_X2;
     y = (x1[0] - x2[0]) %  DSPL_RAND_MOD_X1;
     for(m = 3; m > 0; m--)
     {
       x1[m] = x1[m-1];
       x2[m] = x2[m-1];
     }

     x[k] = (double)y/DSPL_RAND_MOD_X1;
   }

  return RES_OK;
}





/*******************************************************************************
Gaussian random numbers generator
*******************************************************************************/
int DSPL_API randn(double* x, int n, double mu, double sigma)
{
  int k, m;
  double x1[128], x2[128];
  int res;
  if(!x)
    return ERROR_PTR;

   if(n<1)
    return ERROR_SIZE;

  if(sigma < 0.0)
    return ERROR_RAND_SIGMA;

  k=0;
  while(k < n)
  {
    res = randu(x1, 128);
    if(res != RES_OK)
      goto exit_label;

    res = randu(x2, 128);
    if(res != RES_OK)
      goto exit_label;
    m = 0 ;
    while(k<n && m < 128)
    {
      x[k] = sqrt(-2.0*log(x1[m]))*cos(M_2PI*x2[m])*sigma + mu;
      k++;
      m++;
      if(k<n && m < 128)
      {
        x[k] = sqrt(-2.0*log(x1[m]))*sin(M_2PI*x2[m])*sigma + mu;
        k++;
        m++;
      }
    }
  }

  res = RES_OK;
exit_label:
  return res;
}

