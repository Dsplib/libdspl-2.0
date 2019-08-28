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
\ingroup SPEC_MATH_TRANSCEND
\fn int cheby_poly1(double* x, int n, int ord, double* y)
\brief Chebyshev polynomial of the first kind order `ord`

Function calculates Chebyshev polynomial \f$ C_ord(x)\f$ of the first kind 
order `ord` for the real vector `x` (length `n`) by recurrent equation:
\f[
C_ord(x) = 2 x C_{ord-1}(x) - C_{ord-2}(x), 
\f]
where \f$ C_0(x) = 1 \f$, \f$ C_1(x) = x\f$

\param[in] x        Pointer to the real argument vector `x`. \n
                    Vector size is `[n x 1]`. \n \n

\param[in] n        Size of vectors `x` and `y`. \n \n

\param[in] ord      Chebyshev polynomial order. \n \n

\param[out] y       Pointer to the Chebyshev polynomial values, 
                    corresponds to the argument `x`. \n
                    Vector size is `[n x 1]`. \n
                    Memory must be allocated. \n \n

\return 
`RES_OK` if Chebyshev polynomial is calculated successfully. \n
Else \ref ERROR_CODE_GROUP "code error". \n

Example:

\include cheby_poly1_test.c

 \n \n
Text files will be created in `dat` directory: \n

<pre>
cheby_poly1_ord1.txt
cheby_poly1_ord2.txt
cheby_poly1_ord3.txt
cheby_poly1_ord4.txt
</pre>

GNUPLOT package will create Chebyshev polynomials plot from saved text-files:

\image html cheby_poly1.png

GNUPLOT script is follow:
\include cheby_poly1.plt

\author Sergey Bakhurin www.dsplib.org
*******************************************************************************/
int DSPL_API cheby_poly1(double* x, int n, int ord, double* y)
{
  int k, m;
  double t[2];

  if(!x || !y)
    return ERROR_PTR;
  if(n < 1)
    return ERROR_SIZE;
  if(ord<0)
    return ERROR_POLY_ORD;
  if(ord==0)
  {
    for(k = 0; k < n; k++)
    {
      y[k] = 1.0;
    }
    return RES_OK;
  }

  if(ord==1)
  {
    memcpy(y, x, n*sizeof(double));
    return RES_OK;
  }

  for(k = 0; k < n; k++)
  {
    m = 2;
    t[1]  = x[k];
    t[0]  = 1.0;
    while(m <= ord)
    {
      y[k] = 2.0 * x[k] *t[1] - t[0];
      t[0] = t[1];
      t[1] = y[k];
      m++;
    }
  }
  return RES_OK;
}




/******************************************************************************
\ingroup SPEC_MATH_TRANSCEND
\fn int cheby_poly2(double* x, int n, int ord, double* y)
\brief Chebyshev polynomial of the second kind order `ord`

Function calculates Chebyshev polynomial \f$ U_ord(x)\f$ of the first kind 
order `ord` for the real vector `x` (length `n`) by recurrent equation:
\f[
U_ord(x) = 2 x U_{ord-1}(x) - U_{ord-2}(x), 
\f]
where \f$ U_0(x) = 1 \f$, \f$ U_1(x) = 2x\f$

\param[in] x        Pointer to the real argument vector `x`. \n
                    Vector size is `[n x 1]`. \n \n

\param[in] n        Size of vectors `x` and `y`. \n \n

\param[in] ord      Chebyshev polynomial order. \n \n

\param[out] y       Pointer to the Chebyshev polynomial values, 
                    corresponds to the argument `x`. \n
                    Vector size is `[n x 1]`. \n
                    Memory must be allocated. \n \n

\return 
`RES_OK` if Chebyshev polynomial is calculated successfully. \n
Else \ref ERROR_CODE_GROUP "code error". \n

Example:

\include cheby_poly2_test.c

 \n \n
Text files will be created in `dat` directory: \n

<pre>
cheby_poly2_ord1.txt
cheby_poly2_ord2.txt
cheby_poly2_ord3.txt
cheby_poly2_ord4.txt
</pre>

GNUPLOT package will create Chebyshev polynomials plot from saved text-files:

\image html cheby_poly2.png

GNUPLOT script is follow:
\include cheby_poly2.plt

\author Sergey Bakhurin www.dsplib.org
*******************************************************************************/
int DSPL_API cheby_poly2(double* x, int n, int ord, double* y)
{
  int k, m;
  double t[2];

  if(!x || !y)
    return ERROR_PTR;
  if(n < 1)
    return ERROR_SIZE;
  if(ord<0)
    return ERROR_POLY_ORD;
  if(ord==0)
  {
    for(k = 0; k < n; k++)
    {
      y[k] = 1.0;
    }
    return RES_OK;
  }

  if(ord==1)
  {
    for(k = 0; k < n; k++)
    {
      y[k] = 2.0*x[k];
    };
    return RES_OK;
  }

  for(k = 0; k < n; k++)
  {
    m = 2;
    t[1]  = 2.0*x[k];
    t[0]  = 1.0;
    while(m <= ord)
    {
      y[k] = 2.0 * x[k] *t[1] - t[0];
      t[0] = t[1];
      t[1] = y[k];
      m++;
    }
  }
  return RES_OK;
}

