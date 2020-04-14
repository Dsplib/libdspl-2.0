/*
* Copyright (c) 2015-2020 Sergey Bakhurin
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
\ingroup TYPES_GROUP
\fn int cmplx2re(complex_t* x, int n, double* re, double* im)
\brief  Separate complex vector to the real and image vectors

Function fills `re` and `im` vectors corresponds to real and image
parts of the input complex array `x`.  \n


\param[in]  x   Pointer to the real complex vector. \n
                Vector size is `[n x 1]`.  \n \n

\param[in]  n   Size of the input complex vector `x` and real and image
                vectors `re` and `im`. \n \n

\param[out] re  Pointer to the real part  vector. \n
                Vector size is `[n x 1]`.  \n
                Memory must be allocated.  \n \n

\param[out] im  Pointer to the image part vector. \n
                Vector size is `[n x 1]`.  \n
                Memory must be allocated.  \n \n

\return
`RES_OK` if function converts complex vector successfully.  \n
Else \ref ERROR_CODE_GROUP "code error". \n

Example: \n
\code{.cpp}
    complex_t x[3] = {{1.0, 2.0}, {3.0, 4.0}, {5.0, 6.0}};
    double  re[3], im[3];

    cmplx2re(x, 3, re, im);
\endcode

Vectors `re` and `im` will contains:

\verbatim
re[0] = 1.0; im[0] = 2.0;
re[1] = 3.0; im[1] = 4.0;
re[2] = 5.0; im[2] = 6.0;
\endverbatim

\author Sergey Bakhurin. www.dsplib.org
*******************************************************************************/
int DSPL_API cmplx2re(complex_t* x, int n, double* re, double* im)
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
\ingroup TYPES_GROUP
\fn int re2cmplx(double* x, int n, complex_t *y)
\brief  Convert real array to the complex array.

Function copies the vector `x` to the real part of vector `y`.
Image part of the vector `y` sets as zero. \n
So complex vector contains data: \n
`y[i] = x[i] + j0, here i = 0,1,2 ... n-1`


\param[in]  x   Pointer to the real vector `x`. \n
                Vector size is `[n x 1]`.  \n \n

\param[in]  n   Size of the real vector `x` and complex vector `y`. \n \n

\param[out] y   Pointer to the complex vector `y`. \n
                Vector size is `[n x 1]`.  \n
                Memory must be allocated.  \n \n


\return
`RES_OK` if function returns successfully.  \n
Else \ref ERROR_CODE_GROUP "code error": \n



Например при выполнении следующего кода
\code{.cpp}
    double x[3] = {1.0, 2.0, 3.0};
    complex_t y[3];

    re2cmplx(x, 3, y);
\endcode

Vector `y` will keep:

\verbatim
    y[0] = 1+0j;
    y[1] = 2+0j;
    y[2] = 3+0j.
\endverbatim

\author Sergey Bakhurin. www.dsplib.org
*******************************************************************************/
int DSPL_API re2cmplx(double* x, int n, complex_t* y)
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


