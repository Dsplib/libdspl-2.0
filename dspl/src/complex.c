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
\ingroup SPEC_MATH_TRIG_GROUP
\fn int acos_cmplx(complex_t* x, int n, complex_t *y)
\brief  The inverse of the cosine function the complex vector argument `x`

Function calculates the inverse of the cosine function as: \n

\f[
\textrm{Arccos}(x) = \frac{\pi}{2} - \textrm{Arcsin}(x) = 
\frac{\pi}{2} -j \textrm{Ln}\left( j x + \sqrt{1 - x^2} \right)
\f]  


\param[in]  x   Pointer to the argument vector `x`. \n
                Vector size is `[n x 1]`.  \n \n

\param[in]  n   Input vector `x` and the inverse cosine vector `y` size. \n \n
    

\param[out] y   Pointer to the output complex vector `y`,
                corresponds to the input vector `x`. \n
                Vector size is `[n x 1]`.  \n
                Memory must be allocated.  \n \n

\return
`RES_OK` if function calculated successfully.  \n
Else \ref ERROR_CODE_GROUP "code error". \n

Example: \n
\code{.cpp}
  complex_t x[3] = {{1.0, 2.0}, {3.0, 4.0}, {5.0, 6.0}};
  complex_t y[3];
  int k;
  
  acos_cmplx(x, 3, y);
  
  for(k = 0; k < 3; k++)
    printf("acos_cmplx(%.1f%+.1fj) = %.3f%+.3fj\n", 
             RE(x[k]), IM(x[k]), RE(y[k]), IM(y[k]));
\endcode 
 \n

Output is: \n
\verbatim
acos_cmplx(1.0+2.0j) = 1.144-1.529j
acos_cmplx(3.0+4.0j) = 0.937-2.306j
acos_cmplx(5.0+6.0j) = 0.880-2.749j
\endverbatim

\author
Sergey Bakhurin www.dsplib.org
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
\ingroup SPEC_MATH_TRIG_GROUP
\fn int asin_cmplx(complex_t* x, int n, complex_t *y)
\brief  The inverse of the sine function the complex vector argument `x`

Function calculates the inverse of the sine function as: \n

\f[
 \textrm{Arcsin}(x) = j \textrm{Ln}\left( j x + \sqrt{1 - x^2} \right)
\f]  


\param[in]  x   Pointer to the argument vector `x`. \n
                Vector size is `[n x 1]`.  \n \n

\param[in]  n   Input vector `x` and the inverse sine vector `y` size. \n \n
    

\param[out] y   Pointer to the output complex vector `y`,
                corresponds to the input vector `x`. \n
                Vector size is `[n x 1]`.  \n
                Memory must be allocated.  \n \n

\return
`RES_OK` if function calculated successfully.  \n
Else \ref ERROR_CODE_GROUP "code error". \n

Example: \n
\code{.cpp}
  complex_t x[3] = {{1.0, 2.0}, {3.0, 4.0}, {5.0, 6.0}};
  complex_t y[3];
  int k;
  
  asin_cmplx(x, 3, y);  
  for(k = 0; k < 3; k++)
    printf("asin_cmplx(%.1f%+.1fj) = %.3f%+.3fj\n", 
            RE(x[k]), IM(x[k]), RE(y[k]), IM(y[k]));

\endcode 
 \n

Output is: \n
\verbatim
asin_cmplx(1.0+2.0j) = 0.427+1.529j
asin_cmplx(3.0+4.0j) = 0.634+2.306j
asin_cmplx(5.0+6.0j) = 0.691+2.749j
\endverbatim

\author
Sergey Bakhurin www.dsplib.org
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
    IM(tmp) =     - CMIM(x[k], x[k]); // 1-x[k]^2
    sqrt_cmplx(&tmp, 1, y+k);   // sqrt(1 - x[k]^2)
    RE(y[k]) -= IM(x[k]);       // j * x[k] + sqrt(1 - x[k]^2)
    IM(y[k]) += RE(x[k]);       // j * x[k] + sqrt(1 - x[k]^2)
    log_cmplx(y+k, 1, &tmp);    // log( j * x[k] + sqrt(1 - x[k]^2) )
    RE(y[k]) =  IM(tmp);        // -j * log( j * x[k] + sqrt(1 - x[k]^2) )
    IM(y[k]) = -RE(tmp);        // -j * log( j * x[k] + sqrt(1 - x[k]^2) )
  }
  return RES_OK;
}




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
\ingroup SPEC_MATH_TRIG_GROUP
\fn int cos_cmplx(complex_t* x, int n, complex_t *y)
\brief  The cosine function the complex vector argument `x`

Function calculates the cosine function as: \n

\f[
\textrm{cos}(x) = \frac{\exp(jx) + \exp(-jx)}{2} 
\f]  


\param[in]  x   Pointer to the argument vector `x`. \n
                Vector size is `[n x 1]`.  \n \n

\param[in]  n   Input vector `x` and the cosine vector `y` size. \n \n
    

\param[out] y   Pointer to the output complex vector `y`,
                corresponds to the input vector `x`. \n
                Vector size is `[n x 1]`.  \n
                Memory must be allocated.  \n \n

\return
`RES_OK` if function calculated successfully.  \n
Else \ref ERROR_CODE_GROUP "code error". \n

Example: \n
\code{.cpp}
  complex_t x[3] = {{1.0, 2.0}, {3.0, 4.0}, {5.0, 6.0}};
  complex_t y[3];
  int k;
  
  cos_cmplx(x, 3, y);
  
  for(k = 0; k < 3; k++)
    printf("cos_cmplx(%.1f%+.1fj) = %9.3f%+9.3fj\n", 
            RE(x[k]), IM(x[k]), RE(y[k]), IM(y[k]));
  
\endcode 
 \n

Output is: \n
\verbatim
cos_cmplx(1.0+2.0j) =     2.033   -3.052j
cos_cmplx(3.0+4.0j) =   -27.035   -3.851j
cos_cmplx(5.0+6.0j) =    57.219 +193.428j
\endverbatim

\author
Sergey Bakhurin www.dsplib.org
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
\ingroup SPEC_MATH_COMMON_GROUP
\fn int log_cmplx(complex_t* x, int n, complex_t *y)
\brief  The logarithm function the complex vector argument `x`

Function calculates the logarithm function as: \n

\f[
\textrm{Ln}(x) = j \varphi + \ln(|x|), 
\f]  
here \f$\varphi\f$ - the complex number phase.

\param[in]  x   Pointer to the argument vector `x`. \n
                Vector size is `[n x 1]`.  \n \n

\param[in]  n   Input vector `x` and the logarithm vector `y` size. \n \n
    

\param[out] y   Pointer to the output complex vector `y`,
                corresponds to the input vector `x`. \n
                Vector size is `[n x 1]`.  \n
                Memory must be allocated.  \n \n

\return
`RES_OK` if function calculated successfully.  \n
Else \ref ERROR_CODE_GROUP "code error". \n

Example: \n
\code{.cpp}
  complex_t x[3] = {{1.0, 2.0}, {3.0, 4.0}, {5.0, 6.0}};
  complex_t y[3];
  int k;
  
  log_cmplx(x, 3, y);  

  for(k = 0; k < 3; k++)
    printf("log_cmplx(%.1f%+.1fj) = %.3f%+.3fj\n", 
            RE(x[k]), IM(x[k]), RE(y[k]), IM(y[k]));
\endcode 
 \n

Output is: \n
\verbatim
log_cmplx(1.0+2.0j) = 0.805+1.107j
log_cmplx(3.0+4.0j) = 1.609+0.927j
log_cmplx(5.0+6.0j) = 2.055+0.876j
\endverbatim

\author
Sergey Bakhurin www.dsplib.org
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





/******************************************************************************
\ingroup SPEC_MATH_TRIG_GROUP
\fn int sin_cmplx(complex_t* x, int n, complex_t *y)
\brief  The sine function the complex vector argument `x`

Function calculates the sine function as: \n

\f[
\textrm{cos}(x) = \frac{\exp(jx) - \exp(-jx)}{2j} 
\f]  


\param[in]  x   Pointer to the argument vector `x`. \n
                Vector size is `[n x 1]`.  \n \n

\param[in]  n   Input vector `x` and the sine vector `y` size. \n \n
    

\param[out] y   Pointer to the output complex vector `y`,
                corresponds to the input vector `x`. \n
                Vector size is `[n x 1]`.  \n
                Memory must be allocated.  \n \n

\return
`RES_OK` if function calculated successfully.  \n
Else \ref ERROR_CODE_GROUP "code error". \n

Example: \n
\code{.cpp}
  complex_t x[3] = {{1.0, 2.0}, {3.0, 4.0}, {5.0, 6.0}};
  complex_t y[3];
  int k;
  
  sin_cmplx(x, 3, y);
  
  for(k = 0; k < 3; k++)
    printf("sin_cmplx(%.1f%+.1fj) = %9.3f%+9.3fj\n", 
            RE(x[k]), IM(x[k]), RE(y[k]), IM(y[k]));
   
\endcode 
 \n

Output is: \n
\verbatim
sin_cmplx(1.0+2.0j) =     3.166   +1.960j
sin_cmplx(3.0+4.0j) =     3.854  -27.017j
sin_cmplx(5.0+6.0j) =  -193.430  +57.218j
\endverbatim

\author
Sergey Bakhurin www.dsplib.org
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
\ingroup SPEC_MATH_COMMON_GROUP
\fn int sqrt_cmplx(complex_t* x, int n, complex_t *y)
\brief Square root of the complex vector argguument `x`.

Function calculates square root value of vector `x` length `n`:  \n
\f[
y(k) = \sqrt{x(k)}, \qquad k = 0 \ldots n-1. 
\f]  


\param[in]  x   Pointer to the input complex vector `x`. \n
                Vector size is `[n x 1]`.  \n \n

\param[in]  n   Size of input and output vectors `x` and `y`. \n \n
    

\param[out] y   Pointer to the square root vector `y`. \n
                Vector size is `[n x 1]`.  \n
                Memory must be allocated.  \n \n

\return `RES_OK` if function is calculated successfully.  \n
Else \ref ERROR_CODE_GROUP "code error". \n

Example
\code{.cpp}
  complex_t x[3] = {{1.0, 2.0}, {3.0, 4.0}, {5.0, 6.0}};
  complex_t y[3]
  int k;
  
  sqrt_cmplx(x, 3, y);
  
  for(k = 0; k < 3; k++)
    printf("sqrt_cmplx(%.1f%+.1fj) = %.3f%+.3fj\n", 
            RE(x[k]), IM(x[k]), RE(y[k]), IM(y[k]));
 
 \endcode 
 \n

Результатом работы будет

\verbatim
sqrt_cmplx(1.0+2.0j) = 1.272+0.786j
sqrt_cmplx(3.0+4.0j) = 2.000+1.000j
sqrt_cmplx(5.0+6.0j) = 2.531+1.185j
\endverbatim

\author Sergey Bakhurin www.dsplib.org 
*******************************************************************************/
int DSPL_API sqrt_cmplx(complex_t* x, int n, complex_t *y)
{
  int k;
  double r, zr, at;
  complex_t t;
  if(!x || !y)
    return ERROR_PTR;
  if(n < 1)
    return ERROR_SIZE;

  for(k = 0; k < n; k++)
  {
    r = ABS(x[k]);
    if(r == 0.0)
    {
      RE(y[k]) = 0.0;
      IM(y[k]) = 0.0;
    }
    else
    {
      RE(t) = RE(x[k]) + r;
      IM(t) = IM(x[k]);
      at = ABS(t);
      if(at == 0.0)
      {
        RE(y[k]) = 0.0;
        IM(y[k]) = sqrt(r);
      }
      else
      {
        zr = 1.0 / ABS(t);
        r = sqrt(r);
        RE(y[k]) = RE(t) * zr * r;
        IM(y[k]) = IM(t) * zr * r;
      }
    }
  }
  return RES_OK;
}

