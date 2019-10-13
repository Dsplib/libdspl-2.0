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
#include <string.h>
#include <math.h>
#include "dspl.h"


/*****************************************************************************
\ingroup SPEC_MATH_ELLIP_GROUP
\fn int ellip_acd(double* w, int n, double k, double* u)  
\brief  Inverse Jacobi elliptic function  
\f$ u = \textrm{cd}^{-1}(w, k)\f$ of real vector argument

Function calculates inverse Jacobi elliptic function  
\f$ u = \textrm{cd}^{-1}(w, k)\f$  of real vector `w`. \n

\param[in]  w   Pointer to the argument vector \f$ w \f$. \n
                Vector size is `[n x 1]`. \n
                Memory must be allocated. \n \n
    
\param[in]  n   Size of vector `w`. \n 
  
\param[in]  k   Elliptical modulus \f$ k \f$. \n
                Elliptical modulus is real parameter,
                which values can be from  0 to 1. \n \n 
    

\param[out]  u  Pointer to the vector of inverse Jacobi elliptic function
                \f$ u = \textrm{cd}^{-1}(w, k)\f$. \n
                Vector size is  `[n x 1]`. \n
                Memory must be allocated. \n \n


\return
`RES_OK` successful exit, else \ref ERROR_CODE_GROUP "error code". \n

\author  Sergey Bakhurin  www.dsplib.org  
 ******************************************************************************/
int DSPL_API ellip_acd(double* w, int n, double k, double* u)
{
  double lnd[ELLIP_ITER], t;
  int i, m;

  if(!u || !w)
    return ERROR_PTR;
  if(n<1)
    return ERROR_SIZE;
  if(k < 0.0 || k>= 1.0)
    return ERROR_ELLIP_MODULE;

  ellip_landen(k,ELLIP_ITER, lnd);


  for(m = 0; m < n; m++)
  {
    u[m] = w[m];
    for(i = 1; i < ELLIP_ITER; i++)
    {
      t = lnd[i-1]*u[m];
      t *= t;
      t = 1.0 + sqrt(1.0 - t);
      u[m]  = 2.0 * u[m] / (t+t*lnd[i]);
    }
    u[m] = 2.0 * acos(u[m]) / M_PI;
  }
  return RES_OK;
}





/*****************************************************************************
\ingroup SPEC_MATH_ELLIP_GROUP
\fn int ellip_acd_cmplx(complex_t* w, int n, double k, complex_t* u)
\brief  Inverse Jacobi elliptic function  
\f$ u = \textrm{cd}^{-1}(w, k)\f$ of complex vector argument

Function calculates inverse Jacobi elliptic function  
\f$ u = \textrm{cd}^{-1}(w, k)\f$  of complex vector `w`. \n

\param[in]  w   Pointer to the argument vector \f$ w \f$. \n
                Vector size is `[n x 1]`. \n
                Memory must be allocated. \n \n
    
\param[in]  n   Size of vector `w`. \n 
  
\param[in]  k   Elliptical modulus \f$ k \f$. \n
                Elliptical modulus is real parameter,
                which values can be from  0 to 1. \n \n 
    

\param[out]  u  Pointer to the vector of inverse Jacobi elliptic function
                \f$ u = \textrm{cd}^{-1}(w, k)\f$. \n
                Vector size is  `[n x 1]`. \n
                Memory must be allocated. \n \n


\return
`RES_OK` successful exit, else \ref ERROR_CODE_GROUP "error code". \n

\author  Sergey Bakhurin  www.dsplib.org  
 ******************************************************************************/
int DSPL_API ellip_acd_cmplx(complex_t* w, int n, double k, complex_t* u)
{
  double lnd[ELLIP_ITER], t;
  complex_t tmp0, tmp1;
  int i, m;

  if(!u || !w)
    return ERROR_PTR;
  if(n<1)
    return ERROR_SIZE;
  if(k < 0.0 || k>= 1.0)
    return ERROR_ELLIP_MODULE;

  ellip_landen(k,ELLIP_ITER, lnd);


  for(m = 0; m < n; m++)
  {
    RE(u[m]) = RE(w[m]);
    IM(u[m]) = IM(w[m]);
    for(i = 1; i < ELLIP_ITER; i++)
    {
      RE(tmp0) = lnd[i-1]*RE(u[m]);
      IM(tmp0) = lnd[i-1]*IM(u[m]);
      RE(tmp1) = 1.0 - CMRE(tmp0, tmp0);
      IM(tmp1) =     - CMIM(tmp0, tmp0);

      sqrt_cmplx(&tmp1, 1, &tmp0);
      RE(tmp0) += 1.0;

      RE(tmp1) = RE(tmp0) * (1.0 + lnd[i]);
      IM(tmp1) = IM(tmp0) * (1.0 + lnd[i]);

      t = 2.0 / ABSSQR(tmp1);

      RE(tmp0) = t * CMCONJRE(u[m], tmp1);
      IM(tmp0) = t * CMCONJIM(u[m], tmp1);

      RE(u[m]) = RE(tmp0);
      IM(u[m]) = IM(tmp0);

    }
    acos_cmplx(&tmp0, 1, u+m);
    t = 2.0 / M_PI;
    RE(u[m]) *= t;
    IM(u[m]) *= t;
  }
  return RES_OK;
}





/*****************************************************************************
\ingroup SPEC_MATH_ELLIP_GROUP
\fn int ellip_asn(double* w, int n, double k, double* u)  
\brief  Inverse Jacobi elliptic function  
\f$ u = \textrm{sn}^{-1}(w, k)\f$ of real vector argument

Function calculates inverse Jacobi elliptic function  
\f$ u = \textrm{sn}^{-1}(w, k)\f$  of real vector `w`. \n

\param[in]  w   Pointer to the argument vector \f$ w \f$. \n
                Vector size is `[n x 1]`. \n
                Memory must be allocated. \n \n
    
\param[in]  n   Size of vector `w`. \n 
  
\param[in]  k   Elliptical modulus \f$ k \f$. \n
                Elliptical modulus is real parameter,
                which values can be from  0 to 1. \n \n 
    

\param[out]  u  Pointer to the vector of inverse Jacobi elliptic function
                \f$ u = \textrm{sn}^{-1}(w, k)\f$. \n
                Vector size is  `[n x 1]`. \n
                Memory must be allocated. \n \n


\return
`RES_OK` successful exit, else \ref ERROR_CODE_GROUP "error code". \n

\author  Sergey Bakhurin  www.dsplib.org  
 ******************************************************************************/
int DSPL_API ellip_asn(double* w, int n, double k, double* u)
{
  double lnd[ELLIP_ITER], t;
  int i, m;

  if(!u || !w)
    return ERROR_PTR;
  if(n<1)
    return ERROR_SIZE;
  if(k < 0.0 || k>= 1.0)
    return ERROR_ELLIP_MODULE;

  ellip_landen(k,ELLIP_ITER, lnd);


  for(m = 0; m < n; m++)
  {
    u[m] = w[m];
    for(i = 1; i < ELLIP_ITER; i++)
    {
      t = lnd[i-1]*u[m];
      t *= t;
      t = 1.0 + sqrt(1.0 - t);
      u[m]  = 2.0 * u[m] / (t+t*lnd[i]);
    }
    u[m] = 2.0 * asin(u[m]) / M_PI;
  }
  return RES_OK;
}





/*****************************************************************************
\ingroup SPEC_MATH_ELLIP_GROUP
\fn int ellip_asn_cmplx(complex_t* w, int n, double k, complex_t* u)
\brief  Inverse Jacobi elliptic function  
\f$ u = \textrm{sn}^{-1}(w, k)\f$ of complex vector argument

Function calculates inverse Jacobi elliptic function  
\f$ u = \textrm{sn}^{-1}(w, k)\f$  of complex vector `w`. \n

\param[in]  w   Pointer to the argument vector \f$ w \f$. \n
                Vector size is `[n x 1]`. \n
                Memory must be allocated. \n \n
    
\param[in]  n   Size of vector `w`. \n 
  
\param[in]  k   Elliptical modulus \f$ k \f$. \n
                Elliptical modulus is real parameter,
                which values can be from  0 to 1. \n \n 
    

\param[out]  u  Pointer to the vector of inverse Jacobi elliptic function
                \f$ u = \textrm{sn}^{-1}(w, k)\f$. \n
                Vector size is  `[n x 1]`. \n
                Memory must be allocated. \n \n


\return
`RES_OK` successful exit, else \ref ERROR_CODE_GROUP "error code". \n

\author  Sergey Bakhurin  www.dsplib.org  
 ******************************************************************************/
int DSPL_API ellip_asn_cmplx(complex_t* w, int n, double k, complex_t* u)
{
  double lnd[ELLIP_ITER], t;
  complex_t tmp0, tmp1;
  int i, m;

  if(!u || !w)
    return ERROR_PTR;
  if(n<1)
    return ERROR_SIZE;
  if(k < 0.0 || k>= 1.0)
    return ERROR_ELLIP_MODULE;

  ellip_landen(k,ELLIP_ITER, lnd);


  for(m = 0; m < n; m++)
  {
    RE(u[m]) = RE(w[m]);
    IM(u[m]) = IM(w[m]);
    for(i = 1; i < ELLIP_ITER; i++)
    {
      RE(tmp0) = lnd[i-1]*RE(u[m]);
      IM(tmp0) = lnd[i-1]*IM(u[m]);
      RE(tmp1) = 1.0 - CMRE(tmp0, tmp0);
      IM(tmp1) =     - CMIM(tmp0, tmp0);

      sqrt_cmplx(&tmp1, 1, &tmp0);
      RE(tmp0) += 1.0;

      RE(tmp1) = RE(tmp0) * (1.0 + lnd[i]);
      IM(tmp1) = IM(tmp0) * (1.0 + lnd[i]);

      t = 2.0 / ABSSQR(tmp1);

      RE(tmp0) = t * CMCONJRE(u[m], tmp1);
      IM(tmp0) = t * CMCONJIM(u[m], tmp1);

      RE(u[m]) = RE(tmp0);
      IM(u[m]) = IM(tmp0);

    }
    asin_cmplx(&tmp0, 1, u+m);
    t = 2.0 / M_PI;
    RE(u[m]) *= t;
    IM(u[m]) *= t;
  }
  return RES_OK;
}




/*****************************************************************************
\ingroup SPEC_MATH_ELLIP_GROUP
int ellip_cd(double* u, int n, double k, double* y)
\brief  Jacobi elliptic function  
\f$ y = \textrm{cd}(u K(k), k)\f$ of real vector argument

Function calculates Jacobi elliptic function  
\f$ y = \textrm{cd}(u K(k), k)\f$  of real vector `u` and
elliptical modulus `k`. \n

\param[in]  u   Pointer to the argument vector \f$ u \f$. \n
                Vector size is `[n x 1]`. \n
                Memory must be allocated. \n \n
    
\param[in]  n   Size of vector `u`. \n 
  
\param[in]  k   Elliptical modulus \f$ k \f$. \n
                Elliptical modulus is real parameter,
                which values can be from  0 to 1. \n \n 
    

\param[out]  y  Pointer to the vector of Jacobi elliptic function
                \f$ y = \textrm{cd}(u K(k), k)\f$. \n
                Vector size is  `[n x 1]`. \n
                Memory must be allocated. \n \n


\return
`RES_OK` successful exit, else \ref ERROR_CODE_GROUP "error code". \n

\author  Sergey Bakhurin  www.dsplib.org  
 ******************************************************************************/
int DSPL_API ellip_cd(double* u, int n, double k, double* y)
{
  double lnd[ELLIP_ITER];
  int i, m;

  if(!u || !y)
    return ERROR_PTR;
  if(n<1)
    return ERROR_SIZE;
  if(k < 0.0 || k>= 1.0)
    return ERROR_ELLIP_MODULE;

  ellip_landen(k,ELLIP_ITER, lnd);


  for(m = 0; m < n; m++)
  {
    y[m] = cos(u[m] * M_PI * 0.5);
    for(i = ELLIP_ITER-1; i>0; i--)
    {
      y[m] = (1.0 + lnd[i]) / (1.0 / y[m] + lnd[i]*y[m]);
    }
  }
  return RES_OK;
}







/*****************************************************************************
\ingroup SPEC_MATH_ELLIP_GROUP
int ellip_cd_cmplx(complex_t* u, int n, double k, complex_t* y)
\brief  Jacobi elliptic function  
\f$ y = \textrm{cd}(u K(k), k)\f$ of complex vector argument

Function calculates Jacobi elliptic function  
\f$ y = \textrm{cd}(u K(k), k)\f$  of complex vector `u` and
elliptical modulus `k`. \n

\param[in]  u   Pointer to the argument vector \f$ u \f$. \n
                Vector size is `[n x 1]`. \n
                Memory must be allocated. \n \n
    
\param[in]  n   Size of vector `u`. \n 
  
\param[in]  k   Elliptical modulus \f$ k \f$. \n
                Elliptical modulus is real parameter,
                which values can be from  0 to 1. \n \n 
    

\param[out]  y  Pointer to the vector of Jacobi elliptic function
                \f$ y = \textrm{cd}(u K(k), k)\f$. \n
                Vector size is  `[n x 1]`. \n
                Memory must be allocated. \n \n


\return
`RES_OK` successful exit, else \ref ERROR_CODE_GROUP "error code". \n

\author  Sergey Bakhurin  www.dsplib.org  
 ******************************************************************************/
int DSPL_API ellip_cd_cmplx(complex_t* u, int n, double k, complex_t* y)
{
  double lnd[ELLIP_ITER], t;
  int i, m;
  complex_t tmp;

  if(!u || !y)
    return ERROR_PTR;
  if(n<1)
    return ERROR_SIZE;
  if(k < 0.0 || k>= 1.0)
    return ERROR_ELLIP_MODULE;

  ellip_landen(k,ELLIP_ITER, lnd);


  for(m = 0; m < n; m++)
  {
    RE(tmp) = RE(u[m]) * M_PI * 0.5;
    IM(tmp) = IM(u[m]) * M_PI * 0.5;

    cos_cmplx(&tmp, 1, y+m);

    for(i = ELLIP_ITER-1; i>0; i--)
    {
      t = 1.0 / ABSSQR(y[m]);

      RE(tmp) =  RE(y[m]) * t + RE(y[m]) * lnd[i];
      IM(tmp) = -IM(y[m]) * t + IM(y[m]) * lnd[i];

      t = (1.0 + lnd[i]) / ABSSQR(tmp);

      RE(y[m]) =   RE(tmp) * t;
      IM(y[m]) =  -IM(tmp) * t;

    }
  }
  return RES_OK;
}






/*****************************************************************************
\ingroup SPEC_MATH_ELLIP_GROUP
\fn int  ellip_landen(double k, int n, double* y)
\brief  Function calculates complete elliptical integral 
coefficients  \f$ k_i \f$ 

Complete elliptical integral \f$ K(k) \f$ can be described as:

\f[
K(k) = \frac{\pi}{2} \prod_{i = 1}^{\infty}(1+k_i),
\f]

here \f$ k_i \f$ -- coefficients which calculated 
iterative from \f$ k_0 = k\f$: 

\f[
k_i = 
\left( 
\frac{k_{i-1}}
{
1+\sqrt{1-k_{i-1}^2}
}
\right)^2
\f]

This function calculates `n` fist coefficients \f$ k_i \f$, which can
be used for Complete elliptical integral.
  

\param[in]  k   Elliptical modulus \f$ k \f$. \n
                Elliptical modulus is real parameter,
                which values can be from  0 to 1. \n \n

    
\param[in]  n   Number of \f$ k_i \f$ which need to calculate.\n 
                Parameter `n` is size of output vector `y`.\n 

\param[out]  y  pointer to the real vector which keep \f$ k_i \f$. \n
                Vector size is `[n x 1]`. \n
                Memory must be allocated. \n \n


\return
  `RES_OK` -- successful exit, else \ref ERROR_CODE_GROUP "error code". \n
      
Example:

\include ellip_landen_test.c

Result:

\verbatim
 i        k[i]

 1    4.625e-01
 2    6.009e-02
 3    9.042e-04
 4    2.044e-07
 5    1.044e-14
 6    2.727e-29
 7    1.859e-58
 8   8.640e-117
 9   1.866e-233
10    0.000e+00
11    0.000e+00
12    0.000e+00
13    0.000e+00
\endverbatim

\note  Complete elliptical integral converges enough fast
 if modulus \f$ k<1 \f$. There are 10 to 20 coefficients \f$ k_i \f$ ​​
 are sufficient for practical applications
 to ensure  complete elliptic integral precision within EPS.

\author Sergey Bakhurin www.dsplib.org  
 ******************************************************************************/
int DSPL_API ellip_landen(double k, int n, double* y)
{
  int i;
  y[0] = k;

  if(!y)
    return ERROR_PTR;
  if(n < 1)
    return ERROR_SIZE;
  if(k < 0.0 || k>= 1.0)
    return ERROR_ELLIP_MODULE;

  for(i = 1; i < n; i++)
  {
    y[i] = y[i-1] / (1.0 + sqrt(1.0 - y[i-1] * y[i-1]));
    y[i] *= y[i];
  }

  return RES_OK;
}





/*****************************************************************************
 * Elliptic modular equation
 ******************************************************************************/
int DSPL_API ellip_modulareq(double rp, double rs, int ord, double *k)
{
  double ep, es, ke, kp, t, sn = 0.0;
  int i, L, r;

  if(rp < 0 || rp == 0)
    return ERROR_FILTER_RP;
  if(rs < 0 || rs == 0)
    return ERROR_FILTER_RS;
  if(ord < 1)
    return ERROR_FILTER_ORD;
  if(!k)
    return ERROR_PTR;


  ep = sqrt(pow(10.0, rp*0.1)-1.0);
  es = sqrt(pow(10.0, rs*0.1)-1.0);

  ke = ep/es;

  ke = sqrt(1.0 - ke*ke);

  r = ord % 2;
  L = (ord-r)/2;

  kp = 1.0;
  for(i = 0; i < L; i++)
  {
    t = (double)(2*i+1) / (double)ord;
    ellip_sn(&t, 1, ke, &sn);
    sn*=sn;
    kp *= sn*sn;
  }

  kp *= pow(ke, (double)ord);
  *k = sqrt(1.0 - kp*kp);

  return RES_OK;

}



/*****************************************************************************
 * Elliptic rational function
 ******************************************************************************/
int DSPL_API ellip_rat(double* w, int n, int ord, double k, double* u)
{
  double t, xi, w2, xi2, k2;
  int i, m, r, L;

  if(!u || !w)
    return ERROR_PTR;
  if(n<1)
    return ERROR_SIZE;
  if(k < 0.0 || k>= 1.0)
    return ERROR_ELLIP_MODULE;

  r = ord%2;
  L = (ord-r)/2;

  if(r)
    memcpy(u, w, n*sizeof(double));
  else
  {
    for(m = 0; m < n; m++)
    {
      u[m] = 1.0;
    }
  }

  k2 = k*k;
  for(i = 0; i < L; i++)
  {
    t = (double)(2*i+1) / (double)ord;
    ellip_cd(&t, 1, k, &xi);
    xi2 = xi*xi;
    for(m = 0; m < n; m++)
    {
      w2 = w[m]*w[m];
      u[m] *= (w2 - xi2) / (1.0 - w2 * k2 * xi2);
      u[m] *= (1.0 - k2*xi2) / (1.0 - xi2);
    }
  }
  return RES_OK;
}






/*****************************************************************************
\ingroup SPEC_MATH_ELLIP_GROUP
int ellip_sn(double* u, int n, double k, double* y)
\brief  Jacobi elliptic function  
\f$ y = \textrm{sn}(u K(k), k)\f$ of real vector argument

Function calculates Jacobi elliptic function  
\f$ y = \textrm{sn}(u K(k), k)\f$  of real vector `u` and
elliptical modulus `k`. \n

\param[in]  u   Pointer to the argument vector \f$ u \f$. \n
                Vector size is `[n x 1]`. \n
                Memory must be allocated. \n \n
    
\param[in]  n   Size of vector `u`. \n 
  
\param[in]  k   Elliptical modulus \f$ k \f$. \n
                Elliptical modulus is real parameter,
                which values can be from  0 to 1. \n \n 
    

\param[out]  y  Pointer to the vector of Jacobi elliptic function
                \f$ y = \textrm{sn}(u K(k), k)\f$. \n
                Vector size is  `[n x 1]`. \n
                Memory must be allocated. \n \n


\return
`RES_OK` successful exit, else \ref ERROR_CODE_GROUP "error code". \n

\author  Sergey Bakhurin  www.dsplib.org  
 ******************************************************************************/
int DSPL_API ellip_sn(double* u, int n, double k, double* y)
{
  double lnd[ELLIP_ITER];
  int i, m;

  if(!u || !y)
    return ERROR_PTR;
  if(n<1)
    return ERROR_SIZE;
  if(k < 0.0 || k>= 1.0)
    return ERROR_ELLIP_MODULE;

  ellip_landen(k,ELLIP_ITER, lnd);


  for(m = 0; m < n; m++)
  {
    y[m] = sin(u[m] * M_PI * 0.5);
    for(i = ELLIP_ITER-1; i>0; i--)
    {
      y[m] = (1.0 + lnd[i]) / (1.0 / y[m] + lnd[i]*y[m]);
    }
  }
  return RES_OK;
}




/*****************************************************************************
\ingroup SPEC_MATH_ELLIP_GROUP
int ellip_sn_cmplx(complex_t* u, int n, double k, complex_t* y)
\brief  Jacobi elliptic function  
\f$ y = \textrm{sn}(u K(k), k)\f$ of complex vector argument

Function calculates Jacobi elliptic function  
\f$ y = \textrm{sn}(u K(k), k)\f$  of complex vector `u` and
elliptical modulus `k`. \n

\param[in]  u   Pointer to the argument vector \f$ u \f$. \n
                Vector size is `[n x 1]`. \n
                Memory must be allocated. \n \n
    
\param[in]  n   Size of vector `u`. \n 
  
\param[in]  k   Elliptical modulus \f$ k \f$. \n
                Elliptical modulus is real parameter,
                which values can be from  0 to 1. \n \n 
    

\param[out]  y  Pointer to the vector of Jacobi elliptic function
                \f$ y = \textrm{sn}(u K(k), k)\f$. \n
                Vector size is  `[n x 1]`. \n
                Memory must be allocated. \n \n


\return
`RES_OK` successful exit, else \ref ERROR_CODE_GROUP "error code". \n

\author  Sergey Bakhurin  www.dsplib.org  
 ******************************************************************************/
int DSPL_API ellip_sn_cmplx(complex_t* u, int n, double k, complex_t* y)
{
  double lnd[ELLIP_ITER], t;
  int i, m;
  complex_t tmp;

  if(!u || !y)
    return ERROR_PTR;
  if(n<1)
    return ERROR_SIZE;
  if(k < 0.0 || k>= 1.0)
    return ERROR_ELLIP_MODULE;

  ellip_landen(k,ELLIP_ITER, lnd);


  for(m = 0; m < n; m++)
  {
    RE(tmp) = RE(u[m]) * M_PI * 0.5;
    IM(tmp) = IM(u[m]) * M_PI * 0.5;

    sin_cmplx(&tmp, 1, y+m);

    for(i = ELLIP_ITER-1; i>0; i--)
    {
      t = 1.0 / ABSSQR(y[m]);

      RE(tmp) =  RE(y[m]) * t + RE(y[m]) * lnd[i];
      IM(tmp) = -IM(y[m]) * t + IM(y[m]) * lnd[i];

      t = (1.0 + lnd[i]) / ABSSQR(tmp);

      RE(y[m]) =   RE(tmp) * t;
      IM(y[m]) =  -IM(tmp) * t;

    }
  }
  return RES_OK;
}

