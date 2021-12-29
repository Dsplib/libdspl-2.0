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
* MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.    See the
* GNU General Public License for more details.
*
* You should have received a copy of the GNU General Public License
* along with Foobar.    If not, see <http://www.gnu.org/licenses/>.
*/


#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "dspl.h"



#ifdef DOXYGEN_ENGLISH
/*! ****************************************************************************
\ingroup SPEC_MATH_ELLIP_GROUP
\fn int ellip_asn_cmplx(complex_t* w, int n, double k, complex_t* u)
\brief  Inverse Jacobi elliptic function \f$ u = \textrm{sn}^{-1}(w, k)\f$ 
of complex vector argument

Function calculates inverse Jacobi elliptic function  
\f$ u = \textrm{sn}^{-1}(w, k)\f$  of complex vector `w`. \n

\param[in]  w
Pointer to the argument vector \f$ w \f$. \n
Vector size is `[n x 1]`. \n
Memory must be allocated. \n \n

\param[in]  n
Size of vector `w`. \n 

\param[in]  k
Elliptical modulus \f$ k \f$. \n
Elliptical modulus is real parameter,
which values can be from  0 to 1. \n \n 

\param[out]  u
Pointer to the vector of inverse Jacobi elliptic function
\f$ u = \textrm{sn}^{-1}(w, k)\f$. \n
Vector size is  `[n x 1]`. \n
Memory must be allocated. \n \n

\return
`RES_OK` successful exit, else \ref ERROR_CODE_GROUP "error code". \n

\author  Sergey Bakhurin  www.dsplib.org  
***************************************************************************** */
#endif
#ifdef DOXYGEN_RUSSIAN
/*! ****************************************************************************
\ingroup SPEC_MATH_ELLIP_GROUP
\fn int ellip_asn_cmplx(complex_t* w, int n, double k, complex_t* u)  
\brief  Обратная эллиптическая функция Якоби 
  \f$ u = \textrm{sn}^{-1}(w, k)\f$ комплексного аргумента

Функция рассчитывает значения значения обратной эллиптической функции 
\f$ u = \textrm{sn}^{-1}(w, k)\f$  для комплексного вектора `w`. \n

Для расчета используется итерационный алгоритм на основе преобразования
Ландена. \n
  

\param[in]  w
Указатель на массив вектора переменной \f$ w \f$. \n
Размер вектора `[n x 1]`. \n
Память должна быть выделена. \n \n

\param[in]  n
Размер вектора `w`.  \n \n 

\param[in]  k
Значение эллиптического модуля \f$ k \f$. \n
Эллиптический модуль -- вещественный параметр, 
принимающий значения от 0 до 1.  \n \n 

\param[out]  u
Указатель на вектор значений обратной эллиптической
функции \f$ u = \textrm{sn}^{-1}(w, k)\f$. \n
Размер вектора `[n x 1]`. \n
Память должна быть выделена. \n \n

\return
`RES_OK`Расчет произведен успешно. \n
В противном случае \ref ERROR_CODE_GROUP "код ошибки". \n

\author Бахурин Сергей www.dsplib.org  
***************************************************************************** */
#endif
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
            IM(tmp1) =         - CMIM(tmp0, tmp0);

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
