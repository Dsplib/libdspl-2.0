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
\fn int ellip_acd(double* w, int n, double k, double* u)  
\brief  Inverse Jacobi elliptic function \f$ u = \textrm{cd}^{-1}(w, k)\f$ 
of the real vector argument

Function calculates inverse Jacobi elliptic function  
\f$ u = \textrm{cd}^{-1}(w, k)\f$  of the real vector `w`. \n

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
\f$ u = \textrm{cd}^{-1}(w, k)\f$. \n
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
\fn int ellip_acd(double* w, int n, double k, double* u)  
\brief  Обратная эллиптическая функция Якоби 
  \f$ u = \textrm{cd}^{-1}(w, k)\f$ вещественного аргумента

Функция рассчитывает значения обратной эллиптической функции 
\f$ u = \textrm{cd}^{-1}(w, k)\f$  для вещественного вектора `w`. \n

Для расчета используется итерационный алгоритм на основе преобразования
Ландена. \n
  
\param[in]  w
Указатель на массив вектора переменной \f$ w \f$. \n
Размер вектора `[n x 1]`. \n
Память должна быть выделена. \n \n
    
\param[in]  n
Размер вектора `w`. \n 
  
\param[in]  k   Значение эллиптического модуля \f$ k \f$. 
Эллиптический модуль -- вещественный параметр, 
принимающий значения от 0 до 1. \n \n 

\param[out]  u
Указатель на вектор значений обратной эллиптической
функции \f$ u = \textrm{cd}^{-1}(w, k)\f$. \n
Размер вектора `[n x 1]`. \n
Память должна быть выделена. \n \n

\return
`RES_OK`  Расчет произведен успешно. \n
В противном случае \ref ERROR_CODE_GROUP "код ошибки". \n

\author Бахурин Сергей www.dsplib.org
***************************************************************************** */
#endif
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
            u[m] = 2.0 * u[m] / (t+t*lnd[i]);
        }
        u[m] = 2.0 * acos(u[m]) / M_PI;
    }
    return RES_OK;
}

