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
\fn int ellip_sn_cmplx(complex_t* u, int n, double k, complex_t* y)
\brief  Jacobi elliptic function \f$ y = \textrm{sn}(u K(k), k)\f$ of 
complex vector argument

Function calculates Jacobi elliptic function  
\f$ y = \textrm{sn}(u K(k), k)\f$  of complex vector `u` and
elliptical modulus `k`. \n

\param[in]  u
Pointer to the argument vector \f$ u \f$. \n
Vector size is `[n x 1]`. \n
Memory must be allocated. \n \n

\param[in]  n
Size of vector `u`. \n 

\param[in]  k
Elliptical modulus \f$ k \f$. \n
Elliptical modulus is real parameter,
which values can be from  0 to 1. \n \n 

\param[out]  y
Pointer to the vector of Jacobi elliptic function
\f$ y = \textrm{sn}(u K(k), k)\f$. \n
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
\fn int  ellip_sn_cmplx(complex_t* u, int n, double k, complex_t* y)
\brief  Эллиптическая функция Якоби 
\f$ y = \textrm{sn}(u K(k), k)\f$ комплексного аргумента

Функция рассчитывает значения значения эллиптической функции 
\f$ y = \textrm{sn}(u K(k), k)\f$  для комплексного вектора `u` и 
эллиптического  модуля `k`. \n

Для расчета используется итерационный алгоритм на основе преобразования
Ландена. \n

\param[in]  u
Указатель на массив вектора переменной \f$ u \f$. \n
Размер вектора `[n x 1]`. \n
Память должна быть выделена. \n \n

\param[in]  n
Размер вектора `u`.  \n \n 

\param[in]  k
Значение эллиптического модуля \f$ k \f$. \n
Эллиптический модуль -- вещественный параметр, 
принимающий значения от 0 до 1.  \n \n 

\param[out]  y
Указатель на вектор значений эллиптической
функции \f$ y = \textrm{sn}(u K(k), k)\f$. \n
Размер вектора `[n x 1]`. \n
Память должна быть выделена. \n \n

\return
`RES_OK` Расчет произведен успешно. \n
В противном случае \ref ERROR_CODE_GROUP "код ошибки". \n

\author Бахурин Сергей www.dsplib.org
***************************************************************************** */
#endif
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

            RE(y[m]) =  RE(tmp) * t;
            IM(y[m]) = -IM(tmp) * t;

        }
    }
    return RES_OK;
}
