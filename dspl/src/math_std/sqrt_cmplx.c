/*
* Copyright (c) 2015-2023 Sergey Bakhurin
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
#include "dspl.h"



#ifdef DOXYGEN_ENGLISH
/*! ****************************************************************************
\ingroup SPEC_MATH_COMMON_GROUP

\brief Square root of the complex vector argguument `x`.

Function calculates square root value of vector `x` length `n`:    \n
\f[
y(k) = \sqrt{x(k)}, \qquad k = 0 \ldots n-1. 
\f]    


\param[in] x
Pointer to the input complex vector `x`. \n
Vector size is `[n x 1]`.    \n \n

\param[in] n
Size of input and output vectors `x` and `y`. \n \n

\param[out] y
Pointer to the square root vector `y`. \n
Vector size is `[n x 1]`. \n
Memory must be allocated. \n \n

\return `RES_OK` if function is calculated successfully.    \n
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

Result:

\verbatim
sqrt_cmplx(1.0+2.0j) = 1.272+0.786j
sqrt_cmplx(3.0+4.0j) = 2.000+1.000j
sqrt_cmplx(5.0+6.0j) = 2.531+1.185j
\endverbatim

\author Sergey Bakhurin www.dsplib.org 
***************************************************************************** */
#endif
#ifdef DOXYGEN_RUSSIAN
/*! ****************************************************************************
\ingroup SPEC_MATH_COMMON_GROUP

\brief  Квадратный корень из комплексного вектора `x` (поэлементный).

Функция рассчитывает значения квадратного корня комплексного аргумента, 
заданного вектором `x` длины `n`:  \n
\f[
y(k) = \sqrt{x(k)}, \qquad k = 0 \ldots n-1. 
\f]  

\param[in]  x
Указатель на вектор аргумента квадратного корня. \n
Размер вектора `[n x 1]`.  \n \n

\param[in]  n
Размер входного и выходного векторов `x` и `y`. \n \n

\param[out] y
Указатель на вектор значений комплексного корня,
соответствующего входному вектору `x`. \n
Размер массива `[n x 1]`.  \n
Память должна быть выделена.  \n \n

\return
`RES_OK` если значение функции рассчитано успешно   .  \n
В противном случае \ref ERROR_CODE_GROUP "код ошибки": \n

Например при выполнении следующего кода 
\code{.cpp}
    complex_t x[3] = {{1.0, 2.0}, {3.0, 4.0}, {5.0, 6.0}};
    complex_t y[3];
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

\author Бахурин Сергей www.dsplib.org 
***************************************************************************** */
#endif
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