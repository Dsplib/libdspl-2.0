/*
* Copyright (c) 2015-2022 Sergey Bakhurin
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



#if DOXYGEN_ENGLISH
/*! ****************************************************************************
\ingroup SPEC_MATH_TRIG_GROUP
\fn int acos_cmplx(complex_t* x, int n, complex_t *y)
\brief The inverse of the cosine function the complex vector argument `x`.

Function calculates the inverse of the cosine function as: \n
\f[
\textrm{Arccos}(x) = \frac{\pi}{2} - \textrm{Arcsin}(x) = 
\frac{\pi}{2} -j \textrm{Ln}\left( j x + \sqrt{1 - x^2} \right)
\f]


\param[in] x
Pointer to the argument vector `x`. \n
Vector size is `[n x 1]`. \n\n

\param[in] n
Input vector `x` and the inverse cosine vector `y` size. \n\n

\param[out] y
Pointer to the output complex vector `y`, 
corresponds to the input vector `x`. \n 
Vector size is `[n x 1]`. \n
Memory must be allocated. \n\n

\return
`RES_OK` if function calculated successfully. \n
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

Output is: \n
\verbatim
acos_cmplx(1.0+2.0j) = 1.144-1.529j
acos_cmplx(3.0+4.0j) = 0.937-2.306j
acos_cmplx(5.0+6.0j) = 0.880-2.749j
\endverbatim

\author Sergey Bakhurin www.dsplib.org
***************************************************************************** */
#endif
#if DOXYGEN_RUSSIAN
/*! ****************************************************************************
\ingroup SPEC_MATH_TRIG_GROUP
\fn int acos_cmplx(complex_t* x, int n, complex_t *y)
\brief  Арккосинус комплексного аргумента `x`.

Функция рассчитывает значения арккосинуса комплексного аргумента, 
заданного вектором `x` длины `n`:  \n
\f[
\textrm{Arccos}(x) = \frac{\pi}{2} - \textrm{Arcsin}(x) = 
\frac{\pi}{2} -j \textrm{Ln}\left( j x + \sqrt{1 - x^2} \right)
\f]  

\param[in]  x
Указатель на вектор аргумента комплексного арккосинуса. \n
Размер вектора `[n x 1]`.  \n \n

\param[in]  n
Размер входного и выходного векторов `x` и `y`. \n \n

\param[out] y
Указатель на вектор значений комплексного арккосинуса,
соответствующего входному вектору `x`. \n
Размер массива `[n x 1]`.  \n
Память должна быть выделена.  \n \n

\return
`RES_OK` если значение функции рассчитано успешно   .  \n
В противном случае \ref ERROR_CODE_GROUP "код ошибки": \n

\note
Функция может использоваться для расчета арккосинуса аргумента 
большего единицы, когда вещественная функция `acos` не определена.

Например при выполнении следующего кода 
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

Результатом работы будет

\verbatim
acos_cmplx(1.0+2.0j) = 1.144-1.529j
acos_cmplx(3.0+4.0j) = 0.937-2.306j
acos_cmplx(5.0+6.0j) = 0.880-2.749j
\endverbatim

\author Бахурин Сергей www.dsplib.org 
***************************************************************************** */

#endif
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
        IM(y[k]) =     - IM(y[k]);
    }
    return RES_OK;
}

