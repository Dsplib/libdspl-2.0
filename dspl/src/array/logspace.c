/*
* Copyright (c) 2015-2020 Sergey Bakhurin
* Digital Signal Processing Library [http://dsplib.org]
*
* This file is part of libdspl-2.0.
*
* is free software: you can redistribute it and/or modify
* it under the terms of the GNU Lesser    General Public License as published by
* the Free Software Foundation, either version 3 of the License, or
* (at your option) any later version.
*
* DSPL is distributed in the hope that it will be useful,
* but WITHOUT ANY WARRANTY; without even the implied warranty of
* MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.    See the
* GNU General Public License for more details.
*
* You should have received a copy of the GNU Lesser General Public License
* along with Foobar.    If not, see <http://www.gnu.org/licenses/>.
*/


#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "dspl.h"



#ifdef DOXYGEN_ENGLISH
/*! ****************************************************************************
\ingroup ARRAY_GROUP
\fn int logspace(double x0, double x1, int n, int type, double* x)
\brief Function fills a vector with `n` logarithmically spaced elements
between \f$10^{x_0}\f$ and \f$10^{x_1}\f$.


Function supports two kinds of filling according to `type` parameter: \n

Symmetric fill (parameter `type=DSPL_SYMMETRIC`): \n

\f$x(k) = 10^{x_0} \cdot dx^k\f$, here \f$dx = \sqrt[n-1]{10^{x_1 - x_0}}\f$,
\f$k = 0 \ldots n-1.\f$

Periodic fill (parameter `type=DSPL_PERIODIC`): \n

\f$x(k) = 10^{x_0} \cdot dx^k\f$, here \f$dx = \sqrt[n]{10^{x_1 - x_0}}\f$,
\f$k = 0 \ldots n-1.\f$ \n


\param[in] x0
Start exponent value \f$x_0\f$. \n \n

\param[in] x1
End exponent value \f$x_1\f$. \n \n

\param[in] n
Number of points `x` (size of vector `x`). \n \n

\param[in] type
Fill type: \n
`DSPL_SYMMETRIC` --- symmetric, \n
`DSPL_PERIODIC` --- periodic. \n \n

\param[in,out] x
Pointer to the output logarithmically spaced vector `x` . \n
Vector size is `[n x 1]`. \n
Memory must be allocated. \n \n

\return
`RES_OK` if function returns successfully. \n
 Else \ref ERROR_CODE_GROUP "error code".

\note
Difference between symmetric and periodic filling we can
understand from the follow examples.    \n
Example 1. Periodic fill.
\code{.cpp}
    double x[5];
    logspace(-2, 3, 5, DSPL_PERIODIC, x);
\endcode

Values in the vector `x` are:

\verbatim
0.01,    0.1,    1,    10,    100
\endverbatim

\n \n

Example 2. Symmetric fill.
\code{.cpp}
        double x[5];
        logspace(-2, 3, 5, DSPL_SYMMETRIC, x);
\endcode

Values in the vector `x` are:

\verbatim
0.01    0.178    3.162    56.234    1000
\endverbatim

\author Sergey Bakhurin www.dsplib.org
***************************************************************************** */
#endif
#ifdef DOXYGEN_RUSSIAN
/*! ****************************************************************************
\ingroup ARRAY_GROUP
\fn int logspace(double x0, double x1, int n, int type, double* x)
\brief    Функция заполняет массив значениями логарифмической шкале

Заполняет массив `x` длиной `n` значениями в диапазоне
от \f$10^{x_0}\f$ до \f$10^{x_1}\f$. \n
Функция поддерживает два типа заполнения в соответствии с параметром `type`: \n

Симметричное заполнение согласно выражению: \n

\f$x(k) = 10^{x_0} \cdot dx^k\f$, где \f$dx = \sqrt[n-1]{10^{x_1 - x_0}}\f$,
\f$k = 0 \ldots n-1.\f$

Периодическое заполнение согласно выражению:

\f$x(k) = 10^{x_0} \cdot dx^k\f$, где \f$dx = \sqrt[n]{10^{x_1 - x_0}}\f$,
\f$k = 0 \ldots n-1.\f$ \n

\param[in] x0
Начальное значение показателя \f$x_0\f$. \n \n

\param[in] x1
Конечное значение показателя \f$x_1\f$. \n \n

\param[in] n
Количество точек массива `x`. \n \n

\param[in] type
Тип заполнения: \n
`DSPL_SYMMETRIC` --- симметричное заполнение, \n
`DSPL_PERIODIC` --- периодическое заполнение. \n \n

\param[in,out] x
Указатель на вектор значений в логарифмической шкале. \n
Размер вектора `[n x 1]`. \n
Память должна быть выделена. \n \n

\return
`RES_OK` --- функция выполнена успешно.    \n
В противном случае \ref ERROR_CODE_GROUP "код ошибки".

\note
Отличие периодического и симметричного заполнения можно
понять из следующих примеров.    \n

Пример 1. Периодическое заполнение.
\code{.cpp}
    double x[5];
    logspace(-2, 3, 5, DSPL_PERIODIC, x);
\endcode
В массиве `x` будут лежать значения:
\verbatim
0.01,    0.1,    1,    10,    100
\endverbatim
\n \n

Пример 2. Симметричное заполнение.
\code{.cpp}
    double x[5];
    logspace(-2, 3, 5, DSPL_SYMMETRIC, x);
\endcode

В массиве `x` будут лежать значения:

\verbatim
0.01    0.178    3.162    56.234    1000
\endverbatim

\author Бахурин Сергей www.dsplib.org
***************************************************************************** */
#endif
int DSPL_API logspace(double x0, double x1, int n, int type, double* x)
{
    double mx, a, b;
    int k;

    if(n < 2)
        return ERROR_SIZE;
    if(!x)
        return ERROR_PTR;

    a = pow(10.0, x0);
    b = pow(10.0, x1);

    switch (type)
    {
        case DSPL_SYMMETRIC:
            mx = pow(b/a, 1.0/(double)(n-1));
            x[0] = a;
            for(k = 1; k < n; k++)
                x[k] = x[k-1] * mx;
            break;
        case DSPL_PERIODIC:
            mx = pow(b/a, 1.0/(double)n);
            x[0] = a;
            for(k = 1; k < n; k++)
                x[k] = x[k-1] * mx;
            break;
        default:
            return ERROR_SYM_TYPE;
    }
    return RES_OK;
}


