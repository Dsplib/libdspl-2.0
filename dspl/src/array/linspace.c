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
\fn int linspace(double x0, double x1, int n, int type, double* x)
\brief Function fills a vector with `n` linearly spaced elements
between `x0` and `x1`.

Function supports two kinds of filling according to `type` parameter: \n

Symmetric fill (parameter `type=DSPL_SYMMETRIC`): \n

\f$x(k) = x_0 + k \cdot dx\f$,
\f$dx = \frac{x_1 - x_0}{n-1}\f$, \f$k = 0 \ldots n-1.\f$

Periodic fill (parameter `type=DSPL_PERIODIC`): \n

\f$x(k) = x_0 + k \cdot dx\f$,
\f$dx = \frac{x_1 - x_0}{n}\f$, \f$k = 0 \ldots n-1.\f$

\param[in] x0
Start point \f$x_0\f$. \n \n

\param[in] x1
End point \f$x_1\f$. \n \n

\param[in] n
Number of points `x` (size of vector `x`). \n \n

\param[in] type
Fill type: \n
`DSPL_SYMMETRIC` --- symmetric, \n
`DSPL_PERIODIC` --- periodic. \n \n

\param[in,out] x
Pointer to the output linearly spaced vector `x`. \n
Vector size is `[n x 1]`. \n
Memory must be allocated. \n \n

\return
`RES_OK` if function returns successfully. \n
 Else \ref ERROR_CODE_GROUP "error code".

\note
Difference between symmetric and periodic filling we can
understand from the follow examples.    \n
Example 1. Periodic fill.
    double x[5];
    linspace(0, 5, 5, DSPL_PERIODIC, x);
\endcode
Values in the vector `x` are:
\verbatim
0,    1,    2,    3,    4
\endverbatim
 \n \n
Example 2. Symmetric fill.
\code{.cpp}
    double x[5];
    linspace(0, 5, 5, DSPL_SYMMETRIC, x);
\endcode
Values in the vector `x` are:
\verbatim
0,    1.25,    2.5,    3.75,    5
\endverbatim

\author Sergey Bakhurin www.dsplib.org
***************************************************************************** */
#endif
#ifdef DOXYGEN_RUSSIAN
/*! ****************************************************************************
\ingroup ARRAY_GROUP
\fn int linspace(double x0, double x1, int n, int type, double* x)
\brief Функция заполняет массив линейно-нарастающими, равноотстоящими
значениями от `x0` до `x1`

Заполняет массив `x` длиной `n` значениями в диапазоне
от \f$x_0\f$ до \f$x_1\f$. Функция поддерживает два типа заполнения
в соответствии с параметром `type`: \n

Симметричное заполнение согласно выражению (параметр `type=DSPL_SYMMETRIC`): \n

\f$x(k) = x_0 + k \cdot dx\f$,
\f$dx = \frac{x_1 - x_0}{n-1}\f$, \f$k = 0 \ldots n-1.\f$

Периодическое заполнение (параметр `type=DSPL_PERIODIC`) согласно выражению: \n

\f$x(k) = x_0 + k \cdot dx\f$,
\f$dx = \frac{x_1 - x_0}{n}\f$, \f$k = 0 \ldots n-1.\f$

\param[in] x0
Начальное показателя \f$x_0\f$. \n \n

\param[in] x1
Конечное значение \f$x_1\f$. \n \n

\param[in] n
Количество точек массива `x`. \n \n

\param[in] type
Тип заполнения: \n

`DSPL_SYMMETRIC` --- симметричное заполнение, \n
`DSPL_PERIODIC` --- периодическое заполнение. \n \n

\param[in,out] x
Указатель на вектор равноотстоящих значений . \n
Размер вектора `[n x 1]`. \n
Память должна быть выделена. \n \n

\return
`RES_OK` --- функция выполнена успешно.    \n
В противном случае \ref ERROR_CODE_GROUP "код ошибки". \n \n

\note
Отличие периодического и симметричного заполнения можно
понять из следующих примеров.    \n
Пример 1. Периодическое заполнение.
\code{.cpp}
    double x[5];
    linspace(0, 5, 5, DSPL_PERIODIC, x);
\endcode
В массиве `x` будут лежать значения:
\verbatim
0,    1,    2,    3,    4
\endverbatim
 \n \n
Пример 2. Симметричное заполнение.
\code{.cpp}
        double x[5];
        linspace(0, 5, 5, DSPL_SYMMETRIC, x);
\endcode
В массиве `x` будут лежать значения:
\verbatim
0,    1.25,    2.5,    3.75,    5
\endverbatim

\author
Бахурин Сергей
www.dsplib.org
***************************************************************************** */
#endif
int DSPL_API linspace(double x0, double x1, int n, int type, double* x)
{
    double dx;
    int k;

    if(n < 2)
        return ERROR_SIZE;
    if(!x)
    return ERROR_PTR;

    switch (type)
    {
        case DSPL_SYMMETRIC:
            dx = (x1 - x0)/(double)(n-1);
            x[0] = x0;
            for(k = 1; k < n; k++)
                x[k] = x[k-1] + dx;
            break;
        case DSPL_PERIODIC:
            dx = (x1 - x0)/(double)n;
            x[0] = x0;
            for(k = 1; k < n; k++)
                x[k] = x[k-1] + dx;
            break;
        default:
            return ERROR_SYM_TYPE;
    }
    return RES_OK;
}


