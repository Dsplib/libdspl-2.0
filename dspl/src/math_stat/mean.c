/*
* Copyright (c) 2015-2024 Sergey Bakhurin
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
#include "dspl.h"



#ifdef DOXYGEN_ENGLISH
/*! ****************************************************************************
\ingroup SPEC_MATH_STAT_GROUP

\brief Calculates the mean of the input vector `x`

Function calculates the mean value
\f[
m = \frac{1}{n}\sum_{i = 0}^{n-1} x(i)
\f]

\param[in] x
Pointer to the real input vector `x`. \n
Vector size is `[n x 1]`. \n \n

\param[in]  n
Size of input vector `x`. \n \n

\param[out] m
Pointer to the variable which keeps vector `x` mean value.\n
Memory must be allocated. \n \n

`RES_OK` if function calculates successfully,
 else \ref ERROR_CODE_GROUP "code error".

Example:
\code{.cpp}
    double a[5] = {0.0, 1.0, 2.0, 3.0, 4.0};
    double m;
    mean(a, 5, &m);
    printf("\n\n Mean value:    %8.1f\n", m);
\endcode 
As result the variable `m` will keep value `2`.

\author Sergey Bakhurin www.dsplib.org
***************************************************************************** */
#endif
#ifdef DOXYGEN_RUSSIAN
/*! ****************************************************************************
\ingroup SPEC_MATH_STAT_GROUP

\brief Выборочная оценка математического ожидания вещественного вектора `x`

Функция рассчитывает оценку математического ожидания
\f[
m = \frac{1}{n} \sum_{i = 0}^{n-1} x(i)
\f]

\param[in] x
Указатель на вещественный вектор `x`. \n
Размер вектора `[n x 1]`. \n \n

\param[in]  n
Размер входного вектора `x`. \n \n

\param[out] m
Указатель на адрес памяти, в который сохранить
рассчитанное значение математического ожидания вектора `x`. \n
Память должна быть выделена. \n \n

\return
`RES_OK` если функция выполнена успешно. \n
 В противном случае \ref ERROR_CODE_GROUP "код ошибки".

Пример:
\code{.cpp}
    double a[5] = {0.0, 1.0, 2.0, 3.0, 4.0};
    double m;
    mean(a, 5, &m);
    printf("\n\n Mean value:    %8.1f\n", m);
\endcode 
В результате в переменную `m` будет записано значение `2`.
 
\author Бахурин Сергей. www.dsplib.org
***************************************************************************** */
#endif
int DSPL_API mean(double* x, int n, double* m)
{
    int k; 
    double sum;
    if(!x || !m)
        return ERROR_PTR;
    if(n<1)
        return ERROR_SIZE;

    sum = x[0];
    for(k = 1; k < n; k++)
        sum += x[k];

    (*m) = sum / (double)n;
    return RES_OK;
}


