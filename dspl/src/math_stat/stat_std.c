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
#include <string.h>
#include "dspl.h"






#ifdef DOXYGEN_ENGLISH
/*! ****************************************************************************
\ingroup SPEC_MATH_STAT_GROUP
\fn int stat_std(double* x, int n, double* s)
\brief Calculates the standard deviation of the input vector `x`

Function calculates the the standard deviation value
\f[
s = \sqrt{\frac{1}{n-1} \sum_{i = 0}^{n-1} \big(x(i) - \mu \big)^2},
\f]
here \f$\mu\f$ - mean value of the vector `x`:
\f[
\mu = \frac{1}{n} \sum_{i = 0}^{n-1} x(i).
\f]

\param[in] x
Pointer to the real input vector `x`. \n
Vector size is `[n x 1]`. \n \n

\param[in]  n
Size of input vector `x`. \n \n

\param[out] s
Pointer to the variable which keeps vector `x` standard deviation value.\n
Memory must be allocated. \n \n

`RES_OK` if function calculates successfully,
 else \ref ERROR_CODE_GROUP "code error".

Example:
\code{.cpp}
    double a[5] = {0.0, 1.0, 2.0, 3.0, 4.0};
    double s;
    stat_std(a, 5, &s);
    printf("\n\n Standard deviation value:    %8.1f\n", s);
\endcode 
As result the variable `s` will keep value `1.5811`.

\author Sergey Bakhurin www.dsplib.org
***************************************************************************** */
#endif
#ifdef DOXYGEN_RUSSIAN
/*! ****************************************************************************
\ingroup SPEC_MATH_STAT_GROUP
\fn int stat_std(double* x, int n, double* s)
\brief Выборочная оценка стандартного отклонения вещественного вектора `x`

Функция рассчитывает оценку стандартного отклонения
\f[
s = \sqrt{\frac{1}{n-1} \sum_{i = 0}^{n-1} \big(x(i) - \mu \big)^2},
\f]
где \f$\mu\f$ - выборочная оценка математического ожидания вектора `x`:
\f[
\mu = \frac{1}{n} \sum_{i = 0}^{n-1} x(i).
\f]

\param[in] x
Указатель на вещественный вектор `x`. \n
Размер вектора `[n x 1]`. \n \n

\param[in]  n
Размер входного вектора `x`. \n \n

\param[out] s
Указатель на адрес памяти, в который сохранить
рассчитанное значение стандартного отклонения вектора `x`. \n
Память должна быть выделена. \n \n

\return
`RES_OK` если функция выполнена успешно. \n
 В противном случае \ref ERROR_CODE_GROUP "код ошибки".

Пример:
\code{.cpp}
    double a[5] = {0.0, 1.0, 2.0, 3.0, 4.0};
    double s;
    stat_std(a, 5, &s);
    printf("\n\n Standard deviation value:    %8.1f\n", s);
\endcode 
В результате в переменную `s` будет записано значение `1.5811`.
 
\author Бахурин Сергей. www.dsplib.org
***************************************************************************** */
#endif
int DSPL_API stat_std(double* x, int n, double* s)
{
    int k, err; 
    double sum, m;
    err = mean(x, n, &m);
    if(err != RES_OK)
        goto exit_label;
    sum = (x[0] - m) * (x[0] - m);
    for(k = 1; k < n; k++)
        sum += (x[k] - m) * (x[k] - m);
    (*s) = sqrt(sum / (double)(n-1));
exit_label:
    return err;
}

