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
\fn int decimate(double* x, int n, int d, double* y, int* cnt)
\brief
Real vector decimation

Function `d` times decimates real vector `x`. \n
Output vector `y` keeps values corresponds to:
`y(k) = x(k*d), k = 0...n/d-1` \n

\param[in] x
Pointer to the input real vector `x`. \n
Vector `x` size is `[n x 1]`. \n \n

\param[in] n
Size of input vector `x`. \n \n

\param[in] d
Decimation coefficient. \n
Each d-th vector will be copy from vector `x` to the
output vector `y`. \n \n

\param[out] y
Pointer to the output decimated vector `y`. \n
Output vector size is `[n/d x 1]` will be copy
to the address `cnt`. \n

\param[out] cnt
Address which will keep decimated vector `y` size. \n
Pointer can be `NULL`, vector `y` will not return
in this case. \n \n

\return
`RES_OK` if function calculated successfully. \n
Else \ref ERROR_CODE_GROUP "code error".

Two-times decimation example:
\code{.cpp}
double x[10] = {0.0, 1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0, 9.0};
double y[5];
int d = 2;
int cnt;

decimate(x, 10, d, y, &cnt);
\endcode
As result variable `cnt` will be written value 5 and
vector `y` will keep    array:
\verbatim
c = [0.0, 2.0, 4.0, 6.0, 8.0]
\endverbatim

\author Sergey Bakhurin www.dsplib.org
***************************************************************************** */
#endif
#ifdef DOXYGEN_RUSSIAN
/*! ****************************************************************************
\ingroup ARRAY_GROUP
\fn int decimate(double* x, int n, int d, double* y, int* cnt)
\brief Децимация вещественного вектора данных

Функция производит децимацию вещественного вектора `x` в `d` раз. \n
В результате выходной вектор `y` содержит значения:
`y(k) = x(k*d), k = 0...n/d-1` \n

\param[in] x
Указатель на вектор входных данных `x`. \n
Размер вектора `[n x 1]`. \n \n

\param[in] n
Размер входного вектора `x`. \n \n

\param[in] d
Коэффициент децимации. \n
В результате децимации из вектора `x` будет взять каждый
d-й элемент. \n \n

\param[out] y
Указатель на децимированный вектор `y`. \n
Размер выходного вектора равен `[n/d x 1]`
будет сохранен по адресу `cnt`. \n
Память должна быть выделена. \n \n

\param[out] cnt
Указатель переменную, в которую будет сохранен
размер выходного вектора после децимации. \n
Указатель может быть `NULL`, в этом случае
размер вектора `y` не возвращается. \n \n

\return
`RES_OK` если функция выполнена успешно. \n
 В противном случае \ref ERROR_CODE_GROUP "код ошибки".

Пример децимации вещественного массива данных в 2 раза:
\code{.cpp}
double x[10] = {0.0, 1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0, 9.0};
double y[5];
int d = 2;
int cnt;
decimate(x, 10, d, y, &cnt);
\endcode
В результате в переменную `cnt` будет записан размер 5,
а вектор `y` будет хранить массив данных:
\verbatim
c = [0.0, 2.0, 4.0, 6.0, 8.0]
\endverbatim

\author
Бахурин Сергей
www.dsplib.org
**************************************************************************** */
#endif
int DSPL_API decimate(double* x, int n, int d, double* y, int* cnt)
{
    int k = 0, i = 0;
    if(!x || !y)
        return ERROR_PTR;
    if(n < 1)
        return ERROR_SIZE;
    if(d < 1)
        return ERROR_NEGATIVE;

    k = i = 0;
    while(k + d <= n)
    {
        y[i] = x[k];
        k+=d;
        i++;
    }
    if(cnt)
        *cnt = i;

    return RES_OK;
}

