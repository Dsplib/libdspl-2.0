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
\fn int array_scale_lin(double* x, int n,
                        double xmin, double xmax, double dx,
                        double h, double* y)
\brief Vector `x` linear transformation

Function transforms values \f$x(i)\f$, \f$i = 0,1,\ldots n\f$
to the \f$y(i)\f$, accordint to equation:

\f[
y(i) = k_x x(i) + d_x, \qquad k_x =
\frac{h}{x_{\textrm{max}} - x_{\textrm{min}}}.
\f]

All values of the vector `x` between
\f$x_{\textrm{min}}\f$ and \f$x_{\textrm{max}}\f$, transforms to
the vector `y` between \f$d_x\f$ and \f$h + d_x\f$.
Parameter \f$d_x\f$ sets mean shift of the vector `y`.

This function is convenient for translating values ​​
of different dimensions. For example it can be used
to transfer the values ​​of the vector `x`
to the graph of the height of` h`, where the height can
be set in the number of pixels, in centimeters, etc.

\param[in] x
Pointer to the input vector `x`. \n
Vector size is `[n x 1]`. \n
\n

\param[in] n
Size of vector `x`. \n
\n

\param[in] xmin
Parameter \f$x_{\textrm{min}}\f$. \n
\n

\param[in] xmax
Parameter \f$x_{\textrm{min}}\f$. \n
Value `xmax` must be more than `xmin`. \n
\n

\param[in] dx
Displacement after transformation. \n
This parameter must have output vector `y`
dimensions (pixels, centimeters). \n
\n

\param[in] h
Height of vector `y` after transforming between `dx` and `h+dx`. \n
\n

\param[out] y
Pointer to the output vector `y`. \n
Vector size is `[n x 1]`. \n
Memory must be allocated. \n
\note
Pointer    `y` can be equal to `x`.
Velues of vector `x` will be rewritten in this case. \n
\n

\return
`RES_OK` if function returns successfully. \n
Else \ref ERROR_CODE_GROUP "code error".

\author Sergey Bakhurin www.dsplib.org
***************************************************************************** */
#endif

#ifdef DOXYGEN_RUSSIAN
/*! ****************************************************************************
\ingroup ARRAY_GROUP
\fn int array_scale_lin(double* x, int n,
                        double xmin, double xmax, double dx,
                        double h, double* y)
\brief Линейное растяжение вектора данных `x`
Функция производит преобразование значений \f$x(i)\f$, \f$i = 0,1,\ldots n\f$
в значения \f$y(i)\f$, в соответствии с формулой:

\f[
y(i) = k_x x(i) + d_x, \qquad k_x =
\frac{h}{x_{\textrm{max}} - x_{\textrm{min}}}.
\f]

Таким образом, все значения входного вектора `x` в диапазоне от
\f$x_{\textrm{min}}\f$ до \f$x_{\textrm{max}}\f$, линейно растягиваются в
значения вектора `y` в диапазоне от \f$d_x\f$ до \f$h + d_x\f$.
Заметим, что \f$d_x\f$ задает линейное смещение значений вектора `y`.

Данная функция удобна для перевода величин разных размерностей, в частности,
для переноса значений вектора `x` на график высоты `h`, где высота может
быть задана в количестве пикселей, в сантиметрах и т.д.

\param[in] x
Указатель на вектор входных значений `x`. \n
Размер вектора `[n x 1]`. \n
\n

\param[in] n
Размер вектора `x`. \n
\n

\param[in] xmin
Нижняя граница диапазона трансформации. \n
\n

\param[in] xmax
Верхняя граница диапазона трансформации. \n
Значение `xmax` должно быть строго больше значения `xmin`. \n
\n

\param[in] dx
Смещение после трансформации. \n
Данный параметр должен иметь размерность выходного вектора `y`. \n
\n

\param[in] h
Диапазон значений вектора `y` после трансформации от `dx` до `h+dx`. \n
\n

\param[out] y
Указатель на вектора данных после трансформации. \n
Размер вектора `[n x 1]`. \n
Память должна быть выделена. \n
\note
Указатель `y` может совпадать с `x`, в этом случае,
данные вектора `x` будут перезаписаны линейно измененными в соответствии
с формулой выше. \n
\n

\return
`RES_OK` если функция выполнена успешно. \n
 В противном случае \ref ERROR_CODE_GROUP "код ошибки".

\author
Бахурин Сергей
www.dsplib.org
***************************************************************************** */
#endif

int DSPL_API array_scale_lin(double* x,   int n,
                             double xmin, double xmax, double dx,
                             double h,    double* y)
{
    double kx;
    int k;
    if(!x)
        return ERROR_PTR;
    if(n < 1)
        return ERROR_SIZE;
    if(h<0.0)
        return ERROR_NEGATIVE;

    if(xmin >= xmax)
        return ERROR_MIN_MAX;

    kx = h / (xmax - xmin);

    for(k = 0; k < n; k++)
        y[k] = (x[k] - xmin) * kx + dx;

    return RES_OK;
}
