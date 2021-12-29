/*
* Copyright (c) 2015-2019 Sergey Bakhurin
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

#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "dspl.h"


#ifdef DOXYGEN_ENGLISH
/*! ****************************************************************************
\ingroup SPEC_MATH_POLY_GROUP
\fn int cheby_poly2(double* x, int n, int ord, double* y)
\brief Chebyshev polynomial of the second kind order `ord`

Function calculates Chebyshev polynomial \f$ U_ord(x)\f$ of the first kind
order `ord` for the real vector `x` (length `n`) by recurrent equation:
\f[
U_{ord}(x) = 2 x U_{ord-1}(x) - U_{ord-2}(x),
\f]
where \f$ U_0(x) = 1 \f$, \f$ U_1(x) = 2x\f$

\param[in] x
Pointer to the real argument vector `x`. \n
Vector size is `[n x 1]`. \n \n

\param[in] n
Size of vectors `x` and `y`. \n \n

\param[in] ord
Chebyshev polynomial order. \n \n

\param[out] y
Pointer to the Chebyshev polynomial values, corresponds to the argument `x`. \n
Vector size is `[n x 1]`. \n
Memory must be allocated. \n \n

\return
`RES_OK` if Chebyshev polynomial is calculated successfully. \n
Else \ref ERROR_CODE_GROUP "code error". \n

Example:
\include cheby_poly2_test.c

Text files will be created in `dat` directory: \n

\verbatim
cheby_poly2_ord1.txt
cheby_poly2_ord2.txt
cheby_poly2_ord3.txt
cheby_poly2_ord4.txt
\endverbatim

GNUPLOT package will create Chebyshev polynomials plot from saved text-files:

\image html cheby_poly2.png

\author Sergey Bakhurin www.dsplib.org
***************************************************************************** */
#endif
#ifdef DOXYGEN_RUSSIAN
/*! ****************************************************************************
\ingroup SPEC_MATH_POLY_GROUP
\fn int cheby_poly2(double* x, int n, int ord, double* y)
\brief Многочлен Чебышева второго рода порядка `ord`

Функция производит расчет многочлена Чебышева второго рода \f$ U_{ord}(x)\f$
для вещественного вектора `x` длины `n`на основе рекуррентной формулы
\f[
U_{ord}(x) = 2 x U_{ord-1}(x) - U_{ord-2}(x),
\f]
где \f$ U_0(x) = 1 \f$, \f$ U_1(x) = 2x\f$

\param[in] x
Указатель на вектор `x` аргумента полинома Чебышева второго рода. \n
Размер вектора `[n x 1]`. \n \n

\param[in] n
Размер векторов `x` и `y`. \n \n

\param[in] ord
Порядок полинома Чебышева второго рода. \n \n

\param[out] y
Указатель на вектор значений полинома Чебышева,
соответствующих аргументу `x`. \n
Размер вектора `[n x 1]`. \n
Память должна быть выделена. \n \n

\return
`RES_OK`    Расчет произведен успешно. \n
В противном случае \ref ERROR_CODE_GROUP "код ошибки". \n

Пример использования функции:

\include cheby_poly2_test.c

В каталоге `dat` будут созданы текстовые файлы значений
полиномов порядка 1-4: \n

\verbatim
cheby_poly2_ord1.txt
cheby_poly2_ord2.txt
cheby_poly2_ord3.txt
cheby_poly2_ord4.txt
\endverbatim

Кроме того программа GNUPLOT произведет построение следующих графиков
по сохраненным в файлах данным:

\image html cheby_poly2.png

\author Бахурин Сергей www.dsplib.org
***************************************************************************** */
#endif
int DSPL_API cheby_poly2(double* x, int n, int ord, double* y)
{
    int k, m;
    double t[2];

    if(!x || !y)
        return ERROR_PTR;
    if(n < 1)
        return ERROR_SIZE;
    if(ord<0)
        return ERROR_POLY_ORD;
    if(ord==0)
    {
        for(k = 0; k < n; k++)
        {
            y[k] = 1.0;
        }
        return RES_OK;
    }

    if(ord==1)
    {
        for(k = 0; k < n; k++)
        {
            y[k] = 2.0*x[k];
        };
        return RES_OK;
    }

    for(k = 0; k < n; k++)
    {
        m = 2;
        t[1]    = 2.0*x[k];
        t[0]    = 1.0;
        while(m <= ord)
        {
            y[k] = 2.0 * x[k] *t[1] - t[0];
            t[0] = t[1];
            t[1] = y[k];
            m++;
        }
    }
    return RES_OK;
}

