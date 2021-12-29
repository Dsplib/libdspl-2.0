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

#endif
#ifdef DOXYGEN_RUSSIAN
/*! ****************************************************************************
\ingroup SPEC_MATH_POLY_GROUP
\fn int polyval_cmplx(complex_t* a, int ord, complex_t* x, int n, complex_t* y)
\brief Расчет комплексного полинома

Функция рассчитывает полином \f$P_N(x)\f$ \f$N\f$-го порядка 
комплексного аргумента, заданного вектором `x`. \n

\f[
P_N(x) = a_0 + a_1 \cdot x + a_2 \cdot x^2 + a_3 \cdot x^3 + ... a_N \cdot x^N.
\f]

Для расчета используется формула Горнера: \n
\f[
  P_N(x) = a_0 + x \cdot (a_1 + x \cdot (a_2 + \cdot 
  ( \ldots x \cdot (a_{N-1} + x\cdot a_N) \ldots ))) 
\f]

\param[in]  a
Указатель на вектор комплексных коэффициентов полинома. \n
Размер вектора `[ord+1 x 1]`. \n
Коэффициент `a[0]` соответствует коэффициенту полинома \f$a_0\f$. \n \n

\param[in]  ord
Порядок полинома \f$N\f$.  \n \n

\param[in]  x
Указатель на вектор аргумента полинома. \n
Размер вектора `[n x 1]`. \n
Значения полинома будут расчитаны для всех значений аргумента вектора `x`. \n \n

\param[in]  n
Размер вектора агрумента полинома. \n\n

\param[out]  y
Указатель вектор значения полинома для аргумента `x`. \n
Размер вектора `[n x 1]`. \n
Память должна быть выделена. \n \n

\return
`RES_OK` --- полином расчитан успешно. \n
В противном случае \ref ERROR_CODE_GROUP "код ошибки".

\author Бахурин Сергей. www.dsplib.org
***************************************************************************** */
#endif
int DSPL_API polyval_cmplx(complex_t* a, int ord,
                           complex_t* x, int n, complex_t* y)
{
    int k, m;
    complex_t t;

    if(!a || !x || !y)
        return ERROR_PTR;
    if(ord<0)
        return ERROR_POLY_ORD;
    if(n<1)
        return ERROR_SIZE;

    for(k = 0; k < n; k++)
    {
        RE(y[k]) = RE(a[ord]);
        IM(y[k]) = IM(a[ord]);
        for(m = ord-1; m>-1; m--)
        {
            RE(t) = CMRE(y[k], x[k]);
            IM(t) = CMIM(y[k], x[k]);
            RE(y[k]) = RE(t) + RE(a[m]);
            IM(y[k]) = IM(t) + IM(a[m]);
        }
    }
    return RES_OK;
}

