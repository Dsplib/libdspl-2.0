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
#include "dspl.h"

#ifdef DOXYGEN_ENGLISH
/*! ****************************************************************************
\ingroup SPEC_MATH_COMMON_GROUP

\brief Function \f$ \textrm{sinc}(x,a) = \frac{\sin(ax)}{ax}\f$ 
for the real vector `x`.

\param[in] x
Pointer to the input vector \f$ x \f$. \n
Vector size is `[n x 1]`. \n \n

\param[in] n
Input and output vectors size. \n \n

\param[in]    a
Function parameter \f$ \textrm{sinc}(x,a) = \frac{\sin(ax)}{ax}\f$. \n\n

\param[out]    y
Pointer to the `sinc` function output vector. \n
Vector size is `[n x 1]`. \n
Memory must be allocated. \n \n

\return
`RES_OK` if function calculated successfully.    \n
Else \ref ERROR_CODE_GROUP "code error". \n

\author Sergey Bakhurin www.dsplib.org
***************************************************************************** */
#endif
#ifdef DOXYGEN_RUSSIAN
/*! ****************************************************************************
\ingroup SPEC_MATH_COMMON_GROUP

\brief Функция \f$ \textrm{sinc}(x,a) = \frac{\sin(ax)}{ax}\f$.

Функция рассчитывает значения функции для вещественного вектора `x`.

\param[in]  x
Указатель на вектор переменной \f$ x \f$. \n
Размер вектора `[n x 1]`. \n
Память должна быть выделена. \n \n

\param[in]  n
Размер входного вектора `x`. \n \n

\param[in]  a
Параметр функции \f$ \textrm{sinc}(x,a) = \frac{\sin(ax)}{ax}\f$. \n\n


\param[out]  y
Указатель на вектор значений функции. \n
 Размер вектора `[n x 1]`. \n
 Память должна быть выделена. \n  \n


\return
`RES_OK` --- расчёт произведен успешно. \n
В противном случае \ref ERROR_CODE_GROUP "код ошибки". \n

\author Бахурин Сергей www.dsplib.org
***************************************************************************** */
#endif
int DSPL_API sinc(double* x, int n, double a, double* y)
{
    int k;

    if(!x || !y)
        return ERROR_PTR;
    if(n<1)
        return ERROR_SIZE;

    for(k = 0; k < n; k++)
        y[k] = (x[k]==0.0) ? 1.0 : sin(a*x[k])/(a*x[k]);

    return RES_OK;


}


