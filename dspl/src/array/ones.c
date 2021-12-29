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
\fn int ones(double* x, int n)
\brief Function fills all real vector `x` by ones values.

\param[in, out] x
Pointer to the vector `x`. \n
Vector size is `[n x 1]`. \n
All elements on this vector will be set to one. \n
\n

\param[in] n
Size of vector `x`. \n
\n

\return
`RES_OK` if function returns successfully. \n
 Else \ref ERROR_CODE_GROUP "error code".

Example:
\code{.cpp}
double y[5] = {0};
int i;
ones(y, 5);
for(i = 0; i < 5; i++)
    printf("%6.1f%    ", y[i]);
\endcode
 \n
Vector `y` values are:
\verbatim
    1.0    1.0    1.0    1.0    1.0
\endverbatim

\author Sergey Bakhurin www.dsplib.org
***************************************************************************** */
#endif
#ifdef DOXYGEN_RUSSIAN
/*! ****************************************************************************
\ingroup ARRAY_GROUP
\fn int ones(double* x, int n)
\brief Функция заполнения вещественного массива единицами

\param[in, out] x
Указатель на вещественный вектор `x`. \n
Размер вектора `[n x 1]`. \n
Значения данного вектора будут установлены в единицу. \n
\n

\param[in] n
Размер вектора `x`. \n
\n

\return
`RES_OK` если функция выполнена успешно. \n
В противном случае \ref ERROR_CODE_GROUP "код ошибки".

Пример:
\code{.cpp}
double y[5] = {0};
int i;
ones(y, 5);
for(i = 0; i < 5; i++)
    printf("%6.1f%    ", y[i]);
\endcode
 \n
Результат выполнения:
\verbatim
    1.0    1.0    1.0    1.0    1.0
\endverbatim

\author Бахурин Сергей www.dsplib.org
***************************************************************************** */
#endif
int DSPL_API ones(double* x, int n)
{
    int i;
    if(!x)
        return ERROR_PTR;
    if(n<1)
        return ERROR_SIZE;
    for(i = 0; i < n; i++)
        x[i] = 1.0;
    return RES_OK;
}
