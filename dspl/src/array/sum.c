/*
* Copyright (c) 2015-2024 Sergey Bakhurin
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

\author Sergey Bakhurin www.dsplib.org
***************************************************************************** */
#endif
#ifdef DOXYGEN_RUSSIAN
/*! ****************************************************************************
\ingroup ARRAY_GROUP
\brief Расчет суммы элеметов массива

\param[in] x
Указатель на вещественный вектор `x`. \n
Размер вектора `[n x 1]`. \n    
\n

\param[in] n
Размер вектора `x`. \n
\n

\param[in] n
Указатель на переменную суммы элементов массива. \n
По данному указателю будет записан результат работы функции.
\n


\return
`RES_OK` если функция выполнена успешно. \n
В противном случае \ref ERROR_CODE_GROUP "код ошибки".

Пример:
\code{.cpp}
double y[5] = {0.0, 1.0, 2.0, 3.0, 4.0};
double s;
sum(y, 5, &s);
printf("s = %6.1f% \n", s);
\endcode
 \n
Результат выполнения:
\verbatim
    s = 10
\endverbatim

\author Бахурин Сергей www.dsplib.org
***************************************************************************** */
#endif
int DSPL_API sum(double* x, int n, double* s)
{
    int i;
    double z = 0.0;
    if(!x || !s)
        return ERROR_PTR;
    if(n<1)
        return ERROR_SIZE;
      
    for(i = 0; i < n; i++)
        z += x[i];
    *s = z;
    return RES_OK;
}