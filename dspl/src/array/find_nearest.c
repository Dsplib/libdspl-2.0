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
\brief Find nearest to `val` vector `x` element.

\param[in] x
Pointer to the vector `x`. \n
Vector size is `[n x 1]`. \n    
\n

\param[in] n
Size of vector `x`. \n
\n

\param[in] val \n
Value nearest to which vector `x` element need to find. 

\param[out] idx
Pointer to the variable type `int` which will return index of the  
nearest to `val` vector `x` element. \n
Pointer can be `NULL`. 
Index doesn't return in this case.
\n

\param[out] dist
Pointer to the variable type `double` which will return
distance  from the nearest vector `x` element to the `val` . \n
Pointer can be `NULL`. 
Distance doesn't return in this case.
\n

\return
`RES_OK` if function calculated successfully. \n
Else \ref ERROR_CODE_GROUP "code error".

\author Sergey Bakhurin www.dsplib.org
***************************************************************************** */
#endif
#ifdef DOXYGEN_RUSSIAN
/*! ****************************************************************************
\ingroup ARRAY_GROUP
\brief Найти ближайший `val` элемент массива `x`.

\param[in] x
Указатель на вещественный вектор `x`. \n
Размер вектора `[n x 1]`. \n    
\n

\param[in] n
Размер вектора `x`. \n
\n

\param[in] val \n
Значение к которому ищется ближайший элемент массива. 

\param[out] idx
Указатель на переменную типа `int` в которую будет записан 
индекс ближайшего к `val` элемента массива `x`. \n
Указатель может быть `NULL`. 
Индекс в этом случае не возвращается.
\n

\param[out] dist
Указатель на переменную типа `double` в которую будет записано
расстояние от ближайшего элемента массива `x` до `val` . \n
Указатель может быть `NULL`. 
Значение расстояние от ближайшего элемента массива `x` до `val` 
в этом случае не возвращается.
\n

\return
`RES_OK` если функция выполнена успешно. \n
В противном случае \ref ERROR_CODE_GROUP "код ошибки".

\author Бахурин Сергей www.dsplib.org
***************************************************************************** */
#endif
int DSPL_API find_nearest(double* x, int n, double val, int* idx, double* dist)
{
    double mind, dv;
    int iv, i;
    
    if(!x)
        return ERROR_PTR;
    if(n < 1)
        return ERROR_SIZE;
    
    mind = fabs(x[0] - val);
    iv = 0;
    for(i = 1; i < n; i++)
    {
        dv = fabs(x[i] - val);
        if( dv < mind)
        {
            mind =  dv;
            iv = i;
        }
    }
    
    if(idx)
        *idx = iv;
    if(dist)
        *dist = fabs(x[iv] - val);
    
    return RES_OK;
    
}