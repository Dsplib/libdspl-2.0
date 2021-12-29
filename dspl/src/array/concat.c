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
\fn int concat(void* a, size_t na, void* b, size_t nb, void* c)
\brief
Concatenate arrays `a` and `b`

Let's arrays `a` and `b` are vectors: \n
`a = [a(0), a(1), ... a(na-1)]`, \n
`b = [b(0), b(1), ... b(nb-1)]`, \n
concatenation of these arrays will be array `c` size `na+nb`: \n
`c = [a(0), a(1), ... a(na-1), b(0), b(1), ... b(nb-1)]`.


\param[in] a
Pointer to the first array `a`. \n
Array `a` size is `na` bytes. \n
\n

\param[in] na
Array `a` size (bytes). \n
\n

\param[in] b
Pointer to the second array `b`. \n
Array `b` size is `nb` bytes. \n
\n

\param[in] nb
Array `a` size (bytes). \n
\n

\param[out] c
Pointer to the concatenation result array `c`. \n
Array `c` size is `na + nb` bytes. \n
Memory must be allocated. \n
\n

\return
`RES_OK` if function returns successfully. \n
Else \ref ERROR_CODE_GROUP "code error".

Function uses pointer type `void*` and can be useful for an arrays
concatenation with different types. \n
For example two `double` arrays concatenation:
\code{.cpp}
double a[3] = {1.0, 2.0, 3.0};
double b[2] = {4.0, 5.0};
double c[5];

concat((void*)a, 3*sizeof(double), (void*)b, 2*sizeof(double), (void*)c);
\endcode
Vector `c` keeps follow data:
\verbatim
c = [1.0, 2.0, 3.0, 4.0, 5.0]
\endverbatim

\author Sergey Bakhurin www.dsplib.org
***************************************************************************** */
#endif

#ifdef DOXYGEN_RUSSIAN
/*! ****************************************************************************
\ingroup ARRAY_GROUP
\fn int concat(void* a, size_t na, void* b, size_t nb, void* c)
\brief Конкатенация двух массивов данных

Функция производит конкатенацию двух массивов. Пусть массивы `a` и `b`
заданы как векторы: \n
`a = [a(0), a(1), ... a(na-1)]`,    \n
`b = [b(0), b(1), ... b(nb-1)]`,    \n
тогда результатом конкатенации будет вектор размера `na+nb` вида: \n
`c = [a(0), a(1), ... a(na-1), b(0), b(1), ... b(nb-1)]`.


\param[in] a
Указатель на первый вектор `a`. \n
Размер вектора `na` байт. \n \n

\param[in] na
Размер первого вектора `a` в байт. \n \n

\param[in] b
Указатель на второй вектор `b`. \n
Размер памяти вектора `nb` байт. \n \n

\param[in] nb
Размер второго вектора `b` в байт. \n \n

\param[out] c
Указатель на вектор конкатенации `c`. \n
Размер памяти вектора `na + nb` байт. \n
Память должна быть выделена. \n \n

\return
`RES_OK` если функция выполнена успешно. \n
 В противном случае \ref ERROR_CODE_GROUP "код ошибки".

\note
Функция использует указатели типа `void*` и может быть использована для
конкатенации данных различного типа. \n
Например конкатенация массивов типа `double`:
\code{.cpp}
double a[3] = {1.0, 2.0, 3.0};
double b[2] = {4.0, 5.0};
double c[5];
concat((void*)a, 3*sizeof(double), (void*)b, 2*sizeof(double), (void*)c);
\endcode
в результате вектор `c` будет хранить массив данных:
\verbatim
c = [1.0, 2.0, 3.0, 4.0, 5.0]
\endverbatim

\author
Бахурин Сергей
www.dsplib.org
***************************************************************************** */
#endif
int DSPL_API concat(void* a, size_t na, void* b, size_t nb, void* c)
{
    if(!a || !b || !c || c == b)
        return ERROR_PTR;
    if(na < 1 || nb < 1)
        return ERROR_SIZE;

    if(c != a)
        memcpy(c, a, na);

    memcpy((char*)c+na, b, nb);
    return RES_OK;
}


