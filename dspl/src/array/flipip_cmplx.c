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
\fn int flipip_cmplx(complex_t* x, int n)
\brief Flip complex vector `x` in place

Function flips complex vector `x` length `n` in the memory
 \n
For example complex vector `x`    length 6: \n
\verbatim
x = [0+0j, 1+1j, 2+2j, 3+3j, 4+4j, 5+5j]
\endverbatim
After flipping it will be as follow:
\verbatim
x = [5+5j, 4+4j, 3+3j, 2+2j, 1+1j, 0+0j]
\endverbatim

\param[in, out] x
Pointer to the complex vector `x`. \n
Vector size is `[n x 1]`. \n
Flipped vector will be on the same address. \n

\param[in] n
Length of the vector `x`. \n \n

\return
`RES_OK` if function returns successfully. \n
Else \ref ERROR_CODE_GROUP "error code".

Example:
\code{.cpp}
complex_t y[5] = {{0.0, 0.0}, {1.0, 1.0}, {2.0, 2.0}, {3.0, 3.0}, {4.0, 4.0}};
for(i = 0; i < 5; i++)
    printf("%6.1f%+.1fj    ", RE(y[i]), IM(y[i]));
flipip_cmplx(y, 5);
printf("\n");
for(i = 0; i < 5; i++)
    printf("%6.1f%+.1fj    ", RE(y[i]), IM(y[i]));
\endcode
 \n
Program result:
\verbatim
0.0+0.0j         1.0+1.0j         2.0+2.0j         3.0+3.0j         4.0+4.0j
4.0+4.0j         3.0+3.0j         2.0+2.0j         1.0+1.0j         0.0+0.0j
\endverbatim

\author Sergey Bakhurin www.dsplib.org
***************************************************************************** */
#endif
#ifdef DOXYGEN_RUSSIAN
/*! ****************************************************************************
\ingroup ARRAY_GROUP
\fn int flipip_cmplx(complex_t* x, int n)
\brief Функция отражения комплексного вектора `x`

Функция производит отражение комплексного вектора длины `n`
в памяти данных. \n
Например исходный вектор `x`    длины 6: \n
\verbatim
x = [0+0j, 1+1j, 2+2j, 3+3j, 4+4j, 5+5j]
\endverbatim
После отражения вектор `x` будет иметь вид:
\verbatim
x = [5+5j, 4+4j, 3+3j, 2+2j, 1+1j, 0+0j]
\endverbatim

\param[in, out] x
Указатель на комплексный вектор `x`. \n
Размер вектора `[n x 1]`. \n
Результат отражения будет помещен по этому же адресу. \n
\n

\param[in] n
Размер вектора `x`. \n
\n

\return
`RES_OK` если функция выполнена успешно. \n
 В противном случае \ref ERROR_CODE_GROUP "код ошибки".

Пример:
\code{.cpp}
complex_t y[5] = {{0.0, 0.0}, {1.0, 1.0}, {2.0, 2.0}, {3.0, 3.0}, {4.0, 4.0}};
for(i = 0; i < 5; i++)
    printf("%6.1f%+.1fj    ", RE(y[i]), IM(y[i]));
flipip_cmplx(y, 5);
printf("\n");
for(i = 0; i < 5; i++)
    printf("%6.1f%+.1fj    ", RE(y[i]), IM(y[i]));
\endcode
 \n
Результат выполнения:
\verbatim
0.0+0.0j         1.0+1.0j         2.0+2.0j         3.0+3.0j         4.0+4.0j
4.0+4.0j         3.0+3.0j         2.0+2.0j         1.0+1.0j         0.0+0.0j
\endverbatim

\author
Бахурин Сергей
www.dsplib.org
***************************************************************************** */
#endif
int DSPL_API flipip_cmplx(complex_t* x, int n)
{
    int k;
    complex_t tmp;
    if(!x)
        return ERROR_PTR;
    if(n<1)
        return ERROR_SIZE;

    for(k = 0; k < n/2; k++)
    {
        RE(tmp) = RE(x[k]);
        RE(x[k]) = RE(x[n-1-k]);
        RE(x[n-1-k]) = RE(tmp);

        IM(tmp) = IM(x[k]);
        IM(x[k]) = IM(x[n-1-k]);
        IM(x[n-1-k]) = IM(tmp);
    }
    return RES_OK;
}



