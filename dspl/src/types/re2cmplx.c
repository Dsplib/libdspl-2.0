/*
* Copyright (c) 2015-2020 Sergey Bakhurin
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
* MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
* GNU General Public License for more details.
*
* You should have received a copy of the GNU General Public License
* along with Foobar.  If not, see <http://www.gnu.org/licenses/>.
*/


#include <stdio.h>
#include <stdlib.h>
#include "dspl.h"



#ifdef DOXYGEN_ENGLISH
/*! *****************************************************************************
\ingroup TYPES_GROUP
\fn int re2cmplx(double* x, int n, complex_t *y)
\brief  Convert real array to the complex array.

Function copies the vector `x` to the real part of vector `y`.
Image part of the vector `y` sets as zero. \n
So complex vector contains data: \n
`y[i] = x[i] + j0, here i = 0,1,2 ... n-1`

\param[in]  x
Pointer to the real vector `x`. \n
Vector size is `[n x 1]`.  \n \n

\param[in]  n
Size of the real vector `x` and complex vector `y`. \n \n

\param[out] y
Pointer to the complex vector `y`. \n
Vector size is `[n x 1]`.  \n
Memory must be allocated.  \n \n

\return
`RES_OK` if function returns successfully.  \n
Else \ref ERROR_CODE_GROUP "code error": \n

Example:
\code{.cpp}
    double x[3] = {1.0, 2.0, 3.0};
    complex_t y[3];

    re2cmplx(x, 3, y);
\endcode

Vector `y` will keep:

\verbatim
    y[0] = 1+0j;
    y[1] = 2+0j;
    y[2] = 3+0j.
\endverbatim

\author Sergey Bakhurin. www.dsplib.org
***************************************************************************** */
#endif
#ifdef DOXYGEN_RUSSIAN
/*! ****************************************************************************
\ingroup TYPES_GROUP
\fn int re2cmplx(double* x, int n, complex_t *y)
\brief Преобразование массива вещественных данных в массив комплексных данных.

Функция заполняет реальные части массива `y` данных соответсвующими значениями 
исходного вещественного массива `x`.  \n  


\param[in]  x
Указатель на массив вещественных данных. \n
Размер массива `[n x 1]`.  \n \n

\param[in]  n
Размер массивов входных и выходных данных. \n \n

\param[out] y
Указатель на адрес массива комплексных данных. \n
Размер массива `[n x 1]`.  \n
Память должна быть выделена.  \n \n


\return
`RES_OK` если преобразование произведено успешно.  \n
В противном случае \ref ERROR_CODE_GROUP "код ошибки": \n

Например при выполнении следующего кода 
\code{.cpp}
    double x[3] = {1.0, 2.0, 3.0};
    complex_t y[3];

    re2cmplx(x, 3, y);
\endcode 

Значениям `y` будут присвоены значения:

\verbatim
    y[0] = 1+0j;
    y[1] = 2+0j;
    y[2] = 3+0j.
\endverbatim

\author Бахурин Сергей www.dsplib.org 
***************************************************************************** */
#endif
int DSPL_API re2cmplx(double* x, int n, complex_t* y)
{
    int k;
    if(!x || !y)
        return ERROR_PTR;
    if(n < 1)
        return ERROR_SIZE;

    for(k = 0; k < n; k++)
    {
        RE(y[k]) = x[k];
        IM(y[k]) = 0.0;
    }
    return RES_OK;
}


