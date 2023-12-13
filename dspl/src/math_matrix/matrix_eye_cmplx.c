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
#include <float.h>
#include "dspl.h"

#include "blas.h"



#ifdef DOXYGEN_ENGLISH
/*! ****************************************************************************
\ingroup SPEC_MATH_LINALG_GROUP

\brief The complex identity matrix size `n x m` generation.
 
 Function fills matrix `a` by zeros and sets main diagonal as ones.

\param[in]  a
Pointer to the complex matrix size `n x m`. \n
Matrix is stored in the memory as column-major array. \n \n

\param[in]  n
Matrix `a` rows number. \n\n

\param[in]  m
Matrix `a` columns number.. \n\n

\return
`RES_OK` --- function is calculated successfully. \n
Else \ref ERROR_CODE_GROUP "code error". \n

\author Sergey Bakhurin www.dsplib.org
***************************************************************************** */
#endif
#ifdef DOXYGEN_RUSSIAN
/*! ****************************************************************************
\ingroup SPEC_MATH_LINALG_GROUP

\brief Генерирование единичной комплексной матрицы размерности `n x m`.

Данная функция заполняет матрицу нулями 
и записывает единицы на главной диагонали

\param[in]  a
Указатель на комплексную матрицу размерности `n x m`. \n
Матрица должна быть расположена в памяти по столбцам. \n \n

\param[in]  n
Количество строк матрицы. \n\n

\param[in]  m
Количество столбцов матрицы. \n\n

\return
`RES_OK` --- функция выполнена успешно. \n
В противном случае \ref ERROR_CODE_GROUP "код ошибки". \n

\author Бахурин Сергей www.dsplib.org
***************************************************************************** */
#endif
int DSPL_API matrix_eye_cmplx(complex_t* a, int n, int m)
{
    int p, k;
    if(!a)
        return ERROR_PTR;
    if (n < 1 || m < 1)
        return ERROR_MATRIX_SIZE;
        
    k = 0;
    memset(a, 0, n*m*sizeof(complex_t));
    for(p = 0; p < m; p++)
    {
        RE(a[k]) = 1.0;
        k += n+1;
    }

    return RES_OK;
}


