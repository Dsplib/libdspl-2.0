/*
* Copyright (c) 2015-2022 Sergey Bakhurin
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
\fn int matrix_eig_cmplx(complex_t* a, int n, complex_t* v, int* info)

\brief Eigenvalues calculation of the complex square matrix `a`.

Function calculates `n` eigenvalues of the matrix `a` size `n x n`. 

\param[in]  a
Pointer to the complex matrix `a` size `n x n`. \n
Matrix is stored in the memory as column-major array. \n\n

\param[in]  n
Size of matrix `n x n`.\n\n

\param[out]  v
Pointer to the eigenvalues vector. \n 
Vector size is `]n x 1]`. \n 
Memory must be allocated. \n\n

\param[out]  info
Pointer to the `zgees`  LAPACK subroutine output parameter. \n
If an error occurs while calculating the vector of eigenvalues, 
the LAPACK subroutine returns an error code that can be read 
from this pointer. \n\n

\return
`RES_OK` --- function is calculated successfully. \n
Else \ref ERROR_CODE_GROUP "code error". \n
If an `ERROR_LAPACK` error occurs, the LAPACK package 
error code will be written to the `info` address. \n


Eigenvalues calculation example:
\include matrix_eig.c

This program calculates eigenvalues of the matrix size  `3 x 3` 
and print its to display.\n

Result:
\verbatim
A = [ % size [3 x 3] type: complex
1.00   +0.00i,     2.00   +0.00i,     3.00   +0.00i;
1.00   +0.00i,     0.00   +0.00i,     0.00   +0.00i;
0.00   +0.00i,     1.00   +0.00i,     0.00   +0.00i;];

v = [ % size [3 x 1] type: complex
 2.374424 -0.000000i;    
-0.687212 +0.889497i;
-0.687212 -0.889497i;];
\endverbatim

\author Sergey Bakhurin www.dsplib.org
***************************************************************************** */
#endif
#ifdef DOXYGEN_RUSSIAN
/*! ****************************************************************************
\ingroup SPEC_MATH_LINALG_GROUP
\fn int matrix_eig_cmplx(complex_t* a, int n, complex_t* v, int* info)

\brief Расчет собственных значений квадратной комплексной матрицы.

Данная функция производит расчет `n` собственных значений квадратной матрицы 
размером `n x n`. 

\param[in]  a
Указатель на комплексную матрицу размерности `n x n`. \n
Матрица должна быть расположена в памяти по столбцам. \n\n

\param[in]  n
Размерность квадратной матрицы.\n

\param[out]  v
Указатель на вектор собственных значений матрицы. \n 
Размер вектора `n x 1`. \n 
Память должна быть выделена. \n\n

\param[out]  info
Указатель на код возврата функции `zgees` пакета LAPACK. \n
В случае возникновения ошибки при расчете вектора собственных значений,
пакет LAPACK возвращает код ошибки, который может быть прочитан по данному
указателю. \n

\return
`RES_OK` --- функция выполнена успешно. \n
В противном случае \ref ERROR_CODE_GROUP "код ошибки". \n
При возникновении ошибки `ERROR_LAPACK` по адресу  
`info` будет записан код ошибки пакета LAPACK. \n


Пример расчета собственных значений матрицы:
\include matrix_eig.c

Данная программа рассчитывает собственные значения матрицы размерности `3 x 3`
и выводит собственные значения на печать. \n

Результат работы программы:
\verbatim
A = [ % size [3 x 3] type: complex
1.00   +0.00i,     2.00   +0.00i,     3.00   +0.00i;
1.00   +0.00i,     0.00   +0.00i,     0.00   +0.00i;
0.00   +0.00i,     1.00   +0.00i,     0.00   +0.00i;];

v = [ % size [3 x 1] type: complex
 2.374424 -0.000000i;    
-0.687212 +0.889497i;
-0.687212 -0.889497i;];
\endverbatim

\author Бахурин Сергей www.dsplib.org
***************************************************************************** */
#endif
int DSPL_API matrix_eig_cmplx(complex_t* a, int n, complex_t* v, int* info)
{
    int err;
    int sdim = 0;
    int ldvs = 1;
    int lwork = 2*n;
    if(!a || !v)
        return ERROR_PTR;
    
    if(n<1)
        return ERROR_MATRIX_SIZE;
    
    complex_t *work=(complex_t*)malloc(lwork*sizeof(complex_t)); 
    double *rwork = (double*)malloc(n*sizeof(double)); 

    zgees_("N", "N", NULL, &n, a, &n, &sdim, v, NULL, &ldvs, work, &lwork, 
           rwork, NULL, &err);
    
    if(err!=0)
    {
        if(info)
            *info = err;
        err = ERROR_LAPACK;
    }
    else
        err = RES_OK;
    
    free(work);
    free(rwork);
    return err;
}

