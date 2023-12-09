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
\fn int matrix_mul(double* a, int na, int ma, double* b, int nb, int mb,
                   double* c)
\brief Matrix multiplication.

The function calculates the product of matrices \f$\mathbf{C} = \mathbf{AB}\f$,
here \f$\mathbf{A}\f$ -- matrix contains \f$N_A\f$ rows and
\f$M_A\f$ columns, \f$\mathbf{B}\f$ -- matrix contains \f$N_B\f$ rows and 
\f$M_B\f$ columns, product matrix \f$\mathbf{C} = \mathbf{AB}\f$ 
 contains\f$N_A\f$ rows and \f$M_B\f$ columns. 

\note
Matrix multiplication requires the equality \f$M_A = N_B\f$. \n
Function uses BLAS subroutine `dgemm`.


\param[in]  a
Pointer to the input matrix \f$\mathbf{A}\f$  size `na x ma`. \n
Matrix must be located in memory as column-major array. \n \n

\param[in]  na
Matrix `a` rows number. \n\n

\param[in]  ma
Matrix `a` columns number. \n\n

\param[in]  b
Pointer to the input matrix \f$\mathbf{B}\f$  size `nb x mb`. \n
Matrix must be located in memory as column-major array. \n \n

\param[in]  nb
Matrix `b` rows number. \n
Necessary equality `ma = nb`. \n\n

\param[in]  mb
Matrix `a` columns number. \n\n

\param[out] c
Pointer to the output matrix \f$\mathbf{С} = \mathbf{AB}\f$.\n
Matrix size is  `na x mb`. \n
Memory must be allocated. \n\n

\return
`RES_OK` --- function is calculated successfully. \n
Else \ref ERROR_CODE_GROUP "error code". \n

Example:
\include matrix_mul.c

The program forms and multiplies two matrices.
The original matrices and their product are printed.

The result of the program:
\verbatim

A = [ % size [4 x 5] type: real
0.00,     4.00,     8.00,    12.00,    16.00;
1.00,     5.00,     9.00,    13.00,    17.00;
2.00,     6.00,    10.00,    14.00,    18.00;
3.00,     7.00,    11.00,    15.00,    19.00;];

B = [ % size [5 x 3] type: real
0.00,     5.00,    10.00;
1.00,     6.00,    11.00;
2.00,     7.00,    12.00;
3.00,     8.00,    13.00;
4.00,     9.00,    14.00;];

C = [ % size [4 x 3] type: real
120.00,   320.00,   520.00;
130.00,   355.00,   580.00;
140.00,   390.00,   640.00;
150.00,   425.00,   700.00;]; 
\endverbatim

\author Sergey Bakhurin www.dsplib.org
***************************************************************************** */
#endif
#ifdef DOXYGEN_RUSSIAN
/*! ****************************************************************************
\ingroup SPEC_MATH_LINALG_GROUP
\fn int matrix_mul(double* a, int na, int ma, double* b, int nb, int mb,
                   double* c)
\brief Произведение вещественных матриц.

Функция рассчитывает произведение матриц \f$\mathbf{C} = \mathbf{AB}\f$,
где \f$\mathbf{A}\f$ --  матрица размерности \f$N_A\f$ строк и
\f$M_A\f$ столбцов, матрица размерности \f$N_B\f$ строк и
\f$M_B\f$ столбцов, а результирующая матрица \f$\mathbf{C} = \mathbf{AB}\f$ 
имеет размерность \f$N_A\f$ строк и
\f$M_B\f$ столбцов. 

\note
Для умножения матриц необходимо выполнение равенства \f$M_A = N_B\f$. \n
Функция использует подпрограмму `dgemm` пакета BLAS.


\param[in]  a
Указатель на матрицу \f$\mathbf{A}\f$  размерности `na x ma`. \n
Матрица должна быть расположена в памяти по столбцам. \n \n

\param[in]  na
Количество строк матрицы `a`. \n\n

\param[in]  ma
Количество столбцов матрицы `a`. \n\n

\param[in]  b
Указатель на матрицу \f$\mathbf{B}\f$  размерности `nb x mb`. \n
Матрица должна быть расположена в памяти по столбцам. \n \n

\param[in]  nb
Количество строк матрицы `b`. \n
Необходимо выполнение равенства `ma = nb`. \n\n

\param[in]  mb
Количество столбцов матрицы `b`. \n\n

\param[out] c
Указатель на матрицу \f$\mathbf{С} = \mathbf{AB}\f$.\n
Размер матрицы `na x mb`. \n
Память должна быть выделена. \n\n

\return
`RES_OK` --- функция выполнена успешно. \n
В противном случае \ref ERROR_CODE_GROUP "код ошибки". \n

Пример умножения матриц:
\include matrix_mul.c

Программа формирует и умножает две матрицы. 
Исходные матрицы и их произведение выводится на печать.

Результат работы программы:
\verbatim

A = [ % size [4 x 5] type: real
0.00,     4.00,     8.00,    12.00,    16.00;
1.00,     5.00,     9.00,    13.00,    17.00;
2.00,     6.00,    10.00,    14.00,    18.00;
3.00,     7.00,    11.00,    15.00,    19.00;];

B = [ % size [5 x 3] type: real
0.00,     5.00,    10.00;
1.00,     6.00,    11.00;
2.00,     7.00,    12.00;
3.00,     8.00,    13.00;
4.00,     9.00,    14.00;];

C = [ % size [4 x 3] type: real
120.00,   320.00,   520.00;
130.00,   355.00,   580.00;
140.00,   390.00,   640.00;
150.00,   425.00,   700.00;]; 
\endverbatim

\author Бахурин Сергей www.dsplib.org
***************************************************************************** */
#endif
int DSPL_API matrix_mul(double* a, int na, int ma, 
                        double* b, int nb, int mb,
                        double* c)
{
    
    double alpha = 1;
    double beta = 0.0;
    
    if(!a || !b || !c)
        return ERROR_PTR;
    if(na < 1 || ma < 1 || nb < 1 || mb < 1 || ma != nb)
        return ERROR_MATRIX_SIZE;

    /* BLAS DGEMM */
    dgemm_("N", "N", &na, &mb, &ma, &alpha, a, &na, b, &nb, &beta, c, &na);

    return RES_OK;
}

