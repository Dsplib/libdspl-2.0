/*
* Copyright (c) 2015-2019 Sergey Bakhurin
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
#include "dspl_internal.h"
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



#ifdef DOXYGEN_ENGLISH
/*! ****************************************************************************
\ingroup SPEC_MATH_LINALG_GROUP
\fn int matrix_eye(double* a, int n, int m)
\brief The real identity matrix size `n x m` generation.
 
 Function fills matrix `a` by zeros and sets main diagonal as ones.

\param[in]  a
Pointer to the real matrix size `n x m`. \n
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
\fn int matrix_eye(double* a, int n, int m)
\brief Генерирование единичной вещественной матрицы размерности `n x m`.

Данная функция заполняет матрицу нулями 
и записывает единицы на главной диагонали

\param[in]  a
Указатель на вещественную матрицу размерности `n x m`. \n
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
int DSPL_API matrix_eye(double* a, int n, int m)
{
    int    p, k;
    if(!a)
        return ERROR_PTR;
    if (n < 1 || m < 1)
        return ERROR_MATRIX_SIZE;
        
    k = 0;
    memset(a, 0, n*m*sizeof(double));
    for(p = 0; p < m; p++)
    {
        a[k] = 1.0;
        k += n+1;
    }

    return RES_OK;
}


#ifdef DOXYGEN_ENGLISH
/*! ****************************************************************************
\ingroup SPEC_MATH_LINALG_GROUP
\fn int matrix_eye_cmplx(complex_t* a, int n, int m)
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
\fn int matrix_eye_cmplx(complex_t* a, int n, int m)
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




#ifdef DOXYGEN_ENGLISH

#endif
#ifdef DOXYGEN_RUSSIAN

#endif
int DSPL_API matrix_pinv(double* a, int n, int m, double* tol, 
                         double* inv, int* info)
{
    int err, mn, i, j;
    double eps;


    double* u  = NULL;
    double* vt = NULL;
    double* v  = NULL;
    double* ut = NULL;
    double* s  = NULL;
    
    
    mn = (m > n) ? n : m;
    
    u = (double*) malloc(n*n*sizeof(double));
    if(!u)
    {
        err = ERROR_MALLOC;
        goto exit_label;
    }
    
    vt = (double*) malloc(m*m*sizeof(double));
    if(!vt)
    {
        err = ERROR_MALLOC;
        goto exit_label;
    }
    
    s = (double*) malloc(mn*sizeof(double));
    if(!s)
    {
        err = ERROR_MALLOC;
        goto exit_label;
    }
    
    err = matrix_svd(a, n, m, u, s, vt, info);
    if(err != RES_OK)
        goto exit_label;
    
    
      
    if(tol)
        eps = *tol;
    else
    {
        double smax;
        double mx = (n > m) ? (double)n : (double)m;
        err = minmax(s, mn, NULL, &smax);
        eps = DBL_EPSILON * mx * smax;
    }
    
    for(i = 0; i < mn; i++)
        if(s[i] > eps)
            s[i] = 1.0 / s[i];
        else
            s[i] = 0.0;

    v = (double*) malloc(m*m*sizeof(double));
    if(!v)
    {
        err = ERROR_MALLOC;
        goto exit_label;
    }
    err = matrix_transpose(vt, m, m, v);
    if(err != RES_OK)
        goto exit_label;
    
    ut = (double*) malloc(n*n*sizeof(double));
    if(!ut)
    {
        err = ERROR_MALLOC;
        goto exit_label;
    }
    err = matrix_transpose(u, n, n, ut);
    if(err != RES_OK)
        goto exit_label;
      
    for(i = 0; i < mn; i++)
        for(j = 0; j < m; j++)
            v[j + i*m] *= s[i];
    
    if(mn < m)
        memset(v+ mn*m, 0, (m-mn)*sizeof(double));
    
    err = matrix_mul(v, m, n, ut, n, n, inv);

exit_label:
    if(u)
        free(u);
    if(vt)
        free(vt);
    if(s)
        free(s);
    if(v)
        free(v);
    if(ut)
        free(ut);
    
    return err;
}








#ifdef DOXYGEN_ENGLISH

#endif
#ifdef DOXYGEN_RUSSIAN

#endif
int DSPL_API matrix_print(double* a, int n, int m, 
                          const char* name, const char* format)
{
    int p,q;
    
    if(!a)
        return ERROR_PTR;
    if(n < 1 || m < 1)
        return ERROR_SIZE;
        
    printf("\n%s = [ %% size [%d x %d] type: real", name, n, m);
    
    for(p = 0; p < n; p++)
    {
        printf("\n");
        for(q = 0; q < m; q++)
        {
            printf(format, a[q*n + p]);
            if(q == m-1)
                printf(";");
            else
                printf(", ");
        }
    }
    printf("];\n");
    
    return RES_OK;
}


#ifdef DOXYGEN_ENGLISH

#endif
#ifdef DOXYGEN_RUSSIAN

#endif
int DSPL_API matrix_print_cmplx(complex_t* a, int n, int m, 
                                const char* name, const char* format)
{
    int p,q;

    if(!a)
        return ERROR_PTR;
    if(n < 1 || m < 1)
        return ERROR_MATRIX_SIZE;

    if(!a)
        return ERROR_PTR;
    if(n < 1 || m < 1)
        return ERROR_SIZE;
        
    printf("\n%s = [ %% size [%d x %d] type: complex", name, n, m);
    
    for(p = 0; p < n; p++)
    {
        printf("\n");
        for(q = 0; q < m; q++)
        {
            printf(format, RE(a[q*n + p]), IM(a[q*n + p]));
            if(q == m-1)
                printf(";");
            else
                printf(", ");
        }
    }
    printf("];\n");

    return RES_OK;
}



#ifdef DOXYGEN_ENGLISH

#endif
#ifdef DOXYGEN_RUSSIAN

#endif
int DSPL_API matrix_svd(double* a, int n, int m, 
                        double* u, double* s, double* vt, int* info)
{
    char jobz = 'A';
    double* work = NULL;
    int* iwork = NULL;
    int lwork, mn, mx, err;
    int pi;
    
    if(!a || !u || !s || !vt)
        return ERROR_PTR;
    if(n < 1 || m < 1)
        return ERROR_SIZE;

    if(n > m)
    {
        mn = m;
        mx = n;
    }
    else
    { 
        mn = n; 
        mx = m;
    }

    err = RES_OK;
    
    lwork = 4 * mn * mn + 6 * mn + mx;
    work = (double*) malloc(lwork*sizeof(double));
    iwork = (int*) malloc(8*mn*sizeof(int));
    dgesdd_(&jobz, &n, &m, a, &n, s, u, &n, vt, &m, work, &lwork, iwork, &pi);
    
    if(info)
        *info = pi;
      
    if(pi)
        err = ERROR_LAPACK;
    if(work)
        free(work);
    if(iwork)
        free(iwork);
    return err;
}


#ifdef DOXYGEN_ENGLISH

#endif
#ifdef DOXYGEN_RUSSIAN

#endif
int DSPL_API matrix_transpose(double* a, int n, int m, double* b)
{
    int p, q, i, j, aind, bind;
    if(!a || !b)
        return ERROR_PTR;
    if(n < 1 || m < 1)
        return ERROR_MATRIX_SIZE;
        
        
    for(p = 0; p < n - DSPL_MATRIX_BLOCK; p+=DSPL_MATRIX_BLOCK)
    {
        for(q = 0; q < m - DSPL_MATRIX_BLOCK; q+=DSPL_MATRIX_BLOCK)
        {
            for(i = 0; i < DSPL_MATRIX_BLOCK; i++)
            {
                for(j = 0; j < DSPL_MATRIX_BLOCK; j++)
                {
                    aind = (q+j) * n + p + i;
                    bind = (p+i) * m + q + j;
                    b[bind] = a[aind];
                }
            }
        }
    }
    for(i = p; i < n; i++)
        for(j = 0; j < m; j++)
            b[i*m + j] = a[j*n+i];

    for(i = 0; i < p; i++)
        for(j = q; j < m; j++)
            b[i*m + j] = a[j*n+i];
            
    return RES_OK;
}





#ifdef DOXYGEN_ENGLISH

#endif
#ifdef DOXYGEN_RUSSIAN

#endif
int DSPL_API matrix_transpose_cmplx(complex_t* a, int n, int m, complex_t* b)
{
    int p, q, i, j, aind, bind;
    
    if(!a || !b)
        return ERROR_PTR;
    if(n < 1 || m < 1)
        return ERROR_MATRIX_SIZE;
        
    for(p = 0; p < n - DSPL_MATRIX_BLOCK; p+=DSPL_MATRIX_BLOCK)
    {
        for(q = 0; q < m - DSPL_MATRIX_BLOCK; q+=DSPL_MATRIX_BLOCK)
        {
            for(i = 0; i < DSPL_MATRIX_BLOCK; i++)
            {
                for(j = 0; j < DSPL_MATRIX_BLOCK; j++)
                {
                    aind = (q+j) * n + p + i;
                    bind = (p+i) * m + q + j;
                    RE(b[bind]) = RE(a[aind]);
                    IM(b[bind]) = IM(a[aind]);
                }
            }
        }
    }
    for(i = p; i < n; i++)
    {
        for(j = 0; j < m; j++)
        {
            RE(b[i*m + j]) = RE(a[j*n+i]);
            IM(b[i*m + j]) = IM(a[j*n+i]);
        }
    }

    for(i = 0; i < p; i++)
    {
        for(j = q; j < m; j++)
        {
            RE(b[i*m + j]) = RE(a[j*n+i]);
            IM(b[i*m + j]) = IM(a[j*n+i]);
        }
    }
    return RES_OK;
}





#ifdef DOXYGEN_ENGLISH

#endif
#ifdef DOXYGEN_RUSSIAN

#endif
int DSPL_API matrix_transpose_hermite(complex_t* a, int n, int m, complex_t* b)
{
    int p, q, i, j, aind, bind;
    
    if(!a || !b)
        return ERROR_PTR;
    if(n < 1 || m < 1)
        return ERROR_MATRIX_SIZE;
    
    for(p = 0; p < n - DSPL_MATRIX_BLOCK; p+=DSPL_MATRIX_BLOCK)
    {
        for(q = 0; q < m - DSPL_MATRIX_BLOCK; q+=DSPL_MATRIX_BLOCK)
        {
            for(i = 0; i < DSPL_MATRIX_BLOCK; i++)
            {
                for(j = 0; j < DSPL_MATRIX_BLOCK; j++)
                {
                    aind = (q+j) * n + p + i;
                    bind = (p+i) * m + q + j;
                    RE(b[bind]) =    RE(a[aind]);
                    IM(b[bind]) = -IM(a[aind]);
                }
            }
        }
    }
    for(i = p; i < n; i++)
    {
        for(j = 0; j < m; j++)
        {
            RE(b[i*m + j]) = RE(a[j*n+i]);
            IM(b[i*m + j]) = -IM(a[j*n+i]);
        }
    }

    for(i = 0; i < p; i++)
    {
        for(j = q; j < m; j++)
        {
            RE(b[i*m + j]) =    RE(a[j*n+i]);
            IM(b[i*m + j]) = -IM(a[j*n+i]);
        }
    }
    
    return RES_OK;
}




#ifdef DOXYGEN_ENGLISH

#endif
#ifdef DOXYGEN_RUSSIAN

#endif
int DSPL_API vector_dot(double* x, double* y, int n, double* p)
{
    int inc = 1;
    
    if(!x || !y || !p)
        return ERROR_PTR;
    if(n<1)
        return ERROR_SIZE;

    *p = ddot_(&n, x, &inc, y, &inc); 

    return RES_OK;
}

