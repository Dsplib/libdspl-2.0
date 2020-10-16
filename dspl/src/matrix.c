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
#include "dspl.h"
#include "dspl_internal.h"
#include "blas.h"




#ifdef DOXYGEN_ENGLISH

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

#endif
#ifdef DOXYGEN_RUSSIAN
/*! ****************************************************************************
\ingroup SPEC_MATH_LINALG_GROUP
\fn int matrix_eye(double* a, int n, int m)
\brief Генерирование единичной вещественой матрицы размерности `n x m`.

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

#endif
#ifdef DOXYGEN_RUSSIAN

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

#endif
#ifdef DOXYGEN_RUSSIAN

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

