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

#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "dspl.h"



#ifdef DOXYGEN_ENGLISH
/*! ****************************************************************************
\ingroup SPEC_MATH_POLY_GROUP
\fn int polyroots(double* a, int ord, complex_t* r, int* info)
\brief Function calculates real polynomial roots.

Function calculates roots of the real polynomial \f$P_N(x)\f$ order \f$N\f$ 
with `a` coefficient vector size `[(N+1) x 1]`. 
\f[
  P_N(x) = a_0 + a_1 x + a_2  x^2 + a_3  x^3 + ... a_N  x^N.
\f]

The roots of the polynomial are calculated as eigenvalues of the polynomial 
companion matrix. To calculate the eigenvalues, 
a subroutine of the LAPACK package is used.


\param[in]  a
Pointer to the vector of coefficients. \n
Vector size is `[ord+1 x 1]`. \n
Coefficient `a[0]` corresponds to the \f$a_0\f$ polynomial coefficient. \n 
Coefficient `a[ord]` cannot be zero. \n \n

\param[in]  ord
Polynomial order \f$N\f$.  \n \n

\param[out]  r
Pointer to the polynomial roots vector.  \n
Vector size is `[ord x 1]`. \n
Memory must be allocated. \n
The roots of a real polynomial can be either real or form simple 
or multiple complex conjugate pairs of roots. Therefore, the output 
root vector is of a complex data type. \n \n

\param[out]  info
Pointer to the LAPACK subroutine error code. \n 
This code is returned by the LAPACK subroutine and translated through 
this variable for analysis.. \n\n

\return
`RES_OK` --- roots are calculated successfully.  \n  
Else \ref ERROR_CODE_GROUP "code error".

Example:

\include polyroots_test.c

This program calculates the roots of the polynomial
\f[
P(x) = 2 + 2x + x^2
\f]
and prints the calculated roots.
The result of the program:

\verbatim
Error code: 0x00000000
r[0] = -1.00000 1.00000 j
r[1] = -1.00000-1.00000 j
\endverbatim

\author Sergey Bakhurin. www.dsplib.org
***************************************************************************** */
#endif
#ifdef DOXYGEN_RUSSIAN
/*! ****************************************************************************
\ingroup SPEC_MATH_POLY_GROUP
\fn int polyroots(double* a, int ord, complex_t* r, int* info)
\brief Расчет корней вещественного полинома

Функция рассчитывает корни полинома \f$P_N(x)\f$  \f$N-\f$ого 
порядка, заданного вектором коэффициентов `a`. 
\f[
  P_N(x) = a_0 + a_1 x + a_2  x^2 + a_3  x^3 + ... a_N  x^N.
\f]

Корни полинома рассчитываются как собственные числа характеристической 
матрицы полинома. Для расчета собственных чисел используется подпрограмма 
пакета LAPACK.


\param[in]  a
Указатель на вектор вещественных коэффициентов полинома. \n
Размер вектора `[ord+1 x 1]`. \n
Коэффициент `a[0]` соответствует коэффициенту полинома \f$a_0\f$. \n 
Коэффициент `a[ord]` не должен быть равен нулю. \n \n

\param[in]  ord
Порядок полинома \f$N\f$.  \n \n

\param[out]  r
Указатель на вектор комплексных корней полинома.  \n
Размер вектора `[ord x 1]`. \n
Память должна быть выделена. \n
Корни вещественного полинома могут быть как вещественными, 
так и образовывать простые или кратные комплексно-сопряженные пары корней.
Поэтому выходной вектор корней имеет комплексный тип данных.
\n \n

\param[out]  info
Указатель наа код возврата пакета LAPACK. \n 
Данный код возвращается подпрограммой LAPACK и транслируется через данную 
переменную для возможности анализа. \n\n

\return
`RES_OK` --- корни полинома рассчитаны успешно.  \n  
В противном случае \ref ERROR_CODE_GROUP "код ошибки".

Пример расчета корней полинома:

\include polyroots_test.c

Данная программа производит расчет корней полинома
\f[
P(x) = 2 + 2x + x^2
\f]
и выводит рассчитанные корни на печать.
Результат работы программы:

\verbatim
Error code: 0x00000000
r[0] = -1.00000 1.00000 j
r[1] = -1.00000-1.00000 j
\endverbatim

Получили пару комплексно-сопряженных корней полинома.

\author Бахурин Сергей. www.dsplib.org
***************************************************************************** */
#endif
int DSPL_API polyroots(double* a, int ord, complex_t* r, int* info)
{
    complex_t *t = NULL;
    int m;
    int err;
    
    if(!a || !r)
        return ERROR_PTR;
    if(ord<0)
        return ERROR_POLY_ORD;
    if(a[ord] == 0.0)
        return ERROR_POLY_AN;
    
    t = (complex_t*)malloc(ord * ord * sizeof(complex_t));
    if(!t)
        return ERROR_MALLOC;
    
    for(m = 0; m < ord-1; m++)
    {
        RE(t[m * (ord+1) + 1]) = 1.0;
        RE(t[m + ord * (ord - 1)]) = -a[m] / a[ord];
    }
    RE(t[ord * ord - 1]) = -a[ord-1] / a[ord];

    err = matrix_eig_cmplx(t, ord, r, info);
    
    if(t)
      free(t);
    return err;
}
