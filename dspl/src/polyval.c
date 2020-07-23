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
#include "dspl.h"



#ifdef DOXYGEN_ENGLISH

#endif
#ifdef DOXYGEN_RUSSIAN

#endif
int DSPL_API poly_z2a_cmplx(complex_t* z, int nz, int ord, complex_t* a)
{
    int k, ind, res;
    complex_t x[2];

    if(!z || !a)
        return ERROR_PTR;
    if(nz < 0)
        return ERROR_SIZE;
    if(nz > ord || ord < 1)
        return ERROR_POLY_ORD;

    RE(x[1]) = 1.0;
    IM(x[1]) = 0.0;

    memset(a, 0, (ord+1) * sizeof(complex_t));

    RE(a[0]) = 1.0;
    ind = 1;
    for(k = 0; k < nz; k++)
    {
        RE(x[0]) = -RE(z[k]);
        IM(x[0]) = -IM(z[k]);
        res = conv_cmplx(a, ind, x, 2, a);
        if(res!=RES_OK)
            return res;
        ind++;
    }

    return RES_OK;
}





#ifdef DOXYGEN_ENGLISH

#endif
#ifdef DOXYGEN_RUSSIAN
/*! ****************************************************************************
\ingroup SPEC_MATH_COMMON_GROUP
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
Память должна быть выделена. \n \n

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




#ifdef DOXYGEN_ENGLISH

#endif
#ifdef DOXYGEN_RUSSIAN
/*! ****************************************************************************
\ingroup SPEC_MATH_COMMON_GROUP
\fn int polyval(double* a, int ord, double* x, int n, double* y)
\brief Расчет вещественного полинома

Функция рассчитывает полином \f$P_N(x)\f$  \f$N-\f$ого порядка для вещественного
аргумента, заданного вектором `x`. 
\f[
  P_N(x) = a_0 + a_1 \cdot x + a_2 \cdot x^2 + 
  a_3 \cdot x^3 + ... a_N \cdot x^N.
\f]

Для расчета используется формула Горнера:
\f[
  P_N(x) = a_0 + x \cdot (a_1 + x \cdot (a_2 + \cdot 
  ( \ldots x \cdot (a_{N-1} + x\cdot a_N) \ldots ))) 
\f]

\param[in]  a
Указатель на вектор вещественных коэффициентов полинома. \n
Размер вектора `[ord+1 x 1]`. \n
Коэффициент `a[0]` соответствует коэффициенту полинома \f$a_0\f$. \n \n

\param[in]  ord
Порядок полинома \f$N\f$.  \n \n

\param[in]  x
Указатель на вектор аргумента полинома.  \n
Размер вектора `[n x 1]`. \n
Значения полинома будут расчитаны для всех значений аргумента вектора `x`. \n\n

\param[in]  n
Размер вектора агрумента полинома.  \n \n

\param[out]  y
Указатель на значения полинома для аргумента `x`. \n
Размер вектора `[n x 1]`. \n
Память должна быть выделена. \n\n

\return
`RES_OK` --- полином рассчитан успешно.  \n  
В противном случае \ref ERROR_CODE_GROUP "код ошибки".

\author Бахурин Сергей. www.dsplib.org
***************************************************************************** */
#endif
int DSPL_API polyval(double* a, int ord, double* x, int n, double* y)
{
    int k, m;

    if(!a || !x || !y)
        return ERROR_PTR;
    if(ord<0)
        return ERROR_POLY_ORD;
    if(n<1)
        return ERROR_SIZE;

    for(k = 0; k < n; k++)
    {
        y[k] = a[ord];
        for(m = ord-1; m>-1; m--)
            y[k] = y[k]*x[k] + a[m];
    }
    return RES_OK;
}





#ifdef DOXYGEN_ENGLISH

#endif
#ifdef DOXYGEN_RUSSIAN
/*! ****************************************************************************
\ingroup SPEC_MATH_COMMON_GROUP
\fn int polyval_cmplx(complex_t* a, int ord, complex_t* x, int n, complex_t* y)
\brief Расчет комплексного полинома

Функция рассчитывает полином \f$P_N(x)\f$ \f$N\f$-го порядка 
комплексного аргумента, заданного вектором `x`. \n

\f[
P_N(x) = a_0 + a_1 \cdot x + a_2 \cdot x^2 + a_3 \cdot x^3 + ... a_N \cdot x^N.
\f]

Для расчета используется формула Горнера: \n
\f[
  P_N(x) = a_0 + x \cdot (a_1 + x \cdot (a_2 + \cdot 
  ( \ldots x \cdot (a_{N-1} + x\cdot a_N) \ldots ))) 
\f]

\param[in]  a
Указатель на вектор комплексных коэффициентов полинома. \n
Размер вектора `[ord+1 x 1]`. \n
Коэффициент `a[0]` соответствует коэффициенту полинома \f$a_0\f$. \n \n

\param[in]  ord
Порядок полинома \f$N\f$.  \n \n

\param[in]  x
Указатель на вектор аргумента полинома. \n
Размер вектора `[n x 1]`. \n
Значения полинома будут расчитаны для всех значений аргумента вектора `x`. \n \n

\param[in]  n
Размер вектора агрумента полинома. \n\n

\param[out]  y
Указатель вектор значения полинома для аргумента `x`. \n
Размер вектора `[n x 1]`. \n
Память должна быть выделена. \n \n

\return
`RES_OK` --- полином расчитан успешно. \n
В противном случае \ref ERROR_CODE_GROUP "код ошибки".

\author Бахурин Сергей. www.dsplib.org
***************************************************************************** */
#endif
int DSPL_API polyval_cmplx(complex_t* a, int ord,
                           complex_t* x, int n, complex_t* y)
{
    int k, m;
    complex_t t;

    if(!a || !x || !y)
        return ERROR_PTR;
    if(ord<0)
        return ERROR_POLY_ORD;
    if(n<1)
        return ERROR_SIZE;

    for(k = 0; k < n; k++)
    {
        RE(y[k]) = RE(a[ord]);
        IM(y[k]) = IM(a[ord]);
        for(m = ord-1; m>-1; m--)
        {
            RE(t) = CMRE(y[k], x[k]);
            IM(t) = CMIM(y[k], x[k]);
            RE(y[k]) = RE(t) + RE(a[m]);
            IM(y[k]) = IM(t) + IM(a[m]);
        }
    }
    return RES_OK;
}

