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
/*! ****************************************************************************
\ingroup FILTER_CONV_GROUP
\fn int filter_iir(double* b, double* a, int ord, double* x, int n, double* y)
\brief Real IIR filtration

Function calculates real IIR filter output for real signal. The real filter
contains real coefficients of the transfer function \f$H(z)\f$
 numerator and denominator:
\f[
  H(z) = \frac{\sum_{n = 0}^{N} b_n  z^{-n}}
  {1+{\frac{1}{a_0}}\sum_{m = 1}^{M} a_m  z^{-n}},
\f]
here \f$a_0\f$ cannot be equals zeros, \f$N=M=\f$`ord`.

\param[in]  b
Pointer to the vector \f$b\f$ of IIR filter 
transfer function numerator coefficients. \n 
Vector size is `[ord + 1 x 1]`. \n \n 

\param[in]  a
Pointer to the vector \f$a\f$ of IIR filter 
transfer function denominator coefficients. \n 
Vector size is `[ord + 1 x 1]`. \n 
This pointer can be `NULL` if filter is FIR. \n \n 

\param[in]  ord
Filter order. Number of the transfer function 
numerator and denominator coefficients 
(length of vectors `b` and `a`) is `ord + 1`. \n \n 

\param[in]  x
Pointer to the input signal vector. \n 
Vector size is `[n x 1]`. \n \n 

\param[in]  n
Size of the input signal vector `x`. \n \n 

\param[out] y
Pointer to the IIR filter output vector. \n 
Vector size is `[n x  1]`. \n 
Memory must be allocated. \n \n 

\return
`RES_OK` if filter output is calculated successfully. \n
Else \ref ERROR_CODE_GROUP "code error". \n

Example:

\include filter_iir_test.c

Input signal is
\f$s(t) = \sin(2\pi \cdot 0.05 t) + n(t)\f$, here \f$n(t)\f$ white Gaussian
noise with zero mean value and unit standard deviation. \n

Input signal is filtered by elliptic LPF order 6 and output signal and data 
saves in the txt-files  

\verbatim
dat/s.txt  - input signal + noise
dat/sf.txt - filter output.
\endverbatim

Plots:

\image html filter_iir_test.png

GNUPLOT script for make plots is:
\include filter_iir.plt

\author Sergey Bakhurin www.dsplib.org
***************************************************************************** */
#endif
#ifdef DOXYGEN_RUSSIAN
/*! ****************************************************************************
\ingroup FILTER_CONV_GROUP
\fn int filter_iir(double* b, double* a, int ord, double* x, int n, double* y)
\brief Фильтрация вещественного сигнала вещественным БИХ-фильтром

Функция рассчитывает выход фильтра заданного выражением
\f[
  H(z) = \frac{\sum_{n = 0}^{N} b_n  z^{-n}}
  {1+{\frac{1}{a_0}}\sum_{m = 1}^{M} a_m  z^{-m}},
\f]
где \f$a_0\f$ не может быть 0, \f$N=M=\f$`ord`.

\param[in]  b
Указатель на вектор коэффициентов числителя 
передаточной функции \f$H(z)\f$ БИХ-фильтра. \n  
Размер вектора `[ord + 1 x 1]`. \n  \n  

\param[in]  a
Указатель на вектор коэффициентов знаменателя
передаточной функции \f$H(z)\f$ БИХ-фильтра. \n  
Размер вектора `[ord + 1 x 1]`. \n  
Этот указатель может быть `NULL`, тогда фильтрация производится 
без использования рекурсивной части 
(вектор коэффициентов `b` задает КИХ-фильтр). \n \n

\param[in]  ord
Порядок фильтра. Количество коэффициентов числителя  и знаменателя 
передаточной функции \f$H(z)\f$ БИХ-фильтра равно `ord + 1`. \n \n

\param[in]  x
Указатель на вектор отсчетов входного сигнала. \n
Размер вектора `[n x 1]`. \n \n

\param[in]  n
Длина входного сигнала. \n \n 

\param[out] y
Указатель на вектор выходных отсчетов фильтра. \n  
Размер вектора `[n x  1]`. \n  
Память должна быть выделена заранее. \n \n

\return
`RES_OK` Если фильтрация произведена успешно. \n 
В противном случае \ref ERROR_CODE_GROUP "код ошибки". \n 

Пример использования функции `filter_iir`:

\include filter_iir_test.c

На входе цифрового фильтра задан сигнал 
\f$s(t) = \sin(2\pi \cdot 0.05 t) + n(t)\f$, где \f$n(t)\f$ белый гауссовский 
шум, с нулевым средним и единичной дисперсией. \n
Фильтр представляет собой эллиптический ФНЧ 6 порядка. 
Входной сигнал фильтруется данным фильтром, и результат сохраняется в файлы:

\verbatim
dat/s.txt  - исходный зашумленный сигнал
dat/sf.txt - сигнал на выходе фильтра.
\endverbatim

По полученным данным производится построение графиков:

\image html filter_iir_test.png

\author Бахурин Сергей www.dsplib.org
***************************************************************************** */
#endif
int DSPL_API filter_iir(double* b, double* a, int ord,
                        double* x, int n, double* y)
{
    double *buf = NULL;
    double *an  = NULL;
    double *bn  = NULL;
    double  u;
    int     k;
    int     m;
    int     count;

    if(!b || !x || !y)
        return ERROR_PTR;

    if(ord < 1 || n < 1)
        return ERROR_SIZE;

    if(a && a[0]==0.0)
        return ERROR_FILTER_A0;

    count = ord + 1;
    buf = (double*) malloc(count*sizeof(double));
    an =  (double*) malloc(count*sizeof(double));

    memset(buf, 0, count*sizeof(double));

    if(!a)
    {
        memset(an, 0, count*sizeof(double));
        bn = b;
    }
    else
    {
        bn =  (double*) malloc(count*sizeof(double));
        for(k = 0; k < count; k++)
        {
            an[k] = a[k] / a[0];
            bn[k] = b[k] / a[0];
        }
    }

    for(k = 0; k < n; k++)
    {
        for(m = ord; m > 0; m--)
            buf[m] = buf[m-1];
        u = 0.0;
        for(m = ord; m > 0; m--)
            u += buf[m]*an[m];

        buf[0] = x[k] - u;
        y[k] = 0.0;
        for(m = 0; m < count; m++)
            y[k] += buf[m] * bn[m];
    }

    if(buf)
        free(buf);
    if(an)
        free(an);
    if(bn && (bn != b))
        free(bn);
    return RES_OK;
}

