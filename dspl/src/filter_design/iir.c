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


#ifdef DOXYGEN_ENGLISH
/*! ****************************************************************************
\ingroup IIR_FILTER_DESIGN_GROUP
\fn int iir(double rp, double rs, int ord,    double    w0, double    w1, 
            int type, double* b,    double* a)
\brief
Digital IIR filter design.

The function calculates the coefficients of the digital IIR filter 
transfer fucntion \f$ H(z) \f$. 
Filter coeffitients can be used in \ref filter_iir function

\param[in] rp
Magnitude ripple in passband (dB). \n
\n


\param[in] rs
Suppression level in stopband (dB). \n
\n

\param[in] ord
Filter order. \n
Number of \f$H(z)\f$ numerator and denominator coefficients is `ord+1`. \n
For bandpass and bandstop filters `ord` must be even. \n 
\n

\param[in] w0
Normalized cutoff frequency (from 0 to 1) for lowpass or highpass filter. \n
Or left normalized cutoff frequency (from 0 to 1) for 
bandpass and bandstop filter. \n 
\n


\param[in] w1
Right normalized cutoff frequency (from 0 to 1) for 
bandpass and bandstop filter. \n 
This parameter is ingnored for lowpass and highpass filters.
\n

\param[in] type
Filter type. \n
This patameter sets combination of filter type (one of follow): \n
\verbatim
DSPL_FILTER_LPF   - lowpass    filter;
DSPL_FILTER_HPF   - highpass filter;
DSPL_FILTER_BPASS - bandpass filter;
DSPL_FILTER_BSTOP - bandstop filter,
\endverbatim
and of filter approximation type (one of follow):
\verbatim
DSPL_FILTER_BUTTER - Butterworth filter;
DSPL_FILTER_CHEBY1 - Chebyshev of the first kind filter;
DSPL_FILTER_CHEBY2 - Chebyshev of the second kind filter;
DSPL_FILTER_ELLIP  - Elliptic filter.
\endverbatim
\n 
\n

\param[out] b
Pointer to the transfer function \f$H(z)\f$ 
numerator coefficients vector. \n 
Vector size is `ord+1`. \n
Memory must be allocated. \n 
\n

\param[out] a
Pointer to the transfer function \f$H(z)\f$ 
denominator coefficients vector. \n 
Vector size is `ord+1`. \n
\n

\return
`RES_OK` if filter is calculated successfully. \n 
Else \ref ERROR_CODE_GROUP "code error".

Example:

\include iir_test.c

This program calcultes filter coefficients for different flags `type`.

In addition, the filters magnitudes 
is calculated and plotted by GNUPLOT package.

\image html iir_test.png

\author Sergey Bakhurin www.dsplib.org
***************************************************************************** */
#endif
#ifdef DOXYGEN_RUSSIAN
/*! ****************************************************************************
\ingroup IIR_FILTER_DESIGN_GROUP
\fn int iir(double rp, double rs, int ord, double w0, double w1, 
            int type, double* b, double* a)
\brief
Функция расчета коэффициентов передаточной характеристики \f$H(z)\f$
цифрового фильтра БИХ.

Функция рассчитывает коэффициенты передаточной характеристики \f$H(z)\f$
цифрового фильтра, которые могут быть использованы в функции \ref filter_iir 

\param[in] rp
Уровень неравномерности квадрата АЧХ в полосе пропускания фильтра (дБ). \n
Размер вектора `[ord+1 x 1]`. \n 
\n


\param[in] rs
Уровень подавления в полосе заграждения фильтра (дБ).\n 
\n

\param[in] ord
Порядок фильтра. \n
Количество коэффициентов числителя и знаменателя передаточной 
функции \f$H(z)\f$ цифрового фильтров равно `ord+1`. \n
Для полосовых и режекторных фильтров параметр `ord` должен быть чётным. \n 
\n

\param[in] w0
Нормированная частота среза ФНЧ или ФВЧ, или левая частота среза для 
полосового и режекторного фильтра.\n 
\n


\param[in] w1
Правая частота среза полосового и режекторного фильтра. \n
Данный параметр игнорируется для ФНЧ и ФВЧ. \n 
\n

\param[in] type
Тип фильтра. \n
Данный параметр определяет тип фильтра и образуется 
набором флагов типа фильтра: \n
\verbatim
DSPL_FILTER_LPF - фильтр нижних частот;
DSPL_FILTER_HPF - фильтр верхних частот;
DSPL_FILTER_BPASS - полосовой фильтр;
DSPL_FILTER_BSTOP - режекторный фильтр,
\endverbatim
а также флагов типа аппроксимации АЧХ фильтра:
\verbatim
DSPL_FILTER_BUTTER - фильтр Баттерворта;
DSPL_FILTER_CHEBY1 - фильтр Чебышева первого рода;
DSPL_FILTER_CHEBY2 - фильтр Чебышева второго рода;
DSPL_FILTER_ELLIP - эллиптический фильтр.
\endverbatim
\n 
\n

\param[out] b
Указатель на вектор коэффициентов 
числителя передаточной функции \f$H(z)\f$. \n 
Размер вектора    `ord+1`. \n
Память должна быть выделена. \n 
\n

\param[out] a
Указатель на вектор коэффициентов знаменателя передаточной 
функции \f$H(z)\f$. \n 
Размер вектора    `ord+1`. \n
Память должна быть выделена. \n 
\n

\return
`RES_OK` --- Фильтр рассчитан успешно. \n
В противном случае \ref ERROR_CODE_GROUP "код ошибки". \n

Пример использования функции:

\include iir_test.c

Данная программа производит расчет коэффициентов фильтров 
при различном сочетании флагов параметра `type`.

Кроме этого производится расчет АЧХ полученных цифровых фильтров и выводится на 
график АЧХ пакетом GNUPLOT

\image html iir_test.png

\author Бахурин Сергей www.dsplib.org
***************************************************************************** */
#endif
int DSPL_API iir(double rp, double rs, int ord, double w0, double w1,
                 int type, double* b,    double* a)
{
    double *bs = NULL;
    double *as = NULL;
    double *bt = NULL;
    double *at = NULL;
    double wa0, wa1, ws;
    int err, ord_ap = ord;
    int i;
    
    if(((type & DSPL_FILTER_TYPE_MASK) == DSPL_FILTER_LPF) ||
         ((type & DSPL_FILTER_TYPE_MASK) == DSPL_FILTER_HPF))
    {
        bs = (double*)malloc((ord_ap+1)*sizeof(double));
        as = (double*)malloc((ord_ap+1)*sizeof(double));
        bt = (double*)malloc((ord_ap+1)*sizeof(double));
        at = (double*)malloc((ord_ap+1)*sizeof(double));
    }


    if(((type & DSPL_FILTER_TYPE_MASK) == DSPL_FILTER_BPASS) ||
         ((type & DSPL_FILTER_TYPE_MASK) == DSPL_FILTER_BSTOP))
    {
        if(ord % 2)
            return ERROR_FILTER_ORD_BP;
        else
        {
            ord_ap = ord / 2;
            bs = (double*)malloc((ord_ap + 1)*sizeof(double));
            as = (double*)malloc((ord_ap + 1)*sizeof(double));
            bt = (double*)malloc((ord        + 1)*sizeof(double));
            at = (double*)malloc((ord        + 1)*sizeof(double));
        }
    }
    err = iir_ap(rp, rs, ord_ap, type, bs, as);
    if(err != RES_OK)
        goto error_proc;

    /* frequency transformation    */
    wa0 = tan(w0 * M_PI * 0.5);
    wa1 = tan(w1 * M_PI * 0.5);

    switch(type & DSPL_FILTER_TYPE_MASK)
    {

        case DSPL_FILTER_LPF:
            err = low2low(bs, as, ord_ap, 1.0, wa0, bt, at);
            break;

        case DSPL_FILTER_HPF:
            ws    = filter_ws1(ord_ap, rp, rs, type);
            err = low2low( bs, as, ord_ap, 1.0, 1.0 / ws,    bs, as);
            err = low2high(bs, as, ord_ap, 1.0, wa0, bt, at);
            break;

        case DSPL_FILTER_BPASS:
            err = low2bp(bs, as, ord_ap, 1.0, wa0, wa1, bt, at);
            break;

        case DSPL_FILTER_BSTOP:
            /* need frequency transform ws ->    1    rad/s     */

            ws    = filter_ws1(ord_ap, rp, rs, type);
            err = low2low( bs, as, ord_ap, 1.0, 1.0 / ws,    bs, as);
            err = low2bs(bs, as, ord_ap, 1.0, wa0, wa1, bt, at);
            break;

        default:
            err = ERROR_FILTER_TYPE;
            break;
    }
    if(err != RES_OK)
        goto error_proc;


    err = bilinear(bt, at, ord, b, a);
    
    for(i =  1; i <= ord; i++)
    {
        a[i] /= a[0];
        b[i] /= a[0];
    }
    b[0] /= a[0];
    a[0] = 1.0;

error_proc:

    if(bs)
        free(bs);
    if(as)
        free(as);
     if(bt)
        free(bt);
    if(at)
        free(at);

    return err;

}