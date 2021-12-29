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
\ingroup FIR_FILTER_DESIGN_GROUP
\fn int DSPL_API  fir_linphase(int ord, double w0, double w1, int filter_type, 
int win_type, double win_param, double* h)
\brief
Function calculates linear-phase FIR filter coefficients by window method

FIR filter transfer function is
\f[
H(z) = \sum_{n = 0}^{ord} h_n z^{-n}.
\f]

\param[in]  ord
Filter order. \n
Number of FIR filter coefficients is `ord+1`. \n 
\n

\param[in]  w0
Normalized cutoff frequency for lowpass and highpass filter,
or left cutoff frequency for bandpass or bandstop filter. \n 
\n

\param[in]  w1
Right normalized cutoff frequency for bandpass or bandstop filter. \n
This parameter is ignored for lowpass or highpass filters. \n 
Frequecny `w1` must be higher than `w0`. \n 
\n

\param[in]  filter_type
Filter type. \n
This parameter can be one of follow: \n
\verbatim
DSPL_FILTER_LPF   - lowpass  filter;
DSPL_FILTER_HPF   - highpass filter;
DSPL_FILTER_BPASS - bandpass filter;
DSPL_FILTER_BSTOP - bandstop filter.
\endverbatim
\n 
\n

\param [in]  win_type
Window function type. \n
This parameter can be one of follow: \n
\verbatim
-------------------------------------------------------------------------
 win_type                    |  Description
-----------------------------|-------------------------------------------
 DSPL_WIN_BARTLETT           | Nonparametric Bartlett window
-----------------------------|-------------------------------------------
 DSPL_WIN_BARTLETT_HANN      | Nonparametric Bartlett-Hann window
-----------------------------|-------------------------------------------
 DSPL_WIN_BLACKMAN           | Nonparametric  Blackman window 
-----------------------------|-------------------------------------------
 DSPL_WIN_BLACKMAN_HARRIS    | Nonparametric Blackman-Harris window
-----------------------------|-------------------------------------------
 DSPL_WIN_BLACKMAN_NUTTALL   | Nonparametric Blackman-Nuttall
-----------------------------|-------------------------------------------
 DSPL_WIN_CHEBY              | Parametric Dolph-Chebyshev window.
                             | Parametr  `win_param` sets sidelobe attenuation 
                             | level in dB.
-----------------------------|-------------------------------------------
 DSPL_WIN_COS                | Nonparametric Cosine window
-----------------------------|-------------------------------------------
 DSPL_WIN_FLAT_TOP           | Nonparametric maxflat window
-----------------------------|-------------------------------------------
 DSPL_WIN_GAUSSIAN           | Nonparametric Gauss window
-----------------------------|-------------------------------------------
 DSPL_WIN_HAMMING            | Nonparametric Hamming window
-----------------------------|-------------------------------------------
 DSPL_WIN_HANN               | Nonparametric Hann window
-----------------------------|-------------------------------------------
 DSPL_WIN_KAISER             | Parametric Kaiser window
-----------------------------|-------------------------------------------
 DSPL_WIN_LANCZOS            | Nonparametric Lanczos window
-----------------------------|-------------------------------------------
 DSPL_WIN_NUTTALL            | Nonparametric Nuttall window
-----------------------------|-------------------------------------------
 DSPL_WIN_RECT               | Nonparametric rectangular window
-------------------------------------------------------------------------
\endverbatim
\n 
\n

\param [in]  win_param
Parameter value for parametric windows. \n
This parameter is used for parametric windows only and is ignored for
nonparametric windows. \n 
\n

\param[out]  h
Pointer to the linear-phase FIR filter coefficients vector. \n
Vector size is `[ord+1 x 1]`. \n
Memoru must be allocated. \n 
\n

\note
Only symmetric windows can achieve linear-phase FIR filter. \n \n
Bandstop filter type (`filter_type = DSPL_FILTER_BSTOP`) requires 
only even filter order `ord`. 
If `filter_type = DSPL_FILTER_BSTOP` and  `ord` is odd then function 
returns `ERROR_FILTER_ORD` code.
\n

\return
`RES_OK` if filter coefficients is calculated successfully. \n 
Else \ref ERROR_CODE_GROUP "code error".


Example:

\include fir_linphase_test.c

This function calculates coeffictiens of lowpass, highpass, bandpass 
 and bandstop linear-phase FIR filters by using different kind of windows.
 Also program calculates filter magnitudes and plots. \n 

\image html fir_linphase_test.png


\author Sergey Bakhurin www.dsplib.org 
***************************************************************************** */
#endif
#ifdef DOXYGEN_RUSSIAN
/*! ****************************************************************************
\ingroup FIR_FILTER_DESIGN_GROUP
\fn int DSPL_API  fir_linphase(int ord, double w0, double w1, int filter_type, 
                               int win_type, double win_param, double* h)
\brief
Расчет коэффициентов линейно-фазового КИХ-фильтра 
методом оконного взвешивания.

Функция рассчитывает коэффициенты передаточной характеристики 
\f[
H(z) = \sum_{n = 0}^{ord} h_n z^{-n}
\f]
цифрового линейно-фазового КИХ-фильтра фильтра. 

\param[in]  ord
Порядок фильтра (количество элементов задержки). \n
Количество коэффициентов фильтра равно `ord+1`. \n 
\n

\param[in]  w0
Нормированная частота среза ФНЧ или ФВЧ, 
или левая частота среза для полосового и режекторного фильтра. \n 
\n

\param[in]  w1
Правая частота среза полосового и режекторного фильтра. \n
Данный параметр игнорируется для ФНЧ и ФВЧ. \n 
Частота `w1` должна быть больше `w0`. \n 
\n

\param[in]  filter_type
Тип фильтра. \n
Данный параметр определяет тип фильтра 
и может принимать одно из значений: \n
\verbatim
DSPL_FILTER_LPF   - фильтр нижних частот;
DSPL_FILTER_HPF   - фильтр верхних частот;
DSPL_FILTER_BPASS - полосовой фильтр;
DSPL_FILTER_BSTOP - режекторный фильтр.
\endverbatim
\n 
\n

\param [in]  win_type
Тип оконной функции. \n
Может принимать одно из следующих значений: \n
\verbatim
-------------------------------------------------------------------------
Значение  win_type           |  Описание
-----------------------------|-------------------------------------------
 DSPL_WIN_BARTLETT           | Непараметрическое окно Бартлетта
-----------------------------|-------------------------------------------
 DSPL_WIN_BARTLETT_HANN      | Непараметрическое окно Бартлетта-Ханна
-----------------------------|-------------------------------------------
 DSPL_WIN_BLACKMAN           | Непараметрическое окно Блэкмана 
-----------------------------|-------------------------------------------
 DSPL_WIN_BLACKMAN_HARRIS    | Непараметрическое окно Блэкмана-Харриса
-----------------------------|-------------------------------------------
 DSPL_WIN_BLACKMAN_NUTTALL   | Непараметрическое окно Блэкмана-Натталла
-----------------------------|-------------------------------------------
 DSPL_WIN_CHEBY              | Параметрическое окно Дольф-Чебышева.
                             | Параметр  win_param  задает уровень
                             | боковых лепестков в дБ.
-----------------------------|-------------------------------------------
 DSPL_WIN_COS                | Непараметрическое косинус-окно
-----------------------------|-------------------------------------------
 DSPL_WIN_FLAT_TOP           | Непараметрическое окно с максимально 
                             | плоской вершиной
-----------------------------|-------------------------------------------
 DSPL_WIN_GAUSSIAN           | Параметрическое окно Гаусса
-----------------------------|-------------------------------------------
 DSPL_WIN_HAMMING            | Непараметрическое окно Хемминга
-----------------------------|-------------------------------------------
 DSPL_WIN_HANN               | Непараметрическое окно Ханна
-----------------------------|-------------------------------------------
 DSPL_WIN_KAISER             | Параметрическое окно Кайзера
-----------------------------|-------------------------------------------
 DSPL_WIN_LANCZOS            | Непараметрическое окно Ланкзоса
-----------------------------|-------------------------------------------
 DSPL_WIN_NUTTALL            | Непараметрическое окно Натталла
-----------------------------|-------------------------------------------
 DSPL_WIN_RECT               | Непараметрическое прямоугольное окно
-------------------------------------------------------------------------
\endverbatim
\n 
\n

\param [in]  win_param
Параметр окна. \n
Данный параметр применяется только для параметрических оконных функций. \n
Для непараметрических окон игнорируется. \n 
\n

\param[out]  h
Указатель на вектор коэффициентов линейно-фазового КИХ-фильтра \f$H(z)\f$. \n
Размер вектора `[ord+1 x 1]`. \n
Память должна быть выделена. \n 
\n

\note
Для соблюдения условия линейной ФЧХ используются 
только симметричные окна. \n \n
Расчет режекторного линейно-фазового КИХ-фильтра 
(если `filter_type = DSPL_FILTER_BSTOP`) производится только 
для фильтров чётного порядка `ord`. 
В случае нечетного порядка `ord` функция вернет код ошибки `ERROR_FILTER_ORD`.
\n

\return
`RES_OK`
Фильтр рассчитан успешно. \n
В противном случае \ref ERROR_CODE_GROUP "код ошибки". \n


Пример использования функции:

\include fir_linphase_test.c

Программа расчитывает коэффициенты и АЧХ линейно-фазовых КИХ-фильтров нижних, 
верхних частот, полосовых и режекторных с применением различных весовых окон: 
прямоугольное, Хемминга, Блэкмана и Блэкмана-Харриса. \n
Полученные АЧХ выводятся на график 

\image html fir_linphase_test.png


\author  Бахурин Сергей www.dsplib.org
***************************************************************************** */
#endif
int DSPL_API fir_linphase(int ord, double w0, double w1, int filter_type,
                          int win_type, double win_param, double* h)
{
    int n, err;
    double wc, b, del;

    if(ord<1)
        return ERROR_FILTER_ORD;
    if(w0 <= 0.0)
        return ERROR_FILTER_WP;
    if(!h)
        return ERROR_PTR;

    switch(filter_type & DSPL_FILTER_TYPE_MASK)
    {
        /* Lowpass FIR coefficients calculation */
        case DSPL_FILTER_LPF:
            err = fir_linphase_lpf(ord, w0, win_type, win_param, h);
            break;

        /* Highpass FIR coefficients calculation */
        case DSPL_FILTER_HPF:
            err = fir_linphase_lpf(ord, 1.0-w0, win_type, win_param, h);
            if(err == RES_OK)
            {
                /* LPF filter frequency inversion */
                for(n = 0; n < ord+1; n+=2)
                    h[n] = -h[n];
            }
            break;

        /* Bandpass FIR coefficients calculation */
        case DSPL_FILTER_BPASS:
            if(w1 < w0)
            {
                err = ERROR_FILTER_WS;
                break;
            }
            wc = (w0 + w1) * 0.5; /* central frequency */
            b  =  w1 - w0;        /* bandwidth */
            err = fir_linphase_lpf(ord, b*0.5, win_type, win_param, h);
            if(err == RES_OK)
            {
                /* LPF frequency shifting to the central frequency */
                del = 0.5 * (double)ord;
                for(n = 0; n < ord+1; n++)
                    h[n] *= 2.0 * cos(M_PI * ((double)n - del) * wc);
            }
            break;

        /* BandStop FIR coefficients calculation */
        /* ATTENTION! Bandstop filter must be even order only! */
        case DSPL_FILTER_BSTOP:
        {
            double *h0 = NULL;

            /* check filter order. Return error if order is odd. */
            if(ord%2)
                return ERROR_FILTER_ORD;

            /* check frequency (w1 must be higher than w0) */
            if(w1 < w0)
            {
                err = ERROR_FILTER_WS;
                break;
            }
            /* temp coeff vector */
            h0 = (double*)malloc((ord+1) * sizeof(double));

            /* calculate LPF */
            err = fir_linphase(ord, w0, 0.0, DSPL_FILTER_LPF,
                               win_type, win_param, h0);
            if(err!=RES_OK)
            {
                free(h0);
                return err;
            }
            /* calculate HPF */
            err = fir_linphase(ord, w1, 0.0, DSPL_FILTER_HPF,
                               win_type, win_param, h);
            if(err==RES_OK)
            {
                /* Bandstop filter is sum of lowpass and highpass filters */
                for(n = 0; n < ord+1; n++)
                    h[n] += h0[n];
            }
            free(h0);
            break;
        }
        default:
            err = ERROR_FILTER_FT;
    }
    return err;
}
