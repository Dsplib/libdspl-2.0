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

#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <float.h>
#include "dspl.h"



#ifdef DOXYGEN_ENGLISH

#endif
#ifdef DOXYGEN_RUSSIAN
/*! ****************************************************************************
\ingroup PSD_GROUP
\fn int psd_periodogram_cmplx(complex_t* x, int n,
                              int win_type, double win_param,
                              fft_t* pfft, double fs,
                              int flag, double* ppsd, double* pfrq)

\brief Непараметрическая оценка спектральной плотности мощности (СПМ) 
комплексного сигнала методом модифицированной периодограммы.

Функция рассчитывает спектральную плотность мощности \f$ X(f) \f$ 
выборки сигнала длительности \$n \$ отсчетов методом модифицированной 
периодограммы:

\f[
  X(f) = \frac{1}{U F_s} \left| \sum_{m = 0}^{n-1} w(m) x(m)   \exp 
   \left( -j 2\pi f m \right) \right|^2,
\f]
где \f$ w(m) \f$ -- отсчёты оконной функции, \f$ F_s \f$ -- частота 
дискретизации (Гц), \f$ U \f$ нормировочный коэффициент равный
\f[
    U = \sum_{m = 0}^{n-1} w^2(m)
\f]

При использовании прямоугольного окна модифицированная периодограмма переходит
в стандартную периодограмму.

Расчет спектральной плотности мощности ведется при помощи алгоритмов быстрого
преобразования Фурье, для дискретной сетки частот от 0 Гц до \f$ F_s \f$ Гц 
(по умолчанию), или от \f$-F_s /2 \f$ до \f$F_s /2 \f$, если установлен флаг 
расчета двусторонней периодограммы.

\note Периодограмма возвращает асимптотически несмещенную, 
но несостоятельную оценку СПМ (уровень флуктуаций шумовой составляющей СПМ 
не уменьшается с ростом длины выборки `n`).

\param[in]  x
Указатель на входной вектор вещественного сигнала \f$x(m)\f$, 
\f$ m = 0 \ldots n-1 \f$.  \n
Размер вектора `[n x 1]`.  \n \n

\param[in]  n
Размер вектора входного сигнала.
Также размер выходного вектора СПМ и 
вектора частоты также равны `n`.\n\n

\param[in]  win_type
Тип оконной функции, применяемой для модифицированной периодограммы.\n
Подробнее смотри описание функции \ref window. \n\n

\param[in]  win_param
Параметр оконной функциии. \n
Данный параметр применяется только для парамтрических типов окон
(смотри описание функции \ref window).\n
Для непараметрических функций игнорируется. \n\n

\param[in] pfft
Указатель на структуру \ref fft_t.  \n
Указатель может быть `NULL`. В этом случае объект структуры будет 
создан внутри функции и удален перед завершением.\n
Если предполагается многократный вызов функции, то рекомендуется создать 
объект \ref fft_t и передавать в функцию, чтобы не создавать его каждый раз. \n\n

\param[in] fs
частота дискретизации выборки исходного сигнала (Гц). \n\n

\param[in] flag
Комбинация битовых флагов, задающих режим расчета:
\verbatim
DSPL_FLAG_LOGMAG - СПМ считать в логарифмическом масштабе в единицах дБ/Гц
DSPL_FLAG_PSD_TWOSIDED - двусторонняя СПМ (от -Fs/2 до Fs/2)
\endverbatim

\param[in, out] ppsd
Указатель на вектор СПМ рассчитанных по входному сигналу $x$. \n
Размер вектора `[n x 1]`. \n
Память должна быть выделена. \n\n

\param[in, out] pfrq
Указатель на вектор частоты, соответствующей 
значениям рассчитанного вектора СПМ. \n
Размер вектора `[n x 1]`. \n
Указатель может быть `NULL`,в этом случае вектор частоты не 
рассчитывается и не возвращается. \n\n

\return
`RES_OK` если расчет произведен успешно.  \n
В противном случае \ref ERROR_CODE_GROUP "код ошибки". \n \n

Пример периодограммных оценок СПМ для различной длины выборки сигнала:

\include psd_periodogram_cmplx_test.c

Программа производит расчет СПМ сигнала, состоящего из двух гармоник на 
фоне белого гауссова шума. Расчет ведется по выборкам длины 128, 1024 и
8192 отсчетов.

В результате периодограммы (стандартная с прямоугольным окном 
и модифицированная с окном Блэкмана) выводятся на графики:

`n = 8192` точек (черная -- классическая периодограмма с прямоугольным окном,
зеленая -- модифицированная с окном Блэкмана):
\image html psd_perodogram_cmplx_8192.png

`n = 1024` точек (черная -- классическая периодограмма с прямоугольным окном,
зеленая -- модифицированная с окном Блэкмана):
\image html psd_perodogram_cmplx_1024.png

`n = 128` точек (черная -- классическая периодограмма с прямоугольным окном,
зеленая -- модифицированная с окном Блэкмана):
\image html psd_perodogram_cmplx_128.png

Можно видеть, что модифицированная периодограмма позволяет снизить 
растекание СПМ, однако уровень флуктуация шума не уменьшается с увеличением 
размера выборки от 128 до 8192 отсчетов (оценка несостоятельная).

\author Бахурин Сергей www.dsplib.org 
***************************************************************************** */
#endif
int DSPL_API psd_periodogram_cmplx(complex_t* x, int n,
                                   int win_type, double win_param,
                                   fft_t* pfft, double fs,
                                   int flag, double* ppsd, double* pfrq)
{
    double *w = NULL;
    complex_t *s = NULL;
    double u, wn;
    int err, k;
    fft_t *ptr_fft = NULL;
    
    if(!x || !ppsd)
        return  ERROR_PTR;

    if(n<1 )
        return ERROR_SIZE;

    if(fs < 0.0)
        return ERROR_FS;
      
    if(!pfft)
    {
        ptr_fft = (fft_t*)malloc(sizeof(fft_t));
        memset(ptr_fft, 0, sizeof(fft_t));
    }
    else
        ptr_fft = pfft;
      
      
    if(win_type != DSPL_WIN_RECT)
    {
        /* Modified periodogram calculation */
        
        /* window malloc */
        w = (double*)malloc(n*sizeof(double));
        if(!w)
        {
            err = ERROR_MALLOC;
            goto exit_label;
        }
        
        /* create window */
        err = window(w, n, win_type, win_param);
        if(err != RES_OK)
            goto exit_label;
        
        /* window normalization wn = sum(w.^2) */
        wn = 0; 
        for(k = 0; k < n; k++)
            wn += w[k]*w[k];
        
        /* signal buffer malloc */
        s = (complex_t*)malloc(n*sizeof(complex_t));
        if(!s)
        {
            err = ERROR_MALLOC;
            goto exit_label;
        }
        
        /* windowing */
        for(k = 0; k < n; k++)
        {
            RE(s[k]) = RE(x[k]) * w[k];
            IM(s[k]) = IM(x[k]) * w[k];
        }
    }
    else
    {
        /* classic periodogram without windowing */
        s = x;
        wn = (double)n;
    }

    /* calculate FFT */
    err = fft_mag_cmplx(s, n, ptr_fft, fs, flag, ppsd, pfrq);
    if(err != RES_OK)
        goto exit_label;
      

    if(flag & DSPL_FLAG_LOGMAG)
    {
        /* normalization in log scale */
        u = 10.0 * log10(wn * fs); 
        for(k = 0; k < n; k++)
            ppsd[k] -= u;
    }
    else
    {
        /* normalization in linear scale */
        u = 1.0 / (wn * fs);
        for(k = 0; k < n; k++)
            ppsd[k] *= u;
    }

exit_label:
    if(w)
        free(w);
    if(s && s != x)
        free(s);
    if(ptr_fft && (ptr_fft != pfft))
    {
        fft_free(ptr_fft);
        free(ptr_fft);
    }
    return err;
}

