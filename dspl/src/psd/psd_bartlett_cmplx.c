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
\fn int psd_bartlett_cmplx(complex_t* x, int n, int nfft,
                                fft_t* pfft, double fs,
                                int flag, double* ppsd, double* pfrq)
\brief Непараметрическая оценка спектральной плотности мощности (СПМ) 
комплексного сигнала методом Бартлетта.

Функция рассчитывает спектральную плотность мощности \f$ X(f) \f$ 
выборки сигнала длительности \$n \$ отсчетов методом Бартлетта:
\f[
  X(f) = \frac{1}{N F_s }  \sum_{p = 0}^{P-1}\left| \sum_{m = 0}^{n_{FFT}-1} 
   x(m+p \cdot n_{\text{FFT}})   \exp 
   \left( -j 2\pi f m \right) \right|^2,
\f]
где \f$ F_s \f$ -- частота 
дискретизации (Гц), \f$P = n/n_{\text{FFT}}\f$ -- количество сегментов 
смещений выборки сигналов размера \f$n_{FFT}\f$.


При использовании \f$n_{FFT} = n\f$ оценка Бартлетта переходит
в стандартную периодограмму.

Расчет спектральной плотности мощности ведется при помощи алгоритмов быстрого
преобразования Фурье, для дискретной сетки частот от 0 Гц до \f$ F_s \f$ Гц 
(по умолчанию), или от \f$-F_s /2 \f$ до \f$F_s /2 \f$, если установлен флаг 
расчета двусторонней СПМ.

\note Метод Бартлетта возвращает асимптотически несмещенную, 
состоятельную оценку СПМ (уровень флуктуаций шумовой СПМ 
 уменьшается с ростом длины выборки `n` при фиксированной `nfft`).

\param[in]  x
Указатель на входной вектор комплексного сигнала \f$x(m)\f$, 
\f$ m = 0 \ldots n-1 \f$.  \n
Размер вектора `[n x 1]`.  \n \n

\param[in]  n
Размер вектора входного сигнала.
Также размер выходного вектора СПМ и 
вектора частоты также равны `n`.\n\n


\param[in] nfft
Размер сегмента.\n
Размер выходного вектора СПМ, и соответствующего ей вектора частоты.\n\n


\param[in] pfft
Указатель на структуру \ref fft_t.  \n
Указатель может быть `NULL`. В этом случае объект структуры будет 
создан внутри функции и удален перед завершением.\n
Если предполагается многократный вызов функции, то рекомендуется создать 
объект \ref fft_t и передавать в функцию, чтобы не 
создавать его каждый раз. \n\n

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

Пример оценок СПМ методом Бартлетта:

\include psd_bartlett_test_cmplx.c

Программа производит расчет СПМ сигнала, состоящего из двух комплексных 
экспонент на фоне белого гауссова шума. 
Расчет ведется по выборкe длины 8192 отсчета при длине сегмента `nfft` 
128, 1024 и 8192 отсчетов.

Рассчитанные СПМ выводятся на графики:

`n = 8192, nfft = 8192`:
\image html psd_bartlett_cmplx_8192.png

`n = 8192, nfft = 1024`:
\image html psd_bartlett_cmplx_1024.png

`n = 8192, nfft = 128`:
\image html psd_bartlett_cmplx_128.png


Можно видеть, что метод Бартлетта позволяет снизить 
уровень флуктуация шума с увеличением количества сегментов. 
Однако наблюдается эффект растекания спектра, который существенно ухудшает
динамический диапазон анализа. 

Для более качественной оценки СПМ смотри функцию \ref psd_welch_cmplx

\author Бахурин Сергей www.dsplib.org 
***************************************************************************** */
#endif
int DSPL_API psd_bartlett_cmplx(complex_t* x, int n, int nfft,
                                fft_t* pfft, double fs,
                                int flag, double* ppsd, double* pfrq)
{
    int err, pos, k;
    double *pdgr = NULL;
    complex_t *tmp = NULL;
    fft_t *ptr_fft = NULL;
    
    pos = 0;
    
    pdgr = (double*)malloc(nfft * sizeof(double));
    if(!pdgr)
        return ERROR_MALLOC;

    if(!pfft)
    {
        ptr_fft = (fft_t*)malloc(sizeof(fft_t));
        memset(ptr_fft, 0, sizeof(fft_t));
    }
    else
        ptr_fft = pfft;
    
    memset(ppsd, 0, nfft * sizeof(double));
    while(pos + nfft <= n)
    {
        err = fft_mag_cmplx(x + pos, nfft, ptr_fft, fs, 
                            flag & DSPL_FLAG_FFT_SHIFT, pdgr, NULL);
        if(err != RES_OK)
            goto exit_label;
        for(k = 0; k < nfft; k++)
            ppsd[k] += pdgr[k];
        pos += nfft;
    }
    
    if(pos < n)
    {
        tmp = (complex_t*)malloc(nfft * sizeof(complex_t));
        if(!tmp)
        {
            err = ERROR_MALLOC;
            goto exit_label;
        }
        memset(tmp ,0, nfft * sizeof(complex_t));
        memcpy(tmp, x + pos, (n - pos)*sizeof(complex_t));

        err = fft_mag_cmplx(tmp, nfft, ptr_fft, fs, 
                            flag & DSPL_FLAG_FFT_SHIFT, pdgr, NULL);
        if(err != RES_OK)
            goto exit_label;
        
        for(k = 0; k < nfft; k++)
            ppsd[k] += pdgr[k];
    }

    /* fill frequency */
    if(pfrq)
    {
        if(flag & DSPL_FLAG_FFT_SHIFT)
            if(n%2)
                err = linspace(-fs*0.5 + fs*0.5/(double)nfft, 
                                fs*0.5 - fs*0.5/(double)nfft, 
                                n, DSPL_SYMMETRIC, pfrq);
            else
                err = linspace(-fs*0.5, fs*0.5, nfft, DSPL_PERIODIC, pfrq);
        else
            err = linspace(0, fs, nfft, DSPL_PERIODIC, pfrq);
    }
    
    /* scale magnitude */
    if(flag & DSPL_FLAG_LOGMAG)
    {
        for(k = 0; k < nfft; k++)
            ppsd[k] = 10.0 * log10(ppsd[k] / (double)n / fs);
    }
    else
    {
        for(k = 0; k < nfft; k++)
            ppsd[k] /= (double)n * fs;
    }


exit_label:
    if(pdgr)
        free(pdgr);
    if(tmp)
        free(tmp);
    if(ptr_fft && (ptr_fft != pfft))
    {
        fft_free(ptr_fft);
        free(ptr_fft);
    }
    return err;
}


