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
#include <stdio.h>
#include <string.h>
#include <float.h>
#include "dspl.h"






#ifdef DOXYGEN_ENGLISH

#endif
#ifdef DOXYGEN_RUSSIAN

#endif
int DSPL_API psd_bartlett(double* x, int n, int nfft,
                          fft_t* pfft, double fs,
                          int flag, double* ppsd, double* pfrq)
{
    int err, pos, k;
    double *pdgr = NULL;
    double *tmp = NULL;
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
        err = fft_mag(x + pos, nfft, ptr_fft, fs, 
                      flag & DSPL_FLAG_FFT_SHIFT, pdgr, NULL);
        if(err != RES_OK)
            goto exit_label;
        for(k = 0; k < nfft; k++)
            ppsd[k] += pdgr[k];
        pos += nfft;
    }
    
    if(pos < n)
    {
        tmp = (double*)malloc(nfft * sizeof(double));
        if(!tmp)
        {
            err = ERROR_MALLOC;
            goto exit_label;
        }
        memset(tmp ,0, nfft * sizeof(double));
        memcpy(tmp, x + pos, (n - pos)*sizeof(double));

        err = fft_mag(tmp, nfft, ptr_fft, fs, 
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





#ifdef DOXYGEN_ENGLISH

#endif
#ifdef DOXYGEN_RUSSIAN
/*! ****************************************************************************
\ingroup PSD_GROUP
\fn 
int psd_periodogram(double* x, int n,
                             int win_type, double win_param,
                             fft_t* pfft, double fs,
                             int flag, double* ppsd, double* pfrq)
\brief Непараметрическая оценка спектральной плотности сигнала методом 
       модифицированной периодограммы

Функция рассчитывает спектральную плотность мощности \f$ X(f) \f$ 
выборки сигнала длительности \$n \$ отсчетов методом модифицированной 
периодограммы:

\f[
  X(f) = \frac{1}{N U F_s} \Big| \sum_{m = 0}^{n-1} w(m) x(m)   \exp 
   \left( -j \frac{2\pi} f m \right) \Big|^2,
\f]
где \f$ w(m) \f$ -- отсчёты оконной функции, \f$ F_s \f$ -- частота 
дискретизации (Гц), \f$ U \f$ нормировочный коэффициент равный

\f[
    U \sum_{m = 0}^{n-1} w^2(m)
\f]

При использовании прямоугольного окна модифицированная периодограмма переходит
в стандартную периодограмму.

Расчет спектральной плотности мощности ведется при помощи алгоритмов быстрого
преобразования Фурье, для дискретной сетки частот от 0 Гц до \f$ F_s \f$ Гц 
(по умолчанию), или от \f$-F_s /2 \f$ до \f$F_s /2 \f$, если установлен флаг 
расчета двусторонней периодограммы.


\note Периодограмма, как стандартная, так и модифицированная возвращает 
асимптотически несмещенную, но несостоятельную оценку СПМ (уроверь флуктуаций
шумовой составляющей СПМ не уменьшается с ростом длины выборки `n`).



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


\param[in] win_type
Параметр оконной функции.\n 
Данный параметр используется, если задано параметрическая оконная функция.
Для непараметрических окон данный параметр игнорируется.\n
Подробнее смотри описание функции \ref window. \n\n

\param[in] pfft
Указатель на структуру `fft_t`.  \n
Указатель может быть `NULL`. В этом случае объект структуры будет 
создан внутри функции и удален перед завершением.\n
Если предполагается многократный вызов функции, то рекомендуется создать 
объект `fft_t` и передавать в функцию, чтобы не создавать его каждый раз. \n\n

\param[in] fs
частота дискретизации выборки исходного сигнала (Гц). \n\n

\param[in] flag
Комбинация битовых флагов, задающих режим расчета:
\verbatim
DSPL_FLAG_LOGMAG - СПМ считать в логарифмическом масштабе в единицах дБ/Гц
DSPL_FLAG_PSD_TWOSIDED - двусторонняя СПМ (от -Fs/2 до Fs/2)
\endverbatim

\param[in, out] ppsd
Указатель на вектор СПМ расчитанных по входному сигналу $x$. \n
Размер вектора `[n x 1]`. \n
Память должна быть выделена. \n\n

\param[in, out] pfrq
Указатель на вектор частоты, соответсвующей 
значениям рассчитанного вектора СПМ. \n
Размер вектора `[n x 1]`. \n
Указатель можеть быть `NULL`,в этом случае вектор частоты не 
рассчитывается и не возвращается. \n\n


\return
`RES_OK` если расчет произведен успешно.  \n
В противном случае \ref ERROR_CODE_GROUP "код ошибки". \n \n

\author Бахурин Сергей www.dsplib.org 
***************************************************************************** */

#endif
int DSPL_API psd_periodogram(double* x, int n,
                             int win_type, double win_param,
                             fft_t* pfft, double fs,
                             int flag, double* ppsd, double* pfrq)
{
    double *w = NULL;
    double *s = NULL;
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
        s = (double*)malloc(n*sizeof(double));
        if(!s)
        {
            err = ERROR_MALLOC;
            goto exit_label;
        }
        
        /* windowing */
        for(k = 0; k < n; k++)
            s[k] = x[k] * w[k];
    }
    else
    {
        /* classic periodogram without windowing */
        s = x;
        wn = (double)n;
    }

    /* calculate FFT */
    err = fft_mag(s, n, ptr_fft, fs, flag, ppsd, pfrq);
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



#ifdef DOXYGEN_ENGLISH

#endif
#ifdef DOXYGEN_RUSSIAN

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




#ifdef DOXYGEN_ENGLISH

#endif
#ifdef DOXYGEN_RUSSIAN

#endif
int DSPL_API psd_welch(double* x, int n, 
                      int win_type, double win_param,
                      int npsd, int noverlap, fft_t* pfft, double fs,
                      int flag, double* ppsd, double* pfrq)
{
    double *win = NULL;
    double wg;
    complex_t *tmp = NULL;
    double *s = NULL;
    int err, k, pos, cnt;
    fft_t *ptr_fft = NULL;

    if(!x || !ppsd)
        return  ERROR_PTR;

    if(n<1 || npsd < 1)
        return ERROR_SIZE;

    if(noverlap < 1 || noverlap > npsd)
        return ERROR_OVERLAP;

    if(fs < 0.0)
        return ERROR_FS;

    win = (double*)malloc(npsd * sizeof(double));
    if(!win)
    {
        err = ERROR_MALLOC;
        goto exit_label;
    }

    if(!pfft)
    {
        ptr_fft = (fft_t*)malloc(sizeof(fft_t));
        memset(ptr_fft, 0, sizeof(fft_t));
    }
    else
        ptr_fft = pfft;
      

    err = window(win, npsd, win_type, win_param);
    if(err != RES_OK)
        goto exit_label;

    wg = 0.0;
    for(k = 0; k < npsd; k++)
        wg += win[k] * win[k];
    wg = 1.0 / wg;
    
    
    tmp = (complex_t*)malloc(npsd*sizeof(complex_t));
    if(!tmp)
    {
        err = ERROR_MALLOC;
        goto exit_label;
    }

    s = (double*)malloc(npsd*sizeof(double));
    if(!s)
    {
        err = ERROR_MALLOC;
        goto exit_label;
    }

    pos = 0;
    cnt = 0;
    memset(ppsd, 0, npsd*sizeof(double));
    while(pos+npsd <= n)
    {
        for(k = 0; k < npsd; k++)
            s[k] = x[k+pos]*win[k];

        err = fft(s, npsd, ptr_fft, tmp);
        if(err != RES_OK)
            goto exit_label;

        for(k = 0; k < npsd; k++)
            ppsd[k] += wg * ABSSQR(tmp[k]);

        pos += noverlap;
        cnt++;
        
    }

    for(k = 0; k < npsd; k++)
        ppsd[k] /= (double)cnt * fs;

    if(flag & DSPL_FLAG_LOGMAG)
    {
        for(k = 0; k < npsd; k++)
            ppsd[k] = 10.0 * log10(ppsd[k] + DBL_EPSILON);
    }
    


    if(flag & DSPL_FLAG_PSD_TWOSIDED)
    {
        err = fft_shift(ppsd, npsd, ppsd);
        if(err != RES_OK)
            goto exit_label;
    }

    if(pfrq)
    {
        if(flag & DSPL_FLAG_PSD_TWOSIDED)
        {
            err = linspace(-0.5*fs, fs*0.5, npsd, DSPL_PERIODIC, pfrq);
            if(err != RES_OK)
                goto exit_label;
        }
        else
        {
            err = linspace(0, fs, npsd, DSPL_PERIODIC, pfrq);
            if(err != RES_OK)
                goto exit_label;
        }
    }

    err = RES_OK;
exit_label:
    if(win)
        free(win);
    if(tmp)
        free(tmp);
    if(s)
        free(s);
    if(ptr_fft && (ptr_fft != pfft))
        free(ptr_fft);
    return err;
}







#ifdef DOXYGEN_ENGLISH

#endif
#ifdef DOXYGEN_RUSSIAN

#endif
int DSPL_API psd_welch_cmplx(complex_t* x, int n, 
                             int win_type, double win_param,
                             int npsd, int noverlap, fft_t* pfft, double fs,
                             int flag, double* ppsd, double* pfrq)
{
    double    *win     = NULL;
    complex_t *tmp     = NULL;
    complex_t *s       = NULL;
    fft_t     *ptr_fft = NULL;
    
    double wg;
    int err, k, pos, cnt;
    
    if(!x || !ppsd)
        return  ERROR_PTR;

    if(n<1 || npsd < 1)
        return ERROR_SIZE;

    if(noverlap < 1 || noverlap > npsd)
        return ERROR_OVERLAP;

    if(fs < 0.0)
        return ERROR_FS;

    win = (double*)malloc(npsd * sizeof(double));
    if(!win)
    {
        err = ERROR_MALLOC;
        goto exit_label;
    }

    if(!pfft)
    {
        ptr_fft = (fft_t*)malloc(sizeof(fft_t));
        memset(ptr_fft, 0, sizeof(fft_t));
    }
    else
        ptr_fft = pfft;

    err = window(win, npsd, win_type, win_param);
    if(err != RES_OK)
        goto exit_label;

    wg = 0.0;
    for(k = 0; k < npsd; k++)
        wg += win[k] * win[k];
    wg = 1.0 / wg;

    tmp = (complex_t*)malloc(npsd*sizeof(complex_t));
    if(!tmp)
    {
        err = ERROR_MALLOC;
        goto exit_label;
    }

    s = (complex_t*)malloc(npsd*sizeof(complex_t));
    if(!s)
    {
        err = ERROR_MALLOC;
        goto exit_label;
    }

    pos = 0;
    cnt = 0;
    memset(ppsd, 0, npsd*sizeof(double));
    while(pos+npsd <= n)
    {
        for(k = 0; k < npsd; k++)
        {
            RE(s[k]) = RE(x[k+pos])*win[k];
            IM(s[k]) = IM(x[k+pos])*win[k];
        }
        err = fft_cmplx(s, npsd, ptr_fft, tmp);
        if(err != RES_OK)
            goto exit_label;

        for(k = 0; k < npsd; k++)
            ppsd[k] += wg * ABSSQR(tmp[k]);

        pos += noverlap;
        cnt++;
        
    }

    for(k = 0; k < npsd; k++)
        ppsd[k] /= (double)cnt * fs;

    if(flag & DSPL_FLAG_LOGMAG)
    {
        for(k = 0; k < npsd; k++)
            ppsd[k] = 10.0 * log10(ppsd[k] + DBL_EPSILON);
    }

    if(flag & DSPL_FLAG_PSD_TWOSIDED)
    {
        err = fft_shift(ppsd, npsd, ppsd);
        if(err != RES_OK)
            goto exit_label;
    }

    if(pfrq)
    {
        if(flag & DSPL_FLAG_PSD_TWOSIDED)
        {
            err = linspace(-0.5*fs, fs*0.5, npsd, DSPL_PERIODIC, pfrq);
            if(err != RES_OK)
                goto exit_label;
        }
        else
        {
            err = linspace(0, fs, npsd, DSPL_PERIODIC, pfrq);
            if(err != RES_OK)
                goto exit_label;
        }
    }
    err = RES_OK;

exit_label:
    if(win)
        free(win);
    if(tmp)
        free(tmp);
    if(s)
        free(s);
    if(ptr_fft && (ptr_fft != pfft))
        free(ptr_fft);
    return err;
}
