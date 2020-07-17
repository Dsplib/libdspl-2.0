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
\ingroup FILTER_ANALYSIS_GROUP
\fn int filter_freq_resp(double* b, double* a, int ord, double* w, int n, 
        int flag, double* mag, double* phi, double* tau)

\brief
Magnitude, phase response and group delay vectors calculation 
for  digital or analog filter corresponds to \f$H(s)\f$, or \f$H(z)\f$
transfer function.


\param[in]  b
Pointer to the \f$ H(s) \f$ or \f$H(z)\f$  transfer function 
numerator coefficients vector. \n 
Vector size is `[ord+1 x 1]`. \n \n 

\param[in]  a
Pointer to the \f$ H(s) \f$ or \f$H(z)\f$ transfer function 
denominator coefficients vector. \n 
Vector size is `[ord+1 x 1]`. \n \n 

\param[in]  ord
Filter order. \n
Transfer function \f$ H(s) \f$ or \f$H(z)\f$ numerator 
and denominator coefficients number equals `ord+1`. \n \n 

\param[in]  w
Pointer to the angular frequency  \f$ \omega \f$ (rad/s), 
which used for analog filter characteristics calculation
(flag sets as `DSPL_FLAG_ANALOG`). \n
For digital filter (flag sets as `DSPL_FLAG_DIGITAL`),
 parameter `w` describes normalized frequency of
frequency response \f$ H \left(\mathrm{e}^{j\omega} \right) \f$.
Digital filter frequency response is \f$ 2\pi \f$-periodic function,
and vector `w` advisable to set from 0 to \f$ \pi \f$,
or from 0 to \f$ 2\pi \f$, or from \f$ -\pi \f$ to \f$ \pi \f$.
Vector size is `[n x 1]`. \n \n

\param[in]  n
Size of frequency vector `w`. \n  \n 

\param[in]  flag
Binary flags to set calculation rules: \n 
\verbatim
DSPL_FLAG_ANALOG  Coefficients corresponds to analog filter
DSPL_FLAG_DIGITAL Coefficients corresponds to digital filter
DSPL_FLAG_LOGMAG  Calculate magnitude in logarithmic scale (in dB) 
DSPL_FLAG_UNWRAP  Unwrap radian phases by adding multiples of 2*pi 
\endverbatim

\param[out]  mag
Pointer to the filter magnitude vector. \n 
Vector size is `[n x 1]`. \n
If pointer is `NULL`, then magnitude will not calculted. \n \n

\param[out]  phi
Pointer to the phase response vector. \n 
Vector size is  `[n x 1]`. \n
If pointer is `NULL`, then phase response will not calculted. \n \n

\param[out]  tau
Pointer to the group delay vector. \n 
Vector size is  `[n x 1]`. \n
If pointer is `NULL`, then group delay will not calculted. \n \n

\return
\return `RES_OK` if function is calculated successfully. \n
Else \ref ERROR_CODE_GROUP "code error".

Example:

\include butter_ap_test.c

Result:

\verbatim
b[ 0] =   1.002   a[ 0] =   1.002
b[ 1] =   0.000   a[ 1] =   2.618
b[ 2] =   0.000   a[ 2] =   3.418
b[ 3] =   0.000   a[ 3] =   2.615
b[ 4] =   0.000   a[ 4] =   1.000
\endverbatim
\n

In `dat` folder will be created 3 files: \n

\verbatim
butter_ap_test_mag.txt    magnitude
butter_ap_test_phi.txt    phase response
butter_ap_test_tau.txt    group delay
\endverbatim

In addition, GNUPLOT will build the following graphs from data stored in files:

\image html butter_ap_test.png

\author Sergey Bakhurin www.dsplib.org
***************************************************************************** */

#endif
#ifdef DOXYGEN_RUSSIAN
/*! ****************************************************************************
\ingroup FILTER_ANALYSIS_GROUP
\fn int filter_freq_resp(double* b, double* a, int ord, double* w, int n, 
        int flag, double* mag, double* phi, double* tau)

\brief
Расчет амплитудно-частотной (АЧХ), фазочастотной характеристик (ФЧХ), а также 
группового времени запаздывания (ГВЗ) цифрового или аналогового или фильтра.

Функция рассчитывает АЧХ, ФЧХ и ГВЗ аналогового или цифрового фильтра, заданного 
передаточной характеристикой \f$H(s)\f$, или \f$H(z)\f$ соответственно

\param[in]  b
Указатель на вектор коэффициентов числителя 
передаточной функции \f$ H(s) \f$. \n 
Размер вектора `[ord+1 x 1]`. \n \n 

\param[in]  a
Указатель на вектор коэффициентов знаменателя 
передаточной функции \f$ H(s) \f$. \n 
Размер вектора `[ord+1 x 1]`. \n \n 

\param[in]  ord
Порядок фильтра. Количество коэффициентов 
числителя и знаменателя передаточной
функции \f$ H(s) \f$ равно `ord+1`. \n  \n 

\param[in]  w
Указатель на вектор значений циклической частоты  \f$ \omega \f$ (рад/с), 
для которого будет рассчитаны  АЧХ, ФЧХ и ГВЗ аналогового фильтра, 
если установлен  флаг `DSPL_FLAG_ANALOG`. \n
В случае если флаг `DSPL_FLAG_ANALOG` не установлен, то  вектор частоты `w` 
используется как нормированная частота комплексного коэффициента передачи 
\f$ H \left(\mathrm{e}^{j\omega} \right) \f$  цифрового фильтра. \n
В этом случае характеристика цифрового фильтра является 
\f$ 2\pi \f$-периодической, и вектор частоты может содержать 
произвольные значения, однако целесообразно задавать 
его от 0 до \f$ \pi \f$, а такжет от 0 до \f$ 2\pi \f$, или 
от \f$ -\pi \f$ до \f$ \pi \f$. \n
Размер вектора `[n x 1]`. \n \n 
 
\param[in]  n
Размер вектора циклической частоты `w`. \n  \n 

\param[in]  flag
Комбинация флагов, которые задают расчет параметров: \n 
\verbatim
DSPL_FLAG_ANALOG  Коэффициенты относятся к аналоговому фильтру
DSPL_FLAG_LOGMAG  АЧХ рассчитывать в логарифмическом масштабе
DSPL_FLAG_UNWRAP  раскрывать периодичность ФЧХ
\endverbatim

\param[out]  mag
Указатель на вектор АЧХ. \n 
Размер вектора `[n x 1]`. \n
Память должна быть выделена. \n
Если указатель `NULL`, то расчет АЧХ не производится. \n \n

\param[out]  phi
Указатель на вектор ФЧХ. \n 
Размер вектора `[n x 1]`. \n
Память должна быть выделена. \n
Если указатель `NULL`, то расчет ФЧХ не производится. \n \n

\param[out]  tau
Указатель на вектор ГВЗ. \n 
Размер вектора `[n x 1]`. \n
Память должна быть выделена. \n
Если указатель `NULL`, то расчет ГВЗ не производится. \n \n

\return
`RES_OK`  Параметры фильтра рассчитаны успешно. \n
В противном случае \ref ERROR_CODE_GROUP "код ошибки". \n

Пример использования функции `filter_freq_resp`:

\include butter_ap_test.c

Результат работы программы:

\verbatim
b[ 0] =   1.002   a[ 0] =   1.002
b[ 1] =   0.000   a[ 1] =   2.618
b[ 2] =   0.000   a[ 2] =   3.418
b[ 3] =   0.000   a[ 3] =   2.615
b[ 4] =   0.000   a[ 4] =   1.000
\endverbatim
\n

В каталоге `dat` будут созданы три файла: \n

\verbatim
butter_ap_test_mag.txt    АЧХ фильтра
butter_ap_test_phi.txt    ФЧХ фильтра
butter_ap_test_tau.txt    ГВЗ фильтра
\endverbatim

Кроме того программа GNUPLOT произведет построение следующих графиков 
по сохраненным в файлах данным:

\image html butter_ap_test.png

\author Бахурин Сергей www.dsplib.org
***************************************************************************** */
#endif
int DSPL_API filter_freq_resp(double* b, double* a, int ord,
                              double* w, int n, int flag,
                              double* mag, double* phi, double* tau)
{
    int res, k, flag_analog;

    complex_t *hc = NULL;
    double *phi0 = NULL;
    double *phi1 = NULL;
    double *w0     = NULL;
    double *w1     = NULL;

    if(!b || !w)
        return ERROR_PTR;
    if(ord < 1)
        return ERROR_FILTER_ORD;
    if(n < 1)
        return ERROR_SIZE;

    flag_analog = flag & DSPL_FLAG_ANALOG;
    
    hc = (complex_t*) malloc (n*sizeof(complex_t));
    
    res = flag_analog ? 
          freqs(b, a, ord, w, n, hc) : 
          freqz(b, a, ord, w, n, hc);

    if(res != RES_OK)
        goto exit_label;


    if(mag)
    {
        if(flag & DSPL_FLAG_LOGMAG)
        {
            for(k = 0; k < n; k++)
                mag[k] = 10.0 * log10(ABSSQR(hc[k]));
        }
        else
        {
            for(k = 0; k < n; k++)
                mag[k] = sqrt(ABSSQR(hc[k]));
        }
    }


    if(phi)
    {
        for(k = 0; k < n; k++)
            phi[k] = atan2(IM(hc[k]), RE(hc[k]));

        if(flag & DSPL_FLAG_UNWRAP)
        {
            res = unwrap(phi, n, M_2PI, 0.8);
            if(res != RES_OK)
                goto exit_label;
        }
    }


    if(tau)
    {
        phi0 = (double*) malloc(n*sizeof(double));
        phi1 = (double*) malloc(n*sizeof(double));
        w0   = (double*) malloc(n*sizeof(double));
        w1   = (double*) malloc(n*sizeof(double));

        w0[0] = w[0] - (w[1] - w[0])*0.02;
        w1[0] = w[0] + (w[1] - w[0])*0.02;

        for(k = 1; k < n; k++)
        {
            w0[k] = w[k] - (w[k] - w[k-1])*0.02;
            w1[k] = w[k] + (w[k] - w[k-1])*0.02;
        }
        res = filter_freq_resp(b, a, ord, w0, n, 
                               DSPL_FLAG_UNWRAP | flag_analog,
                               NULL, phi0, NULL);
        if(res != RES_OK)
            goto exit_label;
        res = filter_freq_resp(b, a, ord, w1, n, 
                               DSPL_FLAG_UNWRAP | flag_analog, 
                               NULL, phi1, NULL);
        if(res != RES_OK)
            goto exit_label;

        for(k = 0; k < n; k++)
            tau[k] = (phi0[k] - phi1[k])/(w1[k] - w0[k]);
    }


exit_label:
    if(hc)
        free(hc);
    if(phi0)
        free(phi0);
    if(phi1)
        free(phi1);
    if(w0)
        free(w0);
    if(w1)
        free(w1);
    return res;
}




#ifdef DOXYGEN_ENGLISH
/*! ****************************************************************************
\ingroup FILTER_ANALYSIS_GROUP
\fn int freqs(double* b, double* a, int ord, double* w, int n, complex_t *h)
\brief Analog filter frequency response \f$ H(j \omega) \f$ calculation

Function calculates analog filter frequency response \f$ H(j \omega)\f$ 
corresponds to transfer function \f$ H(s) \f$:

\f[
  H(s) = \frac  {\sum_{k = 0}^{N} b_k  s^k}
      {\sum_{m = 0}^{N} a_m  s^m},
\f]
here \f$ N \f$ - filter order (equals to `ord`).

\param[in]  b
Pointer to the transfer function \f$ H(s) \f$ 
numerator coefficients vector. \n 
Vector size is `[ord+1 x 1]`. \n \n 

\param[in]  a
Pointer to the transfer function \f$ H(s) \f$ 
denominator coefficients vector. \n 
Vector size is `[ord+1 x 1]`. \n \n 

\param[in]  ord
Filter order. \n
Transfer function \f$ H(s) \f$ numerator and denominator
coefficients number equals `ord+1`. \n \n 

\param[in]  w
Pointer to the angular frequency  \f$ \omega \f$ (rad/s), 
which used for frequency response \f$ H(j \omega) \f$ calculation. \n 
Vector size is `[n x 1]`. \n \n 

\param[in]  n
The size of the angular frequency vector `w`. \n \n 

\param[out]  h
Pointer to the frequency response vector \f$ H(j \omega) \f$, 
corresponds to angular frequency `w`. \n 
Vector size is `[n x 1]`. \n 
Memory must be allocated. \n  \n

\return `RES_OK` if frequency response vector is calculated successfully. \n
Else \ref ERROR_CODE_GROUP "code error".

\author Sergey Bakhurin www.dsplib.org
***************************************************************************** */

#endif
#ifdef DOXYGEN_RUSSIAN
/*! ****************************************************************************
\ingroup FILTER_ANALYSIS_GROUP
\fn int freqs(double* b, double* a, int ord, double* w, int n, complex_t *h)

\brief Расчет комплексного коэффициента передачи 
\f$ H(j \omega) \f$ аналогового фильтра.

Функция рассчитывает значения комплексного коэффициента передачи 
\f$ H(j \omega)\f$ аналогового фильтра, заданного коэффициентами 
передаточной функции \f$ H(s) \f$:

\f[
  H(s) = \frac  {\sum_{k = 0}^{N} b_k  s^k}
      {\sum_{m = 0}^{N} a_m  s^m},
\f]
где \f$ N \f$ - порядок фильтра (параметр `ord`).

Комплексный коэффициент передачи рассчитывается путем 
подстановки \f$ s = j \omega \f$.

\param[in]  b
Указатель на вектор коэффициентов числителя 
передаточной функции \f$ H(s) \f$. \n 
Размер вектора `[ord+1 x 1]`. \n \n 


\param[in]  a
Указатель на вектор коэффициентов знаменателя 
передаточной функции \f$ H(s) \f$. \n 
Размер вектора `[ord+1 x 1]`. \n \n 


\param[in]  ord
Порядок фильтра. Количество коэффициентов числителя и 
знаменателя передаточной функции \f$ H(s) \f$ 
равно `ord+1`. \n \n 


\param[in]  w
Указатель на вектор значений циклической частоты \f$ \omega \f$ (рад/с), 
для которого будет рассчитан комплексный 
коэффициент передачи \f$ H(j \omega) \f$. \n 
Размер вектора `[n x 1]`. \n \n 


\param[in]  n
Размер вектора циклической частоты `w`. \n \n 


\param[out]  h
Указатель на вектор комплексного коэффициента передачи \f$ H(j \omega) \f$, 
рассчитанного для циклической частоты `w`. \n 
Размер вектора `[n x 1]`. \n 
Память должна быть выделена. \n  \n

\return
`RES_OK` Комплексный коэффициент передачи рассчитан успешно. \n
В противном случае \ref ERROR_CODE_GROUP "код ошибки". \n

\author Бахурин Сергей www.dsplib.org
***************************************************************************** */
#endif
int DSPL_API freqs(double* b, double* a, int ord,
                   double* w, int n, complex_t *h)
{
    complex_t jw;
    complex_t *bc = NULL;
    complex_t *ac = NULL;
    complex_t num, den;
    double mag;
    int k;
    int res;

    if(!b || !a || !w || !h)
        return ERROR_PTR;
    if(ord<0)
        return ERROR_FILTER_ORD;
    if(n<1)
        return ERROR_SIZE;

    RE(jw) = 0.0;

    bc = (complex_t*) malloc((ord+1) * sizeof(complex_t));
    res = re2cmplx(b, ord+1, bc);

    if( res!=RES_OK )
        goto exit_label;

    ac = (complex_t*) malloc((ord+1) * sizeof(complex_t));
    res = re2cmplx(a, ord+1, ac);
    if( res!=RES_OK )
        goto exit_label;

    for(k = 0; k < n; k++)
    {
        IM(jw) = w[k];
        res = polyval_cmplx(bc, ord, &jw, 1, &num);
        if(res != RES_OK)
            goto exit_label;
        res = polyval_cmplx(ac, ord, &jw, 1, &den);
        if(res != RES_OK)
            goto exit_label;
        mag = ABSSQR(den);
        if(mag == 0.0)
        {
            res = ERROR_DIV_ZERO;
            goto exit_label;
        }
        mag = 1.0 / mag;
        RE(h[k]) = CMCONJRE(num, den) * mag;
        IM(h[k]) = CMCONJIM(num, den) * mag;
    }
    res = RES_OK;
exit_label:
    if(bc)
        free(bc);
    if(ac)
        free(ac);
    return res;
}






#ifdef DOXYGEN_ENGLISH

#endif
#ifdef DOXYGEN_RUSSIAN

#endif
int DSPL_API freqs_cmplx(double* b, double* a, int ord,
                         complex_t* s, int n, complex_t *h)
{
    complex_t *bc = NULL;
    complex_t *ac = NULL;
    complex_t num, den;
    double mag;
    int k;
    int res;

    if(!b || !a || !s || !h)
        return ERROR_PTR;
    if(ord<0)
        return ERROR_FILTER_ORD;
    if(n<1)
        return ERROR_SIZE;


    bc = (complex_t*) malloc((ord+1) * sizeof(complex_t));
    res = re2cmplx(b, ord+1, bc);

    if( res!=RES_OK )
        goto exit_label;

    ac = (complex_t*) malloc((ord+1) * sizeof(complex_t));
    res = re2cmplx(a, ord+1, ac);
    if( res!=RES_OK )
        goto exit_label;

    for(k = 0; k < n; k++)
    {
        res = polyval_cmplx(bc, ord, s+k, 1, &num);
        if(res != RES_OK)
            goto exit_label;
        res = polyval_cmplx(ac, ord, s+k, 1, &den);
        if(res != RES_OK)
            goto exit_label;
        mag = ABSSQR(den);
        if(mag == 0.0)
        {
            res = ERROR_DIV_ZERO;
            goto exit_label;
        }
        mag = 1.0 / mag;
        RE(h[k]) = CMCONJRE(num, den) * mag;
        IM(h[k]) = CMCONJIM(num, den) * mag;

    }
    res = RES_OK;
    exit_label:
    if(bc)
        free(bc);
    if(ac)
        free(ac);
    return res;
}







#ifdef DOXYGEN_ENGLISH

#endif
#ifdef DOXYGEN_RUSSIAN

#endif
int DSPL_API freqs2time(double* b, double* a, int ord, double fs,
                        int n, fft_t* pfft, double *t, double *h)
{
    double *w = NULL;
    complex_t *hs = NULL;
    complex_t *ht = NULL;
    int err, k;

    if(!b || !a || !t || !h)
        return ERROR_PTR;
    if(ord<1)
        return ERROR_FILTER_ORD;
    if(n<1)
        return ERROR_SIZE;

    w    = (double*)malloc(n*sizeof(double));
    hs = (complex_t*)malloc(n*sizeof(complex_t));


    err = linspace(-fs*0.5, fs*0.5, n, DSPL_PERIODIC, w);
    if(err != RES_OK)
        goto exit_label;

    err = freqs(b, a, ord, w, n, hs);
    if(err != RES_OK)
        goto exit_label;

    err = fft_shift_cmplx(hs, n, hs);
    if(err != RES_OK)
        goto exit_label;

    ht = (complex_t*)malloc(n*sizeof(complex_t));

    err = ifft_cmplx(hs, n, pfft, ht);
    if(err != RES_OK)
    {
        err = idft_cmplx(hs, n, ht);
        if(err != RES_OK)
            goto exit_label;
    }

    for(k = 0; k < n; k++)
    {
        t[k] = (double)k/fs;
        h[k] = RE(ht[k]) * fs;
    }

exit_label:
    if(w)
        free(w);
    if(hs)
        free(hs);
    if(ht)
        free(ht);
    return err;
}



#ifdef DOXYGEN_ENGLISH
/*! ****************************************************************************
\ingroup FILTER_ANALYSIS_GROUP
\fn int freqz(double* b, double* a, int ord, double* w, int n, complex_t *h)

\brief Function calculates the digital filter frequency response 
\f$ H \left(e^{j \omega} \right)\f$ corresponds to transfer function \f$H(z)\f$.

Digital filter transfer function:
\f[
H(z) = \frac{\sum\limits_{k = 0}^{N} b_k  z^{-k}}
    {\sum\limits_{m = 0}^{N} a_m z^{-m}},
\f]
here \f$N\f$ --- filter order (parameter `ord`). \n

Frequency response \f$ H \left(e^{j \omega} \right)\f$ we can get
if substitute \f$z = e^{j \omega} \f$. \n

\param[in]  b
Pointer to the \f$ H(z) \f$ transfer function 
numerator coefficients vector. \n 
Vector size is `[ord+1 x 1]`. \n \n 

\param[in]  a
Pointer to the \f$H(z)\f$ transfer function 
denominator coefficients vector. \n 
Vector size is `[ord+1 x 1]`. \n \n 

\param[in]  ord
Filter order. \n
Transfer function \f$H(z)\f$ numerator 
and denominator coefficients number equals `ord+1`. \n \n 

\param[in]  w
Pointer to the normalized frequency of digital filter 
frequency response \f$ H \left(\mathrm{e}^{j\omega} \right) \f$. \n
Digital filter frequency response is \f$ 2\pi \f$-periodic function,
and vector `w` advisable to set from 0 to \f$ \pi \f$,
or from 0 to \f$ 2\pi \f$, or from \f$ -\pi \f$ to \f$ \pi \f$.
Vector size is `[n x 1]`. \n \n

\param[in]  n
Size of frequency vector `w`. \n  \n 

\param[out]  h
Pointer to the frequency response vector 
\f$ H \left(\mathrm{e}^{j\omega} \right) \f$, 
corresponds to normalized frequency `w`. \n 
Vector size is `[n x 1]`. \n 
Memory must be allocated. \n  \n

\return `RES_OK` if frequaency response vector is calculated successfully. \n
Else \ref ERROR_CODE_GROUP "code error".

\author Sergey Bakhurin www.dsplib.org
***************************************************************************** */
#endif
#ifdef DOXYGEN_RUSSIAN
/*! ****************************************************************************
\ingroup FILTER_ANALYSIS_GROUP
\fn int freqz(double* b, double* a, int ord, double* w, int n, complex_t *h)

\brief  Расчет комплексного коэффициента передачи 
  \f$ H \left(e^{j \omega} \right)\f$ цифрового фильтра.

Функция рассчитывает значения комплексного коэффициента передачи 
\f$ H \left(e^{j \omega} \right)\f$  цифрового фильтра, заданного 
коэффициентами передаточной функции \f$H(z)\f$:

\f[
H(z) = \frac  {\sum_{k = 0}^{N} b_k  z^{-k}}
    {\sum_{m = 0}^{N} a_m z^{-m}},
\f]

где \f$N\f$ --- порядок фильтра (параметр `ord`). \n

Комплексный коэффициент передачи рассчитывается путем 
подстановки \f$z = e^{j \omega} \f$. \n

  
\param[in]  b
Указатель на вектор коэффициентов числителя 
передаточной функции \f$H(z)\f$. \n
Размер вектора `[ord+1 x 1]`. \n \n

\param[in]  a
Указатель на вектор коэффициентов знаменателя 
передаточной функции \f$H(z)\f$. \n
Размер вектора `[ord+1 x 1]`. \n \n

\param[in]  ord
Порядок фильтра. Количество коэффициентов числителя и знаменателя 
передаточной функции \f$H(z)\f$ равно `ord+1`. \n \n

\param[in]  w
Указатель на вектор значений  нормированной циклической частоты \f$\omega\f$, 
для которого будет рассчитан комплексный коэффициент передачи 
\f$ H \left(e^{j \omega} \right)\f$. \n
Размер вектора `[n x 1]`. \n \n
      
\param[in]  n
Размер вектора нормированной циклической частоты `w`. \n \n

\param[out]  h
Указатель на вектор комплексного коэффициента передачи 
\f$ H \left(e^{j \omega} \right)\f$, рассчитанного для 
циклической частоты `w`. \n
Размер вектора `[n x 1]`. \n
Память должна быть выделена. \n \n


\return
`RES_OK` Комплексный коэффициент передачи рассчитан успешно. \n
В противном случае \ref ERROR_CODE_GROUP "код ошибки". \n

\note
Комплексный коэффициент передачи \f$ H \left(e^{j \omega} \right)\f$ 
цифрового фильтра представляет собой \f$ 2 \pi-\f$периодическую функцию 
нормированной циклической частоты \f$\omega\f$.
Поэтому анализ цифровых фильтров целесообразно вести на одном периоде 
повторения \f$ H \left(e^{j \omega} \right)\f$, т.е. в интервале 
\f$\omega\f$ от 0 до \f$2 \pi\f$, или от \f$-\pi\f$ до \f$ \pi\f$.  \n
Кроме того известно, что для фильтра с вещественными коэффициентами 
\f$ H \left(e^{j \omega} \right) = H^* \left(e^{-j \omega} \right)\f$,
а значит, анализ цифрового фильтра с вещественными коэффициентами 
достаточно вести для нормированной частоты \f$\omega\f$ от 0 до \f$\pi\f$.

\author Бахурин Сергей www.dsplib.org
***************************************************************************** */
#endif
int DSPL_API freqz(double* b, double* a, int ord, double* w,
                   int n, complex_t *h)
{
    complex_t jw;
    complex_t *bc = NULL;
    complex_t *ac = NULL;
    complex_t num, den;
    double mag;
    int k;
    int res;

    if(!b || !w || !h)
        return ERROR_PTR;
    if(ord<0)
        return ERROR_FILTER_ORD;
    if(n<1)
        return ERROR_SIZE;


    bc = (complex_t*) malloc((ord+1) * sizeof(complex_t));
    res = re2cmplx(b, ord+1, bc);
    if( res!=RES_OK )
        goto exit_label;

    if(a)
    {
        /* IIR filter if a != NULL */
        ac = (complex_t*) malloc((ord+1) * sizeof(complex_t));
        res = re2cmplx(a, ord+1, ac);
        if( res!=RES_OK )
            goto exit_label;
        for(k = 0; k < n; k++)
        {
            RE(jw) =  cos(w[k]);
            IM(jw) = -sin(w[k]);
            res = polyval_cmplx(bc, ord, &jw, 1, &num);
            if(res != RES_OK)
                goto exit_label;
            res = polyval_cmplx(ac, ord, &jw, 1, &den);
            if(res != RES_OK)
                goto exit_label;
            mag = ABSSQR(den);
            if(mag == 0.0)
            {
                res = ERROR_DIV_ZERO;
                goto exit_label;
            }
            mag = 1.0 / mag;
            RE(h[k]) = CMCONJRE(num, den) * mag;
            IM(h[k]) = CMCONJIM(num, den) * mag;
        }
    }
    else
    {
        /* FIR filter if a == NULL */
        for(k = 0; k < n; k++)
        {
            RE(jw) =    cos(w[k]);
            IM(jw) = -sin(w[k]);
            res = polyval_cmplx(bc, ord, &jw, 1, h+k);
            if(res != RES_OK)
                goto exit_label;
        }
    }
    res = RES_OK;
exit_label:
    if(bc)
        free(bc);
    if(ac)
        free(ac);
    return res;
}








#ifdef DOXYGEN_ENGLISH

#endif
#ifdef DOXYGEN_RUSSIAN

#endif
int DSPL_API unwrap(double* phi, int n, double lev, double mar)
{
    double a[2] = {0.0, 0.0};
    double d;
    double th;
    int k;
    int flag = 1;

    if(!phi)
        return ERROR_PTR;

    if(n<1)
        return ERROR_SIZE;

    if(lev<=0 || mar <=0)
        return ERROR_UNWRAP;

    th = mar*lev;
    while(flag)
    {
        flag = 0;
        a[0] = a[1] = 0.0;
        for(k = 0; k<n-1; k++)
        {
            d = phi[k+1] - phi[k];
            if( d > th)
            {
                a[0] -= lev;
                flag = 1;
            }
            if( d < -th)
            {
                a[0] += lev;
                flag = 1;
            }
            phi[k]+=a[1];
            a[1] = a[0];
        }
        phi[n-1]+=a[1];
    }

    return RES_OK;
}

