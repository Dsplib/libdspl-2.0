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
#include "dspl_internal.h"


#ifdef DOXYGEN_ENGLISH
/*! ****************************************************************************
\ingroup DFT_GROUP
\fn int ifft_cmplx(complex_t* x, int n, fft_t* pfft, complex_t* y)
\brief Inverse fast Fourier transform

Function calculates \f$ n \f$-point IFFT of complex data 
\f$ x(m) \f$, \f$ m = 0 \ldots n-1 \f$. \n
\f[
  Y(k) = \frac{1}{N} \sum_{m = 0}^{n-1} x(m)   \exp 
   \left( j   \frac{2\pi}{n}    m    k \right),
\f]
here \f$ k = 0 \ldots n-1 \f$.


\param[in]  x
Pointer to the input vector \f$x(m)\f$, 
\f$ m = 0 \ldots n-1 \f$.  \n
Vector size is `[n x 1]`.  \n \n

\param[in]  n
IFFT size \f$n\f$. \n
IFFT size can be composite: 
\f$n = n_0 \times n_1 \times n_2 \times \ldots \times n_p \times m\f$,
here \f$n_i = 2,3,5,7\f$, а \f$m \f$ -- 
simple number less than 46340
(see \ref fft_create function). \n \n

\param[in]  pfft
Pointer to the `fft_t` object.  \n
This pointer cannot be `NULL`.  \n
Structure \ref fft_t should be previously once
filled with the \ref fft_create function, and the memory should be
cleared before exiting by the \ref fft_free function. \n \n

\param[out] y
Pointer to the IFFT result vector \f$Y(k)\f$, 
\f$ k = 0 \ldots n-1 \f$. \n 
Vector size is `[n x 1]`.  \n
Memory must be allocated. \n \n

\return
`RES_OK` if IFFT is calculated successfully. \n
Else \ref ERROR_CODE_GROUP "code error".

IFFT example:

\include ifft_cmplx_test.c

Result:

\verbatim
| x[ 0] =  1.000   0.000 | y[ 0] = -0.517   0.686 | z[ 0] =   1.000   0.000 |
| x[ 1] =  0.540   0.841 | y[ 1] = -0.943   0.879 | z[ 1] =   0.540   0.841 |
| x[ 2] = -0.416   0.909 | y[ 2] = -2.299   1.492 | z[ 2] =  -0.416   0.909 |
| x[ 3] = -0.990   0.141 | y[ 3] = 16.078  -6.820 | z[ 3] =  -0.990   0.141 |
| x[ 4] = -0.654  -0.757 | y[ 4] =  2.040  -0.470 | z[ 4] =  -0.654  -0.757 |
| x[ 5] =  0.284  -0.959 | y[ 5] =  1.130  -0.059 | z[ 5] =   0.284  -0.959 |
| x[ 6] =  0.960  -0.279 | y[ 6] =  0.786   0.097 | z[ 6] =   0.960  -0.279 |
| x[ 7] =  0.754   0.657 | y[ 7] =  0.596   0.183 | z[ 7] =   0.754   0.657 |
| x[ 8] = -0.146   0.989 | y[ 8] =  0.470   0.240 | z[ 8] =  -0.146   0.989 |
| x[ 9] = -0.911   0.412 | y[ 9] =  0.375   0.283 | z[ 9] =  -0.911   0.412 |
| x[10] = -0.839  -0.544 | y[10] =  0.297   0.318 | z[10] =  -0.839  -0.544 |
| x[11] =  0.004  -1.000 | y[11] =  0.227   0.350 | z[11] =   0.004  -1.000 |
| x[12] =  0.844  -0.537 | y[12] =  0.161   0.380 | z[12] =   0.844  -0.537 |
| x[13] =  0.907   0.420 | y[13] =  0.094   0.410 | z[13] =   0.907   0.420 |
| x[14] =  0.137   0.991 | y[14] =  0.023   0.442 | z[14] =   0.137   0.991 |
| x[15] = -0.760   0.650 | y[15] = -0.059   0.479 | z[15] =  -0.760   0.650 |
| x[16] = -0.958  -0.288 | y[16] = -0.161   0.525 | z[16] =  -0.958  -0.288 |
| x[17] = -0.275  -0.961 | y[17] = -0.300   0.588 | z[17] =  -0.275  -0.961 |
\endverbatim

\author Sergey Bakhurin www.dsplib.org
***************************************************************************** */
#endif
#ifdef DOXYGEN_RUSSIAN
/*! ****************************************************************************
\ingroup DFT_GROUP
\fn int ifft_cmplx(complex_t* x, int n, fft_t* pfft, complex_t* y)
\brief Обратное быстрое преобразование Фурье

Функция рассчитывает \f$ n \f$-точечное обратное быстрое преобразование Фурье 
от \f$ x(m) \f$, \f$ m = 0 \ldots n-1 \f$. \n
\f[
  Y(k) = \frac{1}{N} \sum_{m = 0}^{n-1} x(m)   \exp 
   \left( j   \frac{2\pi}{n}    m    k \right),
\f]
где \f$ k = 0 \ldots n-1 \f$.

Для расчета используется алгоритм БПФ составной длины.

\param[in]  x
Указатель на входной комплексный вектор \f$x(m)\f$, 
\f$ m = 0 \ldots n-1 \f$.  \n
Размер вектора `[n x 1]`.  \n \n

\param[in]  n
Размер ОБПФ \f$n\f$. \n
Размер ОБПФ может быть составным вида 
\f$n = n_0 \times n_1 \times n_2 \times \ldots \times n_p \times m\f$,
где \f$n_i = 2,3,5,7\f$, а \f$m \f$ -- 
произвольный простой множитель не превосходящий 46340
(см. описание функции \ref fft_create). \n \n

\param[in]  pfft
Указатель на структуру `fft_t`.  \n
Указатель не должен быть `NULL`.  \n
Структура \ref fft_t должна быть предварительно однократно
заполнена функцией  \ref fft_create, и память должна быть 
очищена перед выходом функцией \ref fft_free. \n \n

\param[out] y
Указатель на вектор результата ОБПФ \f$Y(k)\f$, 
\f$ k = 0 \ldots n-1 \f$. Размер вектора `[n x 1]`.  \n
Память должна быть выделена. \n \n

\return
`RES_OK` если расчет произведен успешно.  \n
В противном случае \ref ERROR_CODE_GROUP "код ошибки". \n \n

Пример использования функции `fft`:

\include ifft_cmplx_test.c

Результат работы программы:

\verbatim
| x[ 0] =  1.000   0.000 | y[ 0] = -0.517   0.686 | z[ 0] =   1.000   0.000 |
| x[ 1] =  0.540   0.841 | y[ 1] = -0.943   0.879 | z[ 1] =   0.540   0.841 |
| x[ 2] = -0.416   0.909 | y[ 2] = -2.299   1.492 | z[ 2] =  -0.416   0.909 |
| x[ 3] = -0.990   0.141 | y[ 3] = 16.078  -6.820 | z[ 3] =  -0.990   0.141 |
| x[ 4] = -0.654  -0.757 | y[ 4] =  2.040  -0.470 | z[ 4] =  -0.654  -0.757 |
| x[ 5] =  0.284  -0.959 | y[ 5] =  1.130  -0.059 | z[ 5] =   0.284  -0.959 |
| x[ 6] =  0.960  -0.279 | y[ 6] =  0.786   0.097 | z[ 6] =   0.960  -0.279 |
| x[ 7] =  0.754   0.657 | y[ 7] =  0.596   0.183 | z[ 7] =   0.754   0.657 |
| x[ 8] = -0.146   0.989 | y[ 8] =  0.470   0.240 | z[ 8] =  -0.146   0.989 |
| x[ 9] = -0.911   0.412 | y[ 9] =  0.375   0.283 | z[ 9] =  -0.911   0.412 |
| x[10] = -0.839  -0.544 | y[10] =  0.297   0.318 | z[10] =  -0.839  -0.544 |
| x[11] =  0.004  -1.000 | y[11] =  0.227   0.350 | z[11] =   0.004  -1.000 |
| x[12] =  0.844  -0.537 | y[12] =  0.161   0.380 | z[12] =   0.844  -0.537 |
| x[13] =  0.907   0.420 | y[13] =  0.094   0.410 | z[13] =   0.907   0.420 |
| x[14] =  0.137   0.991 | y[14] =  0.023   0.442 | z[14] =   0.137   0.991 |
| x[15] = -0.760   0.650 | y[15] = -0.059   0.479 | z[15] =  -0.760   0.650 |
| x[16] = -0.958  -0.288 | y[16] = -0.161   0.525 | z[16] =  -0.958  -0.288 |
| x[17] = -0.275  -0.961 | y[17] = -0.300   0.588 | z[17] =  -0.275  -0.961 |
\endverbatim

\author Бахурин Сергей www.dsplib.org 
***************************************************************************** */
#endif
int DSPL_API ifft_cmplx(complex_t *x, int n, fft_t* pfft, complex_t* y)
{
    int err, k;
    double norm;

    if(!x || !pfft || !y)
        return ERROR_PTR;
    if(n<1)
        return ERROR_SIZE;


    err = fft_create(pfft, n);
    if(err != RES_OK)
        return err;

    memcpy(pfft->t1, x, n*sizeof(complex_t));
    for(k = 0; k < n; k++)
        IM(pfft->t1[k]) = -IM(pfft->t1[k]);

    err = fft_krn(pfft->t1, pfft->t0, pfft, n, 0);

    if(err!=RES_OK)
        return err;

    norm = 1.0 / (double)n;
    for(k = 0; k < n; k++)
    {
         RE(y[k]) =  RE(pfft->t0[k])*norm;
         IM(y[k]) = -IM(pfft->t0[k])*norm;
    }
    return RES_OK;
}


#ifdef DOXYGEN_ENGLISH
/*! ****************************************************************************
\ingroup DFT_GROUP
\fn int fft(double* x, int n, fft_t* pfft, complex_t* y)
\brief Fast Fourier transform for the real vector.

Function calculated \f$ n \f$-points FFT for the real vector 
\f$ x(m) \f$, \f$ m = 0 \ldots n-1 \f$. \n
\f[
  Y(k) = \sum_{m = 0}^{n-1} x(m) \exp 
  \left( -j   \frac{2\pi}{n} m k \right),
\f]
here \f$ k = 0 \ldots n-1 \f$.


\param[in]  x
Pointer to the input real vector \f$x(m)\f$, 
\f$ m = 0 \ldots n-1 \f$.  \n
Vector size is `[n x 1]`.  \n \n

\param[in]  n
FFT size \f$n\f$. \n
FFT size can be composite: 
\f$n = n_0 \times n_1 \times n_2 \times \ldots \times n_p \times m\f$,
here \f$n_i = 2,3,5,7\f$, а \f$m \f$ -- 
simple number less than 46340
(see \ref fft_create function). \n \n

\param[in]  pfft
Pointer to the `fft_t` object.  \n
This pointer cannot be `NULL`.  \n
Structure \ref fft_t should be previously once
filled with the \ref fft_create function, and the memory should be
cleared before exiting by the \ref fft_free function. \n \n

\param[out] y
Pointer to the FFT result complex vector \f$Y(k)\f$, 
\f$ k = 0 \ldots n-1 \f$. \n
Vector size is `[n x 1]`. \n
Memory must be allocated. \n \n

\return
`RES_OK` if FFT is calculated successfully. \n
Else \ref ERROR_CODE_GROUP "code error".

Example:

\include fft_test.c

Result:

\verbatim
y[ 0] =    91.000    0.000
y[ 1] =    -7.000   30.669
y[ 2] =    -7.000   14.536
y[ 3] =    -7.000    8.778
y[ 4] =    -7.000    5.582
y[ 5] =    -7.000    3.371
y[ 6] =    -7.000    1.598
y[ 7] =    -7.000    0.000
y[ 8] =    -7.000   -1.598
y[ 9] =    -7.000   -3.371
y[10] =    -7.000   -5.582
y[11] =    -7.000   -8.778
y[12] =    -7.000  -14.536
y[13] =    -7.000  -30.669
\endverbatim

\author Sergey Bakhurin www.dsplib.org 
***************************************************************************** */
#endif
#ifdef DOXYGEN_RUSSIAN
/*! ****************************************************************************
\ingroup DFT_GROUP
\fn int fft(double* x, int n, fft_t* pfft, complex_t* y)
\brief Быстрое преобразование Фурье вещественного сигнала

Функция рассчитывает \f$ n \f$-точечное  быстрое преобразование Фурье 
вещественного сигнала \f$ x(m) \f$, \f$ m = 0 \ldots n-1 \f$. \n
\f[
  Y(k) = \sum_{m = 0}^{n-1} x(m) \exp 
  \left( -j \frac{2\pi}{n} m k \right),
\f]
где \f$ k = 0 \ldots n-1 \f$.

Для расчета используется алгоритм БПФ составной длины.

\param[in]  x
Указатель на вектор вещественного входного сигнала \f$x(m)\f$, 
\f$ m = 0 \ldots n-1 \f$.  \n
Размер вектора `[n x 1]`.  \n \n

\param[in]  n
Размер БПФ \f$n\f$. \n
Размер БПФ может быть составным вида 
\f$n = n_0 \times n_1 \times n_2 \times \ldots \times n_p \times m\f$,
где \f$n_i = 2,3,5,7\f$, а \f$m \f$ -- 
произвольный простой множитель не превосходящий 46340
(см. описание функции \ref fft_create). \n \n

\param[in]  pfft
Указатель на структуру `fft_t`. \n
Указатель не должен быть `NULL`. \n
Структура \ref fft_t должна быть предварительно однократно
заполнена функцией  \ref fft_create, и память должна быть 
очищена перед выходом функцией \ref fft_free. \n \n

\param[out] y
Указатель на комплексный вектор результата БПФ \f$Y(k)\f$, 
\f$ k = 0 \ldots n-1 \f$. \n
Размер вектора `[n x 1]`. \n
Память должна быть выделена. \n \n

\return
`RES_OK` если расчет произведен успешно.  \n
 В противном случае \ref ERROR_CODE_GROUP "код ошибки". \n \n

Пример использования функции `fft`:

\include fft_test.c

Результат работы программы:

\verbatim
y[ 0] =    91.000    0.000
y[ 1] =    -7.000   30.669
y[ 2] =    -7.000   14.536
y[ 3] =    -7.000    8.778
y[ 4] =    -7.000    5.582
y[ 5] =    -7.000    3.371
y[ 6] =    -7.000    1.598
y[ 7] =    -7.000    0.000
y[ 8] =    -7.000   -1.598
y[ 9] =    -7.000   -3.371
y[10] =    -7.000   -5.582
y[11] =    -7.000   -8.778
y[12] =    -7.000  -14.536
y[13] =    -7.000  -30.669
\endverbatim

\author Бахурин Сергей www.dsplib.org 
***************************************************************************** */
#endif
int DSPL_API fft(double* x, int n, fft_t* pfft, complex_t* y)
{
    int err;

    if(!x || !pfft || !y)
        return ERROR_PTR;
    if(n<1)
        return ERROR_SIZE;


    err = fft_create(pfft, n);
    if(err != RES_OK)
        return err;

    re2cmplx(x, n, pfft->t1);

    return fft_krn(pfft->t1, y, pfft, n, 0);
}



#ifdef DOXYGEN_ENGLISH

#endif
#ifdef DOXYGEN_RUSSIAN

#endif
int DSPL_API fft_abs(double* x, int n, fft_t* pfft, 
                     double fs, int flag,
                     double* mag, double* freq)
{
    int k, err = RES_OK;
    complex_t *X = NULL;
    if(!x || !pfft)
        return ERROR_PTR;

    if(n<1)
        return ERROR_SIZE;

    if(mag)
    {    
        X = (complex_t*)malloc(n*sizeof(complex_t));
        err = fft(x, n, pfft, X);
        if(err!=RES_OK)
            goto error_proc;
        

        for(k = 0; k < n; k++)
            mag[k] = ABS(X[k]);
        if(flag & DSPL_FLAG_FFT_SHIFT)
        {
            err = fft_shift(mag, n, mag);
            if(err!=RES_OK)
                goto error_proc;
        }
    }
    
    if(freq)
    {
        if(flag & DSPL_FLAG_FFT_SHIFT)
            if(n%2)
                err = linspace(-fs*0.5 + fs*0.5/(double)n, 
                                fs*0.5 - fs*0.5/(double)n, 
                                n, DSPL_SYMMETRIC, freq);
            else
                err = linspace(-fs*0.5, fs*0.5, n, DSPL_PERIODIC, freq);
        else
            err = linspace(0, fs, n, DSPL_PERIODIC, freq);
    } 

error_proc:
    if(X)
        free(X);

    return err;
}




#ifdef DOXYGEN_ENGLISH

#endif
#ifdef DOXYGEN_RUSSIAN

#endif
int DSPL_API fft_abs_cmplx(complex_t* x, int n, fft_t* pfft, 
                           double fs, int flag,
                           double* mag, double* freq)
{
    int k, err = RES_OK;
    complex_t *X = NULL;
    
    if(!x || !pfft)
        return ERROR_PTR;
    
    if(n<1)
        return ERROR_SIZE;
    
    if(mag)
    {    
        X = (complex_t*)malloc(n*sizeof(complex_t));
        err = fft_cmplx(x, n, pfft, X);
        if(err!=RES_OK)
            goto error_proc;
        
        
        for(k = 0; k < n; k++)
            mag[k] = ABS(X[k]);
          
        if(flag & DSPL_FLAG_FFT_SHIFT)
        {
            err = fft_shift(mag, n, mag);
            if(err!=RES_OK)
                goto error_proc;
        }
    }
    
    if(freq)
    {
        if(flag & DSPL_FLAG_FFT_SHIFT)
            if(n%2)
                err = linspace(-fs*0.5 + fs*0.5/(double)n, 
                                fs*0.5 - fs*0.5/(double)n, 
                                n, DSPL_SYMMETRIC, freq);
            else
                err = linspace(-fs*0.5, fs*0.5, n, DSPL_PERIODIC, freq);
        else
            err = linspace(0, fs, n, DSPL_PERIODIC, freq);
    }         
error_proc:    
    if(X)
        free(X);
        
    return err;
}



#ifdef DOXYGEN_ENGLISH
/*! ****************************************************************************
\ingroup DFT_GROUP
\fn int fft_cmplx(complex_t* x, int n, fft_t* pfft, complex_t* y)
\brief Fast Fourier transform for the complex vector.

Function calculated \f$ n \f$-points FFT for the complex vector 
\f$ x(m) \f$, \f$ m = 0 \ldots n-1 \f$. \n
\f[
  Y(k) = \sum_{m = 0}^{n-1} x(m) \exp \left( -j \frac{2\pi}{n} m k \right),
\f]
here \f$ k = 0 \ldots n-1 \f$.

\param[in]  x
Pointer to the input complex vector \f$x(m)\f$, 
\f$ m = 0 \ldots n-1 \f$.  \n
Vector size is `[n x 1]`.  \n \n

\param[in]  n
FFT size \f$n\f$. \n
FFT size can be composite: 
\f$n = n_0 \times n_1 \times n_2 \times \ldots \times n_p \times m\f$,
here \f$n_i = 2,3,5,7\f$, а \f$m \f$ -- 
simple number less than 46340
(see \ref fft_create function). \n \n

\param[in]  pfft
Pointer to the `fft_t` object.  \n
This pointer cannot be `NULL`.  \n
Structure \ref fft_t should be previously once
filled with the \ref fft_create function, and the memory should be
cleared before exiting by the \ref fft_free function. \n \n

\param[out] y
Pointer to the FFT result complex vector \f$Y(k)\f$, 
\f$ k = 0 \ldots n-1 \f$. \n
Vector size is `[n x 1]`. \n
Memory must be allocated. \n \n

\return
`RES_OK` if FFT is calculated successfully. \n
Else \ref ERROR_CODE_GROUP "code error".

Example:

\include fft_cmplx_test.c

Result:

\verbatim
y[ 0] =    -0.517    0.686
y[ 1] =    -0.943    0.879
y[ 2] =    -2.299    1.492
y[ 3] =    16.078   -6.820
y[ 4] =     2.040   -0.470
y[ 5] =     1.130   -0.059
y[ 6] =     0.786    0.097
y[ 7] =     0.596    0.183
y[ 8] =     0.470    0.240
y[ 9] =     0.375    0.283
y[10] =     0.297    0.318
y[11] =     0.227    0.350
y[12] =     0.161    0.380
y[13] =     0.094    0.410
y[14] =     0.023    0.442
y[15] =    -0.059    0.479
y[16] =    -0.161    0.525
y[17] =    -0.300    0.588
\endverbatim

\author Sergey Bakhurin www.dsplib.org 
***************************************************************************** */
#endif
#ifdef DOXYGEN_RUSSIAN
/*! ****************************************************************************
\ingroup DFT_GROUP
\fn int fft_cmplx(complex_t* x, int n, fft_t* pfft, complex_t* y)
\brief Быстрое преобразование Фурье комплексного сигнала

Функция рассчитывает \f$ n \f$-точечное  быстрое преобразование Фурье 
комплексного сигнала \f$ x(m) \f$, \f$ m = 0 \ldots n-1 \f$. \n
\f[
  Y(k) = \sum_{m = 0}^{n-1} x(m) \exp \left( -j \frac{2\pi}{n} m k \right),
\f]
где \f$ k = 0 \ldots n-1 \f$.

Для расчета используется алгоритм БПФ составной длины.

\param[in]  x
Указатель на вектор комплексного 
входного сигнала \f$x(m)\f$, \f$ m = 0 \ldots n-1 \f$.  \n
Размер вектора `[n x 1]`. \n \n

\param[in]  n
Размер БПФ \f$n\f$. \n
Размер БПФ может быть составным вида 
\f$ n = n_0 \times n_1 \times n_2 \times n_3 \times \ldots 
\times n_p \times m \f$,
где \f$n_i = 2,3,5,7\f$, а \f$m \f$ -- 
произвольный простой множитель не превосходящий 46340
(см. описание функции \ref fft_create). \n \n

\param[in]  pfft
Указатель на структуру `fft_t`. \n
Указатель не должен быть `NULL`.  \n
Структура \ref fft_t должна быть предварительно однократно
заполнена функцией  \ref fft_create, и память должна быть 
очищена перед выходом функцией \ref fft_free. \n \n

\param[out] y
Указатель на комплексный вектор 
результата БПФ \f$Y(k)\f$, 
\f$ k = 0 \ldots n-1 \f$. 
Размер вектора `[n x 1]`.  \n
Память должна быть выделена. \n \n

\return
`RES_OK` если расчет произведен успешно.  \n
 В противном случае \ref ERROR_CODE_GROUP "код ошибки". \n \n

Пример использования функции `fft`:

\include fft_cmplx_test.c

Результат работы программы:

\verbatim
y[ 0] =    -0.517    0.686
y[ 1] =    -0.943    0.879
y[ 2] =    -2.299    1.492
y[ 3] =    16.078   -6.820
y[ 4] =     2.040   -0.470
y[ 5] =     1.130   -0.059
y[ 6] =     0.786    0.097
y[ 7] =     0.596    0.183
y[ 8] =     0.470    0.240
y[ 9] =     0.375    0.283
y[10] =     0.297    0.318
y[11] =     0.227    0.350
y[12] =     0.161    0.380
y[13] =     0.094    0.410
y[14] =     0.023    0.442
y[15] =    -0.059    0.479
y[16] =    -0.161    0.525
y[17] =    -0.300    0.588
\endverbatim

\author Бахурин Сергей www.dsplib.org 
***************************************************************************** */
#endif
int DSPL_API fft_cmplx(complex_t* x, int n, fft_t* pfft, complex_t* y)
{
    int err;

    if(!x || !pfft || !y)
        return ERROR_PTR;
    if(n<1)
        return ERROR_SIZE;

    err = fft_create(pfft, n);
    if(err != RES_OK)
        return err;

    memcpy(pfft->t1, x, n*sizeof(complex_t));

    return fft_krn(pfft->t1, y, pfft, n, 0);
}



#ifdef DOXYGEN_ENGLISH

#endif
#ifdef DOXYGEN_RUSSIAN

#endif
int fft_krn(complex_t* t0, complex_t* t1, fft_t* p, int n, int addr)
{
    int n1, n2, k, m, i;
    complex_t *pw = p->w+addr;
    complex_t tmp;
    
    n1 = 1;
    if(n % 4096 == 0) { n1 = 4096; goto label_size; }
    if(n % 2048 == 0) { n1 = 2048; goto label_size; }
    if(n % 1024 == 0) { n1 = 1024; goto label_size; }
    if(n %  512 == 0) { n1 =  512; goto label_size; }
    if(n %  256 == 0) { n1 =  256; goto label_size; }
    if(n %  128 == 0) { n1 =  128; goto label_size; }
    if(n %   64 == 0) { n1 =   64; goto label_size; }
    if(n %   32 == 0) { n1 =   32; goto label_size; }
    if(n %   16 == 0) { n1 =   16; goto label_size; }
    if(n %    7 == 0) { n1 =    7; goto label_size; }
    if(n %    8 == 0) { n1 =    8; goto label_size; }
    if(n %    5 == 0) { n1 =    5; goto label_size; }
    if(n %    4 == 0) { n1 =    4; goto label_size; }
    if(n %    3 == 0) { n1 =    3; goto label_size; }
    if(n %    2 == 0) { n1 =    2; goto label_size; }

label_size:
    if(n1 == 1)
    {
        for(k = 0; k < n; k++)
        {
            RE(t1[k]) = IM(t1[k]) = 0.0;
            for(m = 0; m < n; m++)
            {
                i = (k*m) % n;
                RE(tmp) = CMRE(t0[m], pw[i]);
                IM(tmp) = CMIM(t0[m], pw[i]);
                RE(t1[k]) += RE(tmp);
                IM(t1[k]) += IM(tmp);
            }
        }
    }
    else
    {
        n2 = n / n1;
        
        if(n2>1)
        {
            memcpy(t1, t0, n*sizeof(complex_t));
            matrix_transpose_cmplx(t1, n2, n1, t0);
        }
        
        if(n1 == 4096)
            for(k = 0; k < n2; k++)
                dft4096(t0+4096*k, t1+4096*k, p->w4096, p->w256);
              
        if(n1 == 2048)
            for(k = 0; k < n2; k++)
                dft2048(t0+2048*k, t1+2048*k, p->w2048, p->w32, p->w64);
        
        if(n1 == 1024)
            for(k = 0; k < n2; k++)
                dft1024(t0+1024*k, t1+1024*k, p->w1024, p->w32);
              
        if(n1 == 512)
            for(k = 0; k < n2; k++)
                dft512(t0+512*k, t1+512*k, p->w512, p->w32);
              
        if(n1 == 256)
            for(k = 0; k < n2; k++)
                dft256(t0+256*k, t1+256*k, p->w256);
              
        if(n1 == 128)
            for(k = 0; k < n2; k++)
                dft128(t0+128*k, t1+128*k, p->w128);
              
        if(n1 == 64)
            for(k = 0; k < n2; k++)
                dft64(t0+64*k, t1+64*k, p->w64);

        if(n1 == 32)
            for(k = 0; k < n2; k++)
                dft32(t0+32*k, t1+32*k, p->w32);
        
        if(n1 == 16)
            for(k = 0; k < n2; k++)
                dft16(t0+16*k, t1+16*k);
                
        if(n1 == 7)
            for(k = 0; k < n2; k++)
                dft7(t0+7*k, t1+7*k);
            
        if(n1 == 8)
            for(k = 0; k < n2; k++)
                dft8(t0+8*k, t1+8*k); 
                
        if(n1 == 5)
            for(k = 0; k < n2; k++)
                dft5(t0+5*k, t1+5*k);
     
        if(n1 == 4)
            for(k = 0; k < n2; k++)
                dft4(t0+4*k, t1+4*k);
        
        if(n1 == 3)
            for(k = 0; k < n2; k++)
                dft3(t0+3*k, t1+3*k);
        
        if(n1 == 2)
            for(k = 0; k < n2; k++)
                dft2(t0+2*k, t1+2*k);

        if(n2 > 1)
        {

            for(k =0; k < n; k++)
            {
                RE(t0[k]) = CMRE(t1[k], pw[k]);
                IM(t0[k]) = CMIM(t1[k], pw[k]);
            }

            matrix_transpose_cmplx(t0, n1, n2, t1);
            
            for(k = 0; k < n1; k++)
            {
                fft_krn(t1+k*n2, t0+k*n2, p, n2, addr+n);
            }
            
            matrix_transpose_cmplx(t0, n2, n1, t1);
        }
    }
    return RES_OK;
}




#ifdef DOXYGEN_ENGLISH
/*! ****************************************************************************
\ingroup DFT_GROUP
\fn int fft_create(fft_t* pfft, int n)
\brief Function creates and fill `fft_t` structure.

The function allocates memory and calculates twiddle factors 
 of the `n`-point FFT for the structure` fft_t`.

\param[in,out]  pfft
Pointer to the `fft_t` object.  \n
Pointer cannot be `NULL`.  \n \n

\param[in]  n
FFT size \f$n\f$. \n
FFT size can be composite 
\f$n = n_0 \times n_1 \times n_2 \ldots \times n_p \times m\f$,
here \f$n_i = 2,3,5,7\f$, and \f$m \f$ -- 
 arbitrary prime factor not exceeding 46340. \n
Thus, the FFT algorithm supports arbitrary integer lengths.
degrees of numbers 2,3,5,7, as well as their various combinations.  \n
For example, with \f$ n = 725760 \f$ the structure will be successfully filled, 
because 
\f$ 725760 = 2 \cdot 3 \cdot 4 \cdot 5 \cdot 6 \cdot 7 \cdot 9 \cdot 16 \f$. \n
If \f$ n = 172804 = 43201 \cdot 4 \f$ then the structure will also be 
successfully filled, because the simple factor in \f$ n \f$ does not 
exceed 46340. \n
For size \f$ n = 13 \cdot 17 \cdot 23 \cdot 13 = 66079 \f$
the function will return an error since 66079 is greater than 46340 and is 
not the result of the product of numbers 2,3,5,7. \n \n

\return
`RES_OK` if FFT structure is created and filled successfully. \n
Else \ref ERROR_CODE_GROUP "code error".

\note
Some compilers do not nullify its contents when creating a structure.
Therefore, it is recommended to reset the structure after its declaration:
\code{.cpp}
fft_t pfft = {0};     // fill and fields of fft_t as zeros
int n = 64;           // FFT size

int err;

// Create fft_t object for 64-points FFT 

err = fft_create(&pfft, n);

// ................................... 

// Clear fft_t structure

fft_free(&pfft);
\endcode

Before exiting the program, the memory allocated in the structure
need to clear by  \ref fft_free function. \n \n

\note
The "magic number" 46340 because \f$\sqrt{2^{31}} = 46340.95\f$. \n

\author Sergey Bakhurin www.dsplib.org 
***************************************************************************** */
#endif
#ifdef DOXYGEN_RUSSIAN
/*! ****************************************************************************
\ingroup DFT_GROUP
\fn int fft_create(fft_t* pfft, int n)
\brief Заполнение структуры `fft_t` для алгоритма БПФ

Функция производит выделение памяти и рассчет векторов 
поворотных коэффициентов `n`-точечного БПФ для структуры `fft_t`.

\param[in,out]  pfft
Указатель на структуру `fft_t`.  \n
Указатель не должен быть `NULL`.  \n \n

\param[in]  n
Размер БПФ \f$n\f$. \n
Размер БПФ может быть составным вида 
\f$n = n_0 \times n_1 \times n_2 \ldots \times n_p \times m\f$,
где \f$n_i = 2,3,5,7\f$, а \f$m \f$ -- 
произвольный простой множитель не превосходящий 46340. \n
Таким образом алгоритм БПФ поддерживает произвольные длины, равные целой 
степени чисел 2,3,5,7, а также различные их комбинации.  \n
Так например, при \f$ n = 725760 \f$ структура будет успешно заполнена, 
потому что 
\f$725760 = 2 \cdot 3 \cdot 4 \cdot 5 \cdot 6 \cdot 7 \cdot 9 \cdot 16 \f$, 
т.е. получается как произведение множителей 2,3,5,7. \n
При \f$ n = 172804 = 43201 \cdot 4 \f$ структура также будет успешно заполнена, 
потому что простой множитель входящий в \f$n\f$ не превосходит 46340. \n
Для размера \f$ n = 13 \cdot 17 \cdot 23 \cdot 13 = 66079 \f$ 
функция вернет ошибку, поскольку 66079 больше 46340 и не является результатом 
произведения чисел 2,3,5,7. \n \n

\return
`RES_OK` если структура заполнена успешно. \n
В противном случае \ref ERROR_CODE_GROUP "код ошибки". \n \n

\note
Некоторые компиляторы при создании структуры не обнуляют ее содержимое. 
Поэтому рекомендуется произвести обнуление структуры после ее объявления:
\code{.cpp}
fft_t pfft = {0};     // объявляем объект fft_t
int n = 64;           // Размер БПФ

int err;

// создаем объект для 64-точечного БПФ 

err = fft_create(&pfft, n);

// ................................... 

// очистить память объекта БПФ

fft_free(&pfft);
\endcode

Перед выходом из программы выделенную в структуре память 
необходимо очистить функцией \ref fft_free . \n \n

\note
Магия числа 46340 заключается в том, что \f$\sqrt{2^{31}} = 46340.95\f$. \n

\author Бахурин Сергей www.dsplib.org 
***************************************************************************** */
#endif
int DSPL_API fft_create(fft_t* pfft, int n)
{

    int n1, n2, addr, s, k, m, nw, err;
    double phi;
    s = n;
    nw = addr = 0;

    if(pfft->n == n)
        return RES_OK;

    while(s > 1)
    {
        n2 = 1;
        if(s%4096 == 0)  { n2 = 4096; goto label_size; }
        if(s%2048 == 0)  { n2 = 2048; goto label_size; }
        if(s%1024 == 0)  { n2 = 1024; goto label_size; }
        if(s%512  == 0)  { n2 =  512; goto label_size; }
        if(s%256  == 0)  { n2 =  256; goto label_size; }
        if(s%128  == 0)  { n2 =  128; goto label_size; }
        if(s% 64  == 0)  { n2 =   64; goto label_size; }
        if(s% 32  == 0)  { n2 =   32; goto label_size; }
        if(s% 16  == 0)  { n2 =   16; goto label_size; }
        if(s%  7  == 0)  { n2 =    7; goto label_size; }
        if(s%  8  == 0)  { n2 =    8; goto label_size; }
        if(s%  5  == 0)  { n2 =    5; goto label_size; }
        if(s%  4  == 0)  { n2 =    4; goto label_size; }
        if(s%  3  == 0)  { n2 =    3; goto label_size; }
        if(s%  2  == 0)  { n2 =    2; goto label_size; }


label_size:
        if(n2 == 1)
        {
            if(s > FFT_COMPOSITE_MAX)
            {
                err = ERROR_FFT_SIZE;
                goto error_proc;
            }
            
            nw += s;
            pfft->w = pfft->w ? 
                      (complex_t*) realloc(pfft->w,  nw*sizeof(complex_t)):
                      (complex_t*) malloc(           nw*sizeof(complex_t));
            for(k = 0; k < s; k++)
            {
                phi = - M_2PI * (double)k / (double)s;
                RE(pfft->w[addr]) = cos(phi);
                IM(pfft->w[addr]) = sin(phi);
                addr++;
            }
            s = 1;
        }
        else
        {
            n1 = s / n2;
            nw += s;
            pfft->w = pfft->w ? 
                      (complex_t*) realloc(pfft->w,    nw*sizeof(complex_t)):
                      (complex_t*) malloc(             nw*sizeof(complex_t));

            for(k = 0; k < n1; k++)
            {
                for(m = 0; m < n2; m++)
                {
                    phi = - M_2PI * (double)(k*m) / (double)s;
                    RE(pfft->w[addr]) = cos(phi);
                    IM(pfft->w[addr]) = sin(phi);
                    addr++;
                }
            }
        }
        s /= n2;
    }

    pfft->t0 = pfft->t0 ? (complex_t*) realloc(pfft->t0, n*sizeof(complex_t)):
                          (complex_t*) malloc(           n*sizeof(complex_t));

    pfft->t1 = pfft->t1 ? (complex_t*) realloc(pfft->t1, n*sizeof(complex_t)):
                          (complex_t*) malloc(           n*sizeof(complex_t));
    pfft->n = n;
    
    /* w32 fill */
    addr = 0;
    for(k = 0; k < 4; k++)
    {
        for(m = 0; m < 8; m++)
        {
            phi = - M_2PI * (double)(k*m) / 32.0;
            RE(pfft->w32[addr]) = cos(phi);
            IM(pfft->w32[addr]) = sin(phi);
            addr++;
        }
    }
    
    
    /* w64 fill */
    addr = 0;
    for(k = 0; k < 8; k++)
    {
        for(m = 0; m < 8; m++)
        {
            phi = - M_2PI * (double)(k*m) / 64.0;
            RE(pfft->w64[addr]) = cos(phi);
            IM(pfft->w64[addr]) = sin(phi);
            addr++;
        }
    }
    
    /* w128 fill */
    addr = 0;
    for(k = 0; k < 8; k++)
    {
        for(m = 0; m < 16; m++)
        {
            phi = - M_2PI * (double)(k*m) / 128.0;
            RE(pfft->w128[addr]) = cos(phi);
            IM(pfft->w128[addr]) = sin(phi);
            addr++;
        }
    }
    
    /* w256 fill */
    addr = 0;
    for(k = 0; k < 16; k++)
    {
        for(m = 0; m < 16; m++)
        {
            phi = - M_2PI * (double)(k*m) / 256.0;
            RE(pfft->w256[addr]) = cos(phi);
            IM(pfft->w256[addr]) = sin(phi);
            addr++;
        }
    }
    
    /* w512 fill */
    addr = 0;
    for(k = 0; k < 16; k++)
    {
        for(m = 0; m < 32; m++)
        {
            phi = - M_2PI * (double)(k*m) / 512.0;
            RE(pfft->w512[addr]) = cos(phi);
            IM(pfft->w512[addr]) = sin(phi);
            addr++;
        }
    }
    
    /* w1024 fill */
    if(pfft->w1024 == NULL)
    {
        pfft->w1024 = (complex_t*) malloc(1024 * sizeof(complex_t));
        addr = 0;
        for(k = 0; k < 32; k++)
        {
            for(m = 0; m < 32; m++)
            {
                phi = - M_2PI * (double)(k*m) / 1024.0;
                RE(pfft->w1024[addr]) = cos(phi);
                IM(pfft->w1024[addr]) = sin(phi);
                addr++;
            }
        }
    }
    
    /* w2048 fill */
    if(pfft->w2048 == NULL)
    {
        pfft->w2048= (complex_t*) malloc(2048 * sizeof(complex_t));
        addr = 0;
        for(k = 0; k < 32; k++)
        {
            for(m = 0; m < 64; m++)
            {
                phi = - M_2PI * (double)(k*m) / 2048.0;
                RE(pfft->w2048[addr]) = cos(phi);
                IM(pfft->w2048[addr]) = sin(phi);
                addr++;
            }
        }
    }
    
    /* w4096 fill */
    if(pfft->w4096 == NULL)
    {
        pfft->w4096= (complex_t*) malloc(4096 * sizeof(complex_t));
        addr = 0;
        for(k = 0; k < 16; k++)
        {
            for(m = 0; m < 256; m++)
            {
                phi = - M_2PI * (double)(k*m) / 4096.0;
                RE(pfft->w4096[addr]) = cos(phi);
                IM(pfft->w4096[addr]) = sin(phi);
                addr++;
            }
        }
    }

    return RES_OK;
error_proc:
    if(pfft->t0) free(pfft->t0);
    if(pfft->t1) free(pfft->t1);
    if(pfft->w)    free(pfft->w);
    pfft->n = 0;
    return err;
}





#ifdef DOXYGEN_ENGLISH
/*! ****************************************************************************
\ingroup DFT_GROUP
\fn void fft_free(fft_t *pfft)
\brief Free `fft_t` structure.

The function clears the intermediate data memory
and vectors of FFT twiddle factors of the structure `fft_t`.

\param[in] pfft
Pointer to the `fft_t` object. \n

\author Sergey Bakhurin www.dsplib.org 
***************************************************************************** */
#endif
#ifdef DOXYGEN_RUSSIAN
/*! ****************************************************************************
\ingroup DFT_GROUP
\fn void fft_free(fft_t *pfft)
\brief Очистить структуру `fft_t` алгоритма БПФ

Функция производит очищение памяти промежуточных данных 
и векторов поворотных коэффициентов структуры `fft_t`.

\param[in] pfft
Указатель на структуру `fft_t`. \n

\author Бахурин Сергей www.dsplib.org 
***************************************************************************** */
#endif
void DSPL_API fft_free(fft_t *pfft)
{
    if(!pfft)
        return;
    if(pfft->w)
        free(pfft->w);
    if(pfft->t0)
        free(pfft->t0);
    if(pfft->t1)
        free(pfft->t1);
      
    if(pfft->w1024)
        free(pfft->w1024);
      
    if(pfft->w2048)
        free(pfft->w2048);
      
    if(pfft->w4096)
        free(pfft->w4096);
      
    memset(pfft, 0, sizeof(fft_t));
}





#ifdef DOXYGEN_ENGLISH

#endif
#ifdef DOXYGEN_RUSSIAN

#endif
int DSPL_API fft_mag(double* x, int n, fft_t* pfft, 
                     double fs, int flag,
                     double* mag, double* freq)
{
    int k, err;
    err = fft_abs(x, n, pfft, fs, flag, mag, freq);
    if(err != RES_OK)
        return err;
    if(mag)
    {
        if(flag & DSPL_FLAG_LOGMAG)
            for(k = 0; k < n; k++)
                mag[k] = 20.0 * log10(mag[k] + DBL_EPSILON);
        else
            for(k = 0; k < n; k++)
                mag[k] *= mag[k];
    }
    
    return err;
}







#ifdef DOXYGEN_ENGLISH

#endif
#ifdef DOXYGEN_RUSSIAN

#endif
int DSPL_API fft_mag_cmplx(complex_t* x, int n, fft_t* pfft, 
                           double fs, int flag,
                           double* mag, double* freq)
{
    int k, err;
    err = fft_abs_cmplx(x, n, pfft, fs, flag, mag, freq);
    if(err != RES_OK)
        return err;
    if(mag)
    {
        if(flag & DSPL_FLAG_LOGMAG)
            for(k = 0; k < n; k++)
                mag[k] = 20.0 * log10(mag[k] + DBL_EPSILON);
        else
            for(k = 0; k < n; k++)
                mag[k] *= mag[k];
    }
    
    return err;
}





#ifdef DOXYGEN_ENGLISH
/*! ****************************************************************************
\ingroup DFT_GROUP
\fn int fft_shift(double* x, int n, double* y)
\brief Perform a shift of the vector `x`, for use with the `fft` and `ifft` 
functions, in order 
 <a href="http://en.dsplib.org/content/dft_freq/dft_freq.html">
 to move the frequency 0 to the center
</a> of the vector `y`.

\param[in]  x
Pointer to the input vector (FFT or IFFT result).  \n
Vector size is `[n x 1]`.  \n \n

\param[in]  n
Input and output vector size. \n \n

\param[out] y
Pointer to the output vector with frequency 0 in the center. \n
Vector size is `[n x 1]`.  \n
Memory must be allocated. \n \n


\return
`RES_OK` if function is calculated successfully. \n
Else \ref ERROR_CODE_GROUP "code error".

\author Sergey Bakhurin www.dsplib.org 
***************************************************************************** */
#endif
#ifdef DOXYGEN_RUSSIAN
/*! ****************************************************************************
\ingroup DFT_GROUP
\fn int fft_shift(double* x, int n, double* y)
\brief Перестановка спектральных отсчетов дискретного преобразования Фурье

Функция производит 
<a href="http://ru.dsplib.org/content/dft_freq/dft_freq.html">
перестановку спектральных отсчетов ДПФ
</a> и переносит нулевую частоту в центр вектора ДПФ.  \n
Данная функция обрабатывает вещественные входные и выходные вектора 
и может применяться для перестановки 
амплитудного или фазового спектра.

\param[in]  x
Указатель на исходный вектор ДПФ.  \n
Размер вектора `[n x 1]`.  \n \n

\param[in]  n
Размер ДПФ \f$n\f$ (размер векторов до и после перестановки). \n \n

\param[out] y
Указатель на вектор результата перестановки. \n
Размер вектора `[n x 1]`.  \n
Память должна быть выделена. \n \n

\return
`RES_OK` если перестановка произведена успешно.  \n
В противном случае \ref ERROR_CODE_GROUP "код ошибки". \n

\author Бахурин Сергей www.dsplib.org 
***************************************************************************** */
#endif
int DSPL_API fft_shift(double* x, int n, double* y)
{
    int n2, r;
    int k;
    double tmp;
    double *buf;

    if(!x || !y)
        return ERROR_PTR;

    if(n<1)
        return ERROR_SIZE;

    r = n%2;
    if(!r)
    {
        n2 = n>>1;
        for(k = 0; k < n2; k++)
        {
            tmp = x[k];
            y[k] = x[k+n2];
            y[k+n2] = tmp;
        }
    }
    else
    {
        n2 = (n+1) >> 1;
        buf = (double*) malloc(n2*sizeof(double));
        memcpy(buf, x, n2*sizeof(double));
        memcpy(y, x+n2, (n2-1)*sizeof(double));
        memcpy(y+n2-1, buf, n2*sizeof(double));
        free(buf);
    }
    return RES_OK;
}



#ifdef DOXYGEN_ENGLISH

#endif
#ifdef DOXYGEN_RUSSIAN

#endif
int DSPL_API fft_shift_cmplx(complex_t* x, int n, complex_t* y)
{
    int n2, r;
    int k;
    complex_t tmp;
    complex_t *buf;

    if(!x || !y)
        return ERROR_PTR;

    if(n<1)
        return ERROR_SIZE;

    r = n%2;
    if(!r)
    {
        n2 = n>>1;
        for(k = 0; k < n2; k++)
        {
            RE(tmp) = RE(x[k]);
            IM(tmp) = IM(x[k]);

            RE(y[k]) = RE(x[k+n2]);
            IM(y[k]) = IM(x[k+n2]);

            RE(y[k+n2]) = RE(tmp);
            IM(y[k+n2]) = IM(tmp);
        }
    }
    else
    {
        n2 = (n+1) >> 1;
        buf = (complex_t*) malloc(n2*sizeof(complex_t));
        memcpy(buf, x, n2*sizeof(complex_t));
        memcpy(y, x+n2, (n2-1)*sizeof(complex_t));
        memcpy(y+n2-1, buf, n2*sizeof(complex_t));
        free(buf);
    }
    return RES_OK;
}

