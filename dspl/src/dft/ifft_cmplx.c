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