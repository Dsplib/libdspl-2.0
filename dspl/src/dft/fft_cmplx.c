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

