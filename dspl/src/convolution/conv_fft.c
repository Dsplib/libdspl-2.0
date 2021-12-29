/*
* Copyright (c) 2015-2020 Sergey Bakhurin
* Digital Signal Processing Library [http://dsplib.org]
*
* This file is part of DSPL.
*
* is free software: you can redistribute it and/or modify
* it under the terms of the GNU General Public License as published by
* the Free Software Foundation, either version 3 of the License, or
* (at your option) any later version.
*
* DSPL is distributed in the hope that it will be useful,
* but WITHOUT ANY WARRANTY; without even the implied warranty of
* MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
* GNU General Public License for more details.
*
* You should have received a copy of the GNU General Public License
* along with Foobar.  If not, see <http://www.gnu.org/licenses/>.
*/


#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "dspl.h"





#ifdef DOXYGEN_ENGLISH
/*! ****************************************************************************
\ingroup FILTER_CONV_GROUP
\fn int conv_fft(double* a, int na, double* b, int nb,
                 fft_t* pfft, int nfft, double* c) 
\brief Real vectors fast linear convolution by using fast Fourier
transform algorithms

Function convolves two real vectors \f$ c = a * b\f$ length `na` and `nb`
in the frequency domain by using FFT algorithms. This approach provide 
high-performance convolution which increases with `na` and `nb` increasing.
The output convolution is a vector `c` with length equal to  `na + nb - 1`. 

\param[in]  a
Pointer to the first vector `a`. \n
Vector size is `[na x 1]`. \n \n

\param[in]  na
Size of the first vector `a`. \n \n

\param[in]  b
Pointer to the second vector `b`. \n
Vector size is `[nb x 1]`. \n \n

\param[in]  nb
Size of the second vector `b`. \n \n

\param[in]  pfft
Pointer to the structure `fft_t`. \n
Function changes `fft_t` structure fields so `fft_t` must
be clear before program returns. \n \n

\param[in] nfft
FFT size.  \n
This parameter set which FFT size will be used 
for overlapped frequency domain convolution. \n
FFT size must be more of minimal `na` and `nb` value.
For example if `na = 10`, `nb = 4` then `nfft` parameter must 
be more than 4.  \n

\param[out] c
Pointer to the convolution output vector  \f$ c = a * b\f$. \n
Vector size is `[na + nb - 1  x  1]`. \n
Memory must be allocated. \n \n

\return `RES_OK` if convolution is calculated successfully. \n
Else \ref ERROR_CODE_GROUP "code error".  \n \n

Example:
\include conv_fft_test.c

Program output:

\verbatim
conv_fft error: 0x00000000
conv error:     0x00000000
c[  0] =     -0.00    d[  0] =      0.00
c[  1] =     -0.00    d[  1] =      0.00
c[  2] =      1.00    d[  2] =      1.00
c[  3] =      4.00    d[  3] =      4.00
c[  4] =     10.00    d[  4] =     10.00
c[  5] =     20.00    d[  5] =     20.00
c[  6] =     35.00    d[  6] =     35.00
c[  7] =     56.00    d[  7] =     56.00
c[  8] =     77.00    d[  8] =     77.00
c[  9] =     98.00    d[  9] =     98.00
c[ 10] =    119.00    d[ 10] =    119.00
c[ 11] =    140.00    d[ 11] =    140.00
c[ 12] =    161.00    d[ 12] =    161.00
c[ 13] =    182.00    d[ 13] =    182.00
c[ 14] =    190.00    d[ 14] =    190.00
c[ 15] =    184.00    d[ 15] =    184.00
c[ 16] =    163.00    d[ 16] =    163.00
c[ 17] =    126.00    d[ 17] =    126.00
c[ 18] =     72.00    d[ 18] =     72.00
\endverbatim

\author Sergey Bakhurin www.dsplib.org
***************************************************************************** */
#endif
#ifdef DOXYGEN_RUSSIAN
/*! ****************************************************************************
\ingroup FILTER_CONV_GROUP
\fn int conv_fft(double* a, int na, double* b, int nb,
                 fft_t* pfft, double* c) 
\brief Линейная свертка двух вещественных векторов с использованием алгоритмов
быстрого преобразования Фурье

Функция рассчитывает линейную свертку двух векторов \f$ c = a * b\f$ используя
секционную обработку с перекрытием в частотной области. Это позволяет сократить 
вычислительные операции при расчете длинных сверток.

\param[in]  a
Указатель на первый вектор  \f$a\f$. \n 
Размер вектора `[na x 1]`. \n  \n 

\param[in]  na
Размер первого вектора. \n  \n 

\param[in]  b
Указатель на второй вектор \f$b\f$. \n 
Размер вектора `[nb x 1]`. \n  \n 

\param[in]  nb
Размер второго вектора. \n  \n 

\param[in]  pfft
Указатель на структуру `fft_t` алгоритма 
быстрого преобразования Фурье. \n 
Функция изменит состояние полей структуры `fft_t`,
поэтому структура должна быть очищена перед выходом из 
программы для исключения утечек памяти. \n 

\param[in]  nfft
Размер алгоритма БПФ который будет использован для расчета
секционной свертки с перекрытием. \n 
Данный параметр должен быть больше чем минимальное значение
размеров сворачиваемых векторов. \n 
Например если `na=10`, а `nb=4`, то параметр `nfft` должен быть больше 4. \n 
Библиотека поддерживает алгоритмы БПФ составной длины
\f$n = n_0 \times n_1 \times n_2 \times \ldots \times n_p \times m\f$,
где \f$n_i = 2,3,5,7\f$, а \f$m \f$ --- произвольный простой множитель 
не превосходящий 46340 (см. описание функции \ref fft_create).
Однако, максимальное быстродействие достигается при использовании длин равных 
степени двойки.

\param[out] c
Указатель на вектор свертки \f$ c = a * b\f$. \n 
Размер вектора `[na + nb - 1  x  1]`. \n 
Память должна быть выделена. \n  \n 

\return
`RES_OK` если свертка рассчитана успешно. \n 
 В противном случае \ref ERROR_CODE_GROUP "код ошибки".

\note
Данная функция наиболее эффективна при вычислении длинных сверток.

Пример использования функции:

\include conv_fft_test.c

Результат работы:
\verbatim

conv_fft error: 0x00000000
conv error:     0x00000000
c[  0] =     -0.00    d[  0] =      0.00
c[  1] =     -0.00    d[  1] =      0.00
c[  2] =      1.00    d[  2] =      1.00
c[  3] =      4.00    d[  3] =      4.00
c[  4] =     10.00    d[  4] =     10.00
c[  5] =     20.00    d[  5] =     20.00
c[  6] =     35.00    d[  6] =     35.00
c[  7] =     56.00    d[  7] =     56.00
c[  8] =     77.00    d[  8] =     77.00
c[  9] =     98.00    d[  9] =     98.00
c[ 10] =    119.00    d[ 10] =    119.00
c[ 11] =    140.00    d[ 11] =    140.00
c[ 12] =    161.00    d[ 12] =    161.00
c[ 13] =    182.00    d[ 13] =    182.00
c[ 14] =    190.00    d[ 14] =    190.00
c[ 15] =    184.00    d[ 15] =    184.00
c[ 16] =    163.00    d[ 16] =    163.00
c[ 17] =    126.00    d[ 17] =    126.00
c[ 18] =     72.00    d[ 18] =     72.00
\endverbatim

\author Бахурин Сергей. www.dsplib.org
***************************************************************************** */
#endif
int DSPL_API conv_fft(double* a, int na, double* b, int nb,
                      fft_t* pfft, int nfft, double* c)
{
    complex_t *pa = NULL, *pb = NULL, *pc = NULL;
    int err;
    
    if(!a || !b || !c || !pfft)
        return ERROR_PTR;
    if(na<1 || nb < 1)
        return ERROR_SIZE;
    if(nfft<2)
        return ERROR_FFT_SIZE;
    
    pa = (complex_t*) malloc(na*sizeof(complex_t));
    pb = (complex_t*) malloc(nb*sizeof(complex_t));
    pc = (complex_t*) malloc((na+nb-1)*sizeof(complex_t));
    
    re2cmplx(a, na, pa);
    re2cmplx(b, nb, pb);
    
    err = conv_fft_cmplx(pa, na, pb, nb, pfft, nfft, pc);
    if(err != RES_OK)
        goto exit_label;
    
    err = cmplx2re(pc, na+nb-1, c, NULL);
    
exit_label:
    if(pa) free(pa);
    if(pb) free(pb);
    if(pc) free(pc);
    
    return err;
}

