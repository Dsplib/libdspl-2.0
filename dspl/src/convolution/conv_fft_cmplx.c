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
#include "dspl.h"


#ifdef DOXYGEN_ENGLISH
/*! ****************************************************************************
\ingroup FILTER_CONV_GROUP
\fn int conv_fft_cmplx(complex_t* a, int na, complex_t* b, int nb,
                       fft_t* pfft, int nfft, complex_t* c) 
\brief Complex vectors fast linear convolution by using fast Fourier
transform algorithms

Function convolves two complex vectors \f$ c = a * b\f$ length `na` and `nb`
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
\include conv_fft_cmplx_test.c

Program output:

\verbatim
c[  0] =     -1.00    -0.00j    d[  0] =     -1.00    +0.00j
c[  1] =     -6.00    +4.00j    d[  1] =     -6.00    +4.00j
c[  2] =    -15.00   +20.00j    d[  2] =    -15.00   +20.00j
c[  3] =    -28.00   +56.00j    d[  3] =    -28.00   +56.00j
c[  4] =    -45.00  +120.00j    d[  4] =    -45.00  +120.00j
c[  5] =    -55.00  +210.00j    d[  5] =    -55.00  +210.00j
c[  6] =    -65.00  +300.00j    d[  6] =    -65.00  +300.00j
c[  7] =    -75.00  +390.00j    d[  7] =    -75.00  +390.00j
c[  8] =    -85.00  +480.00j    d[  8] =    -85.00  +480.00j
c[  9] =    -95.00  +570.00j    d[  9] =    -95.00  +570.00j
c[ 10] =   -105.00  +660.00j    d[ 10] =   -105.00  +660.00j
c[ 11] =   -115.00  +750.00j    d[ 11] =   -115.00  +750.00j
c[ 12] =   -125.00  +840.00j    d[ 12] =   -125.00  +840.00j
c[ 13] =   -135.00  +930.00j    d[ 13] =   -135.00  +930.00j
c[ 14] =   -145.00 +1020.00j    d[ 14] =   -145.00 +1020.00j
c[ 15] =   -124.00 +1080.00j    d[ 15] =   -124.00 +1080.00j
c[ 16] =    -99.00 +1016.00j    d[ 16] =    -99.00 +1016.00j
c[ 17] =    -70.00  +820.00j    d[ 17] =    -70.00  +820.00j
c[ 18] =    -37.00  +484.00j    d[ 18] =    -37.00  +484.00j
\endverbatim

\author Sergey Bakhurin www.dsplib.org
***************************************************************************** */
#endif
#ifdef DOXYGEN_RUSSIAN
/*! ****************************************************************************
\ingroup FILTER_CONV_GROUP
\fn int conv_fft_cmplx(complex_t* a, int na, complex_t* b, int nb,
                       fft_t* pfft, complex_t* c) 
\brief Линейная свертка двух комплексных векторов с использованием алгоритмов
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

\include conv_fft_cmplx_test.c

Результат работы:
\verbatim
c[  0] =     -1.00    -0.00j    d[  0] =     -1.00    +0.00j
c[  1] =     -6.00    +4.00j    d[  1] =     -6.00    +4.00j
c[  2] =    -15.00   +20.00j    d[  2] =    -15.00   +20.00j
c[  3] =    -28.00   +56.00j    d[  3] =    -28.00   +56.00j
c[  4] =    -45.00  +120.00j    d[  4] =    -45.00  +120.00j
c[  5] =    -55.00  +210.00j    d[  5] =    -55.00  +210.00j
c[  6] =    -65.00  +300.00j    d[  6] =    -65.00  +300.00j
c[  7] =    -75.00  +390.00j    d[  7] =    -75.00  +390.00j
c[  8] =    -85.00  +480.00j    d[  8] =    -85.00  +480.00j
c[  9] =    -95.00  +570.00j    d[  9] =    -95.00  +570.00j
c[ 10] =   -105.00  +660.00j    d[ 10] =   -105.00  +660.00j
c[ 11] =   -115.00  +750.00j    d[ 11] =   -115.00  +750.00j
c[ 12] =   -125.00  +840.00j    d[ 12] =   -125.00  +840.00j
c[ 13] =   -135.00  +930.00j    d[ 13] =   -135.00  +930.00j
c[ 14] =   -145.00 +1020.00j    d[ 14] =   -145.00 +1020.00j
c[ 15] =   -124.00 +1080.00j    d[ 15] =   -124.00 +1080.00j
c[ 16] =    -99.00 +1016.00j    d[ 16] =    -99.00 +1016.00j
c[ 17] =    -70.00  +820.00j    d[ 17] =    -70.00  +820.00j
c[ 18] =    -37.00  +484.00j    d[ 18] =    -37.00  +484.00j
\endverbatim

\author Бахурин Сергей. www.dsplib.org
***************************************************************************** */
#endif
int DSPL_API conv_fft_cmplx(complex_t* a, int na, complex_t* b, int nb,
                            fft_t* pfft,    int nfft, complex_t* c)
{
    
    int La, Lb, Lc, Nz, n, p0, p1, ind, err;
    complex_t *pa, *pb;
    complex_t *pt, *pA, *pB, *pC;
    
    if(!a || !b || !c)
        return ERROR_PTR;
    if(na < 1 || nb < 1)
        return ERROR_SIZE;
    
    if(na >= nb)
    {
        La = na;
        Lb = nb;
        pa = a; 
        pb = b;
    }
    else
    {
        La = nb;
        pa = b;
        Lb = na;
        pb = a;
    }
        
    Lc = La + Lb - 1;
    Nz = nfft - Lb;

    if(Nz <= 0)
        return ERROR_FFT_SIZE;

    pt = (complex_t*)malloc(nfft*sizeof(complex_t));
    pB = (complex_t*)malloc(nfft*sizeof(complex_t));    
    pA = (complex_t*)malloc(nfft*sizeof(complex_t));    
    pC = (complex_t*)malloc(nfft*sizeof(complex_t));

    memset(pt, 0, nfft*sizeof(complex_t));
    memcpy(pt+Nz, pb, Lb*sizeof(complex_t));

    err = fft_cmplx(pt, nfft, pfft, pB);
    if(err != RES_OK)
        goto exit_label;

    p0 = -Lb;
    p1 = p0 + nfft;
    ind = 0;
    while(ind < Lc)
    {
        if(p0 >=0)
        {
            if(p1 < La)
                err = fft_cmplx(pa + p0, nfft, pfft, pA);
            else
            {
                memset(pt, 0, nfft*sizeof(complex_t));
                memcpy(pt, pa+p0, (nfft+La-p1)*sizeof(complex_t));
                err = fft_cmplx(pt, nfft, pfft, pA);
            }
        }
        else
        {
            memset(pt, 0, nfft*sizeof(complex_t));
            if(p1 < La)                
                memcpy(pt - p0, pa, (nfft+p0)*sizeof(complex_t));
            else
                memcpy(pt - p0, pa, La * sizeof(complex_t));
            err = fft_cmplx(pt, nfft, pfft, pA);
        }
        
        if(err != RES_OK)
            goto exit_label;

        for(n = 0; n < nfft; n++)
        {
            RE(pC[n]) = CMRE(pA[n], pB[n]);
            IM(pC[n]) = CMIM(pA[n], pB[n]);
        }


        if(ind+nfft < Lc)
            err = ifft_cmplx(pC, nfft, pfft, c+ind);
        else
        {
            err = ifft_cmplx(pC, nfft, pfft, pt);
            memcpy(c+ind, pt, (Lc-ind)*sizeof(complex_t));
        }
        if(err != RES_OK)
            goto exit_label;
        
        p0  += Nz;
        p1  += Nz;
        ind += Nz;
    }
 
exit_label: 
    if(pt) free(pt);
    if(pB) free(pB);
    if(pA) free(pA);
    if(pC) free(pC);
    
    return err;
}
