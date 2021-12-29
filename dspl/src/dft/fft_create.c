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
