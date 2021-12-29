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

