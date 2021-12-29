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
