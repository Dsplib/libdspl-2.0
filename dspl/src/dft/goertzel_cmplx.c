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
\ingroup DFT_GROUP
\fn int goertzel_cmplx(complex_t *x, int n, int *ind, int k, complex_t *y)
\brief <a href = "http://en.dsplib.org/content/goertzel/goertzel.html">
Goertzel algorithm </a> individual DFT samples calculation for 
the complex input vector `x`.

Goertzel algorithm calculates `k` samples of `n`-point DFT, according to 
`ind` indexes vector.

\param[in]  x
Pointer to the complex input vector `x` \n
Vector size is `[n x 1]`. \n \n

\param[in]  n
Size of vector `x`. \n \n

\param[in]  ind
Pointer to the DFT samples indexes which need 
to calculate by Goertzel algorithm. \n
Vector size is `[k x 1]`. \n \n

\param[in]  k
Size of vector `ind`. \n \n

\param[out]  y
Pointer to the DFT samples vector corresponds to indexes `ind`. \n
Vector size is `[k x 1]`. \n
Memory must be allocated. \n \n

\return
`RES_OK` if function is calculated successfully. \n 
Else \ref ERROR_CODE_GROUP "code error".

\note
Goertzel's algorithm is effective when it is necessary to calculate 
several DFT samples of a signal of long duration. \n
However, the size `k` of the vector of indices` ind` can be arbitrary, 
including more than the length of the signal `n`.
In this case, some DFT samples will be repeated, but this will not entail 
a runtime error. \n
The values of the indices of the DFT spectral samples `ind` 
can also be arbitrary integers, including negative ones.
In this case, the DFT samples will be calculated.
with indices modulo `n`. \n

\author Sergey Bakhurin www.dsplib.org
***************************************************************************** */
#endif
#ifdef DOXYGEN_RUSSIAN
/*! ****************************************************************************
\ingroup DFT_GROUP
\fn int goertzel_cmplx(complex_t *x, int n, int *ind, int k, complex_t *y)
\brief <a href = "http://ru.dsplib.org/content/goertzel/goertzel.html">
Алгоритм Гёрцеля</a> для расчета отдельных спектральных отсчетов дискретного 
преобразования Фурье комплексного сигнала `x`.

Данный алгоритм позволяет рассчитать `k` спектральных отсчетов 
`n`-точечного ДПФ, заданных вектором индексов `ind`.

\param[in]  x
Указатель на вектор комплексного входного сигнала. \n
Размер вектора `[n x 1]`. \n \n

\param[in]  n
Размер вектора входного сигнала. \n \n

\param[in]  ind
Указатель на вектор индексов спектральных отсчетов для расчета которых
будет использоваться алгоритм Герцеля. \n
Размер вектора `[k x 1]`. \n \n

\param[in]  k
Размер вектора индексов спектральных отсчетов `ind`. \n \n

\param[out]  y
Указатель на вектор спектральных отсчетов, соответствующих номерам `ind`. \n
Размер вектора `[k x 1]`. \n
Память должна быть выделена. \n \n

\return
`RES_OK` --- функция выполнена успешно.  \n
В противном случае \ref ERROR_CODE_GROUP "код ошибки". \n

\note
Алгоритм Герцеля эффективен когда необходимо рассчитать несколько 
спектральных отсчетов сигнала большой длительности.  \n
Однако, размер `k` вектора индексов `ind` может быть произвольным,
в том числе больше длины сигнала `n`. 
В этом случае некоторые спектральные отсчеты 
будут повторяться, но это не повлечет за собой ошибки выполнения. \n
Значения индексов спектральных отсчетов `ind` также могут быть 
произвольными целыми, в том числе и отрицательными. 
В этом случае будут рассчитаны спектральные отсчеты с индексами
по модулю `n`. \n \n

\author Бахурин Сергей www.dsplib.org
***************************************************************************** */
#endif
int DSPL_API goertzel_cmplx(complex_t *x, int n, int *ind, int k, complex_t *y)
{

    int m, p;
    complex_t w;
    double alpha;
    complex_t v[3];

    if(!x || !y || !ind)
        return ERROR_PTR;

    if(n < 1 || k < 1)
        return ERROR_SIZE;

    for(p = 0; p < k; p++)
    {
        RE(w) = cos(M_2PI * (double)ind[p] / (double)n);
        IM(w) = sin(M_2PI * (double)ind[p] / (double)n);

        alpha = 2.0 * RE(w);
        memset(v, 0, 3*sizeof(complex_t));

        for(m = 0; m < n; m++)
        {
            RE(v[2]) = RE(v[1]);
            RE(v[1]) = RE(v[0]);
            RE(v[0]) = RE(x[m]) + alpha * RE(v[1]) - RE(v[2]);

            IM(v[2]) = IM(v[1]);
            IM(v[1]) = IM(v[0]);
            IM(v[0]) = IM(x[m]) + alpha * IM(v[1]) - IM(v[2]);
        }

        RE(y[p]) = CMRE(w, v[0]) - RE(v[1]);
        IM(y[p]) = CMIM(w, v[0]) - IM(v[1]);
    }

    return RES_OK;
}

