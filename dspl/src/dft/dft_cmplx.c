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
#include <math.h>
#include "dspl.h"



#ifdef DOXYGEN_ENGLISH
/*! ****************************************************************************
\ingroup DFT_GROUP
\fn int dft_cmplx(complex_t* x, int n, complex_t* y)
\brief Discrete Fourier transform of a complex signal.

The function calculates the \f$ n \f$ -point discrete Fourier transform
complex signal \f$ x (m) \f$, \f$ m = 0 \ldots n-1 \f$. \n
\f[
  Y (k) = \sum_ {m = 0} ^ {n-1} x (m)
  \exp \left (-j \frac {2 \pi} {n} m k \right),
\f]
where \f$ k = 0 \ldots n-1 \f$.

\param [in] x
Pointer to a vector of complex
input signal \f$ x (m) \f$, \f$ m = 0 \ldots n-1 \f$. \n
The size of the vector is `[n x 1]`. \n \n

\param [in] n
The size of the DFT \f$ n \f$ 
(the size of the vectors of the input signal and the result of the DFT). \n \n

\param [out] y
Integrated Vector Pointer
DFT result \f$ Y (k) \f$, \f$ k = 0 \ldots n-1 \f$. \n
The size of the vector is `[n x 1]`. \n
Memory must be allocated. \n \n


\return
`RES_OK` if the DFT is calculated successfully. \n
Otherwise, \ref ERROR_CODE_GROUP "error code".

An example of using the `dft_cmplx` function:

\include dft_cmplx_test.c

The result of the program:

\verbatim
y [0] = 120.000 0.000
y [1] = -8.000 40.219
y [2] = -8.000 19.314
y [3] = -8.000 11.973
y [4] = -8.000 8.000
y [5] = -8.000 5.345
y [6] = -8.000 3.314
y [7] = -8.000 1.591
y [8] = -8.000 0.000
y [9] = -8.000 -1.591
y [10] = -8.000 -3.314
y [11] = -8.000 -5.345
y [12] = -8.000 -8.000
y [13] = -8.000 -11.973
y [14] = -8.000 -19.314
y [15] = -8.000 -40.219
\endverbatim

\note
This function performs the calculation of the DFT by the naive method
and requires \f$ n ^ 2 \f$ complex multiplications. \n
To increase the calculation speed, it is recommended
use fast Fourier transform algorithms.

\author Bakhurin Sergey www.dsplib.org
***************************************************************************** */
#endif
#ifdef DOXYGEN_RUSSIAN
/*! ****************************************************************************
\ingroup DFT_GROUP
\fn int dft_cmplx(complex_t* x, int n, complex_t* y)
\brief Дискретное преобразование Фурье комплексного сигнала.

Функция рассчитывает \f$ n \f$-точечное  дискретное преобразование Фурье 
комплексного сигнала \f$ x(m) \f$, \f$ m = 0 \ldots n-1 \f$. \n
\f[
    Y(k) = \sum_{m = 0}^{n-1} x(m) 
    \exp \left( -j   \frac{2\pi}{n}   m   k \right),
\f]
где \f$ k = 0 \ldots n-1 \f$.

\param[in]  x
Указатель на вектор комплексного 
входного сигнала \f$x(m)\f$, \f$ m = 0 \ldots n-1 \f$.  \n
Размер вектора `[n x 1]`.  \n \n

\param[in]  n
Размер ДПФ \f$n\f$ (размер векторов входного сигнала и результата ДПФ). \n \n

\param[out]  y
Указатель на комплексный вектор 
результата ДПФ \f$Y(k)\f$, \f$ k = 0 \ldots n-1 \f$. \n
Размер вектора `[n x 1]`.  \n
Память должна быть выделена. \n \n


\return
`RES_OK` если ДПФ рассчитана успешно.  \n
В противном случае \ref ERROR_CODE_GROUP "код ошибки".

Пример использования функции `dft_cmplx`:

\include dft_cmplx_test.c

Результат работы программы:

\verbatim
y[ 0] =   120.000    0.000
y[ 1] =    -8.000   40.219
y[ 2] =    -8.000   19.314
y[ 3] =    -8.000   11.973
y[ 4] =    -8.000    8.000
y[ 5] =    -8.000    5.345
y[ 6] =    -8.000    3.314
y[ 7] =    -8.000    1.591
y[ 8] =    -8.000    0.000
y[ 9] =    -8.000   -1.591
y[10] =    -8.000   -3.314
y[11] =    -8.000   -5.345
y[12] =    -8.000   -8.000
y[13] =    -8.000  -11.973
y[14] =    -8.000  -19.314
y[15] =    -8.000  -40.219
\endverbatim

\note
Данная функция выполняет расчет ДПФ наивным методом 
и требует \f$ n^2 \f$ комплексных умножений. \n
Для увеличения скорости расчета рекомендуется 
использовать алгоритмы быстрого преобразования Фурье.

\author Бахурин Сергей www.dsplib.org 
***************************************************************************** */
#endif
int DSPL_API dft_cmplx(complex_t* x, int n, complex_t* y)
{
    int k;
    int m;
    double divn;
    double phi;
    complex_t e;

    if(!x || !y)
        return ERROR_PTR;

    if(n<1)
        return ERROR_SIZE;

    divn = 1.0 / (double)n;

    for(k = 0; k < n; k++)
    {
        RE(y[k]) = IM(y[k]) = 0.0;
        for(m = 0; m < n; m++)
        {
            phi    = -M_2PI * divn * (double)k * (double)m;
            RE(e) = cos(phi);
            IM(e) = sin(phi);
            RE(y[k]) += CMRE(x[m], e);
            IM(y[k]) += CMIM(x[m], e);
        }
    }
    return RES_OK;
}


