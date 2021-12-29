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
\fn int idft_cmplx(complex_t* x, int n, complex_t* y)
\brief Inverse discrete Fourier transform of the complex spectrum.

The function calculates the \f$ n \f$ -point inverse discrete transform
Fourier complex spectrum \f$ x (m) \f$, \f$ m = 0 \ldots n-1 \f$. \n
\f[
  y (k) = \sum_ {m = 0} ^ {n-1} x (m)
  \exp \left (j \frac {2 \pi} {n} m k \right),
\f]
where \f$ k = 0 \ldots n-1 \f$.

\param [in] x
Pointer to the vector of the input complex signal spectrum \f$ x (m) \f$,
\f$ m = 0 \ldots n-1 \f$. \n
The size of the vector is `[n x 1]`. \n \n

\param [in] n
The size of the ODPF \f$ n \f$ 
(the size of the vectors of the input spectrum and the result of the ODPF). \n
\n

\param [out] y
Pointer to the complex vector of the ODPF result \f$ y (k) \f$,
\f$ k = 0 \ldots n-1 \f$.
The size of the vector is `[n x 1]`. \n
Memory must be allocated. \n \n


\return
`RES_OK` if the ODPF is calculated successfully. \n
Otherwise, \ref ERROR_CODE_GROUP "error code".

An example of using the `dft_cmplx` function:

\include idft_cmplx_test.c

The result of the program:

\verbatim
x [0] = 0.000 + 0.000j, z [0] = 0.000 -0.000
x [1] = 1.000 + 0.000j, z [1] = 1.000 -0.000
x [2] = 2.000 + 0.000j, z [2] = 2.000 -0.000
x [3] = 3.000 + 0.000j, z [3] = 3.000 -0.000
x [4] = 4.000 + 0.000j, z [4] = 4.000 -0.000
x [5] = 5.000 + 0.000j, z [5] = 5.000 -0.000
x [6] = 6.000 + 0.000j, z [6] = 6.000 -0.000
x [7] = 7.000 + 0.000j, z [7] = 7.000 -0.000
x [8] = 8.000 + 0.000j, z [8] = 8.000 -0.000
x [9] = 9.000 + 0.000j, z [9] = 9.000 -0.000
x [10] = 10.000 + 0.000j, z [10] = 10.000 -0.000
x [11] = 11.000 + 0.000j, z [11] = 11.000 +0.000
x [12] = 12.000 + 0.000j, z [12] = 12.000 +0.000
x [13] = 13.000 + 0.000j, z [13] = 13.000 +0.000
x [14] = 14.000 + 0.000j, z [14] = 14.000 +0.000
x [15] = 15.000 + 0.000j, z [15] = 15.000 -0.000
\endverbatim

\note
This function performs the calculation of the DFT using the naive method.
and requires \f$ n ^ 2 \f$ complex multiplications. \n
To increase the calculation speed, it is recommended
use fast Fourier transform algorithms.

\author Bakhurin Sergey www.dsplib.org
***************************************************************************** */
#endif
#ifdef DOXYGEN_RUSSIAN
/*! ****************************************************************************
\ingroup DFT_GROUP
\fn int idft_cmplx(complex_t* x, int n, complex_t* y)
\brief Обратное дискретное преобразование Фурье комплексного спектра.

Функция рассчитывает \f$ n \f$-точечное обратное дискретное преобразование
Фурье комплексного спектра \f$ x(m) \f$, \f$ m = 0 \ldots n-1 \f$. \n
\f[
  y(k) = \sum_{m = 0}^{n-1} x(m) 
  \exp \left( j   \frac{2\pi}{n}   m   k \right),
\f]
где \f$ k = 0 \ldots n-1 \f$.

\param[in]  x
Указатель на вектор входного комплексного спектра сигнала \f$x(m)\f$, 
\f$ m = 0 \ldots n-1 \f$.  \n
Размер вектора `[n x 1]`.  \n \n

\param[in]  n
Размер ОДПФ \f$n\f$ (размер векторов входного спектра и результата ОДПФ). \n \n

\param[out] y
Указатель на комплексный вектор результата ОДПФ \f$y(k)\f$, 
\f$ k = 0 \ldots n-1 \f$. 
Размер вектора `[n x 1]`.  \n
Память должна быть выделена. \n \n


\return
`RES_OK` если ОДПФ рассчитана успешно.  \n
В противном случае \ref ERROR_CODE_GROUP "код ошибки".

Пример использования функции `dft_cmplx`:

\include idft_cmplx_test.c

Результат работы программы:

\verbatim
x[ 0] =     0.000   +0.000j,    z[ 0] =     0.000   -0.000
x[ 1] =     1.000   +0.000j,    z[ 1] =     1.000   -0.000
x[ 2] =     2.000   +0.000j,    z[ 2] =     2.000   -0.000
x[ 3] =     3.000   +0.000j,    z[ 3] =     3.000   -0.000
x[ 4] =     4.000   +0.000j,    z[ 4] =     4.000   -0.000
x[ 5] =     5.000   +0.000j,    z[ 5] =     5.000   -0.000
x[ 6] =     6.000   +0.000j,    z[ 6] =     6.000   -0.000
x[ 7] =     7.000   +0.000j,    z[ 7] =     7.000   -0.000
x[ 8] =     8.000   +0.000j,    z[ 8] =     8.000   -0.000
x[ 9] =     9.000   +0.000j,    z[ 9] =     9.000   -0.000
x[10] =    10.000   +0.000j,    z[10] =    10.000   -0.000
x[11] =    11.000   +0.000j,    z[11] =    11.000   +0.000
x[12] =    12.000   +0.000j,    z[12] =    12.000   +0.000
x[13] =    13.000   +0.000j,    z[13] =    13.000   +0.000
x[14] =    14.000   +0.000j,    z[14] =    14.000   +0.000
x[15] =    15.000   +0.000j,    z[15] =    15.000   -0.000
\endverbatim

\note
Данная функция выполняет расчет ОДПФ наивным методом 
и требует \f$ n^2 \f$ комплексных умножений. \n
Для увеличения скорости расчета рекомендуется 
использовать алгоритмы быстрого преобразования Фурье.

\author Бахурин Сергей www.dsplib.org 
***************************************************************************** */
#endif
int DSPL_API idft_cmplx(complex_t* x, int n, complex_t* y)
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
            phi    =    M_2PI * divn * (double)k * (double)m;
            RE(e) = cos(phi);
            IM(e) = sin(phi);
            RE(y[k]) += CMRE(x[m], e);
            IM(y[k]) += CMIM(x[m], e);
        }
        RE(y[k]) /= (double)n;
        IM(y[k]) /= (double)n;
    }
    return RES_OK;
}

