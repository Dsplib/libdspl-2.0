/*
* Copyright (c) 2015-2022 Sergey Bakhurin
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
* MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.    See the
* GNU General Public License for more details.
*
* You should have received a copy of the GNU General Public License
* along with Foobar.    If not, see <http://www.gnu.org/licenses/>.
*/


#include <stdio.h>
#include <stdlib.h>
#include "dspl.h"



#ifdef DOXYGEN_ENGLISH
/*! ****************************************************************************
\ingroup SPEC_MATH_TRANSCEND
\fn int bessel_i0(double* x, int n, double* y)
\brief Modified Bessel Function of the First Kind – \f$ I_0(x)\f$ [1].

\param[in] x
Pointer to the function argument vector \f$ x \f$. \n
Vector size is `[n x 1]`. \n
Input vector must contain nonnegative values. \n \n

\param[in] n
Input vector size `x`. \n \n

\param[out] y
Pointer to \f$ I_0(x)\f$ function vector. \n
Vector size is `[n x 1]`. \n
Memory must be allocated. \n \n

\return
`RES_OK` if function calculated successfully. \n
Else \ref ERROR_CODE_GROUP "code error". \n

\note
[1] Rational Approximations for the Modified Bessel Function 
    of the First Kind – I0(x) for Computations with Double Precision
    by PAVEL HOLOBORODKO on NOVEMBER 11, 2015

Example:

\include bessel_i0.c

Program calcultes \f$ I_0(x)\f$ function for `x` 
in    \f$[0 \ 3]\f$ interval. 
Data saved if `dat/dat0.txt` file and shows on the plot

\image html bessel_i0.png

\author Sergey Bakhurin www.dsplib.org
***************************************************************************** */
#endif
#ifdef DOXYGEN_RUSSIAN
/*! ****************************************************************************
\ingroup SPEC_MATH_TRANSCEND
\fn int bessel_i0(double* x, int n, double* y)
\brief
Модифицированная функция Бесселя первого рода \f$ I_0(x)\f$.

Функция рассчитывает значения функции для вещественного вектора `x`, 
который должен принимать неотрицательные значения. \n

\param[in]  x
Указатель на вектор переменной \f$ x \f$. \n
Размер вектора `[n x 1]`. \n
Память должна быть выделена. \n \n

\param[in]  n
Размер входного вектора `x`. \n \n

\param[out] y
Указатель на вектор значений функции \f$ I_0(x)\f$. \n
Размер вектора `[n x 1]`. \n
Память должна быть выделена. \n \n

\return
`RES_OK` --- расчёт произведен успешно. \n
В противном случае \ref ERROR_CODE_GROUP "код ошибки". \n

\note
Используемый алгоритм описан в статье:
Rational Approximations for the Modified Bessel Function 
of the First Kind – I0(x) for Computations with Double Precision 
by PAVEL HOLOBORODKO on NOVEMBER 11, 2015 

Пример использования функции `bessel_i0`:

\include bessel_i0.c

Данная программа рассчитывает значения функции \f$ I_0(x)\f$ переменной `x` 
в интервале \f$[0 \ 3]\f$. 
Рассчитанные данные сохраняются в текстовый файл `dat/dat0.txt` 
и выводятся на график `img/bessel_i0.png`

\image html bessel_i0.png

\author Бахурин Сергей www.dsplib.org
***************************************************************************** */
#endif
int DSPL_API bessel_i0(double* x, int n, double* y)
{
    double P16[17] = { 1.0000000000000000000000801e+00,
                       2.4999999999999999999629693e-01,
                       2.7777777777777777805664954e-02,
                       1.7361111111111110294015271e-03,
                       6.9444444444444568581891535e-05,
                       1.9290123456788994104574754e-06,
                       3.9367598891475388547279760e-08,
                       6.1511873265092916275099070e-10,
                       7.5940584360755226536109511e-12,
                       7.5940582595094190098755663e-14,
                       6.2760839879536225394314453e-16,
                       4.3583591008893599099577755e-18,
                       2.5791926805873898803749321e-20,
                       1.3141332422663039834197910e-22,
                       5.9203280572170548134753422e-25,
                       2.0732014503197852176921968e-27,
                       1.1497640034400735733456400e-29};

    double P22[23] = { 3.9894228040143265335649948e-01,
                       4.9867785050353992900698488e-02,
                       2.8050628884163787533196746e-02,
                       2.9219501690198775910219311e-02,
                       4.4718622769244715693031735e-02,
                       9.4085204199017869159183831e-02,
                      -1.0699095472110916094973951e-01,
                       2.2725199603010833194037016e+01,
                      -1.0026890180180668595066918e+03,
                       3.1275740782277570164423916e+04,
                      -5.9355022509673600842060002e+05,
                       2.6092888649549172879282592e+06,
                       2.3518420447411254516178388e+08,
                      -8.9270060370015930749184222e+09,
                       1.8592340458074104721496236e+11,
                      -2.6632742974569782078420204e+12,
                       2.7752144774934763122129261e+13,
                      -2.1323049786724612220362154e+14,
                       1.1989242681178569338129044e+15,
                      -4.8049082153027457378879746e+15,
                       1.3012646806421079076251950e+16,
                      -2.1363029690365351606041265e+16,
                       1.6069467093441596329340754e+16};

    double x2;
    int k;

    if(!x || !y)
        return ERROR_PTR;
    if(n < 1)
        return ERROR_SIZE;

    for(k =0; k < n; k++)
    {
        if(x[k] < 0.0)
            return ERROR_NEGATIVE;

        if(x[k] < 7.75)
        {
            x2 = x[k] * x[k] * 0.25;
            polyval(P16, 16, &x2, 1, y+k);
            y[k] = x2 * y[k] + 1.0;
        }
        else
        {
            x2 = 1.0 / x[k];
            polyval(P22, 22, &x2, 1, y+k);
            y[k] *= exp(x[k]) / sqrt(x[k]);
        }
    }
    return RES_OK;
}
