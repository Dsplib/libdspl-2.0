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

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "dspl.h"
#include "dspl_internal.h"



#ifdef DOXYGEN_ENGLISH
/*! ****************************************************************************
\ingroup IIR_FILTER_DESIGN_GROUP
\fn int bilinear(double* bs, double* as, int ord, double* bz, double* az)

\brief
Transform a s-plane analog filter transfer function \f$H(s)\f$ to the
digital filter transfer function \f$H(z)\f$.

Bilinear transform is rational composition:

\f[
s \leftarrow \frac{1 - z^{-1}}{1 - z^{-1}}.
\f]

Digital filter order, passband magnitude ripple and stopband suppression 
still the same after bilinear transform as analog filter. 

Frequency \f$\Omega\f$ of analog filter and frequency    
\f$\omega\f$ of digital filter relations:

\f[
\Omega = \tan(\omega / 2).
\f]


\param[in]    bs
Pointer to the vector of analog filter \f$H(s)\f$ 
numerator coefficients.
Vector size is `[ord+1 x 1]`. \n
\n

\param[in]    as
Pointer to the vector of analog filter \f$H(s)\f$ 
denominator coefficients vector.
Vector size is `[ord+1 x 1]`. \n
\n

\param[in]    ord
Analog and digital filters order. \n 
\n

\param[out]    bz
Pointer to the vector of digital filter \f$H(z)\f$ 
numerator coefficients after bilinear transform.
Vector size is `[ord+1 x 1]`. \n
Memory must be allocated. \n 
\n

\param[out]    az
Pointer to the vector of digital filter \f$H(z)\f$ 
denominator coefficients after bilinear transform.
Vector size is `[ord+1 x 1]`. \n
Memory must be allocated. \n 
\n

\return
`RES_OK` if bilinear transform is calculated successfully. \n 
Else \ref ERROR_CODE_GROUP "code error".


Example:

\include bilinear_test.c

This program calculates the transfer function \f$H(s)\f$ of analog 
Chebyshev filter of the first kind, with a cutoff frequency of 1 rad/s, 
and produces bilinear trandform to digital filter, 
with a normilized cutoff frequency equals 0.5.

Result:

\verbatim
bz[0] =     0.246        az[0] =     4.425
bz[1] =     0.983        az[1] =    -3.318
bz[2] =     1.474        az[2] =     4.746
bz[3] =     0.983        az[3] =    -2.477
bz[4] =     0.246        az[4] =     1.034
err = 0
\endverbatim

In addition, the frequency response of the resulting digital filter 
is calculated and plotted by GNUPLOT package.

\image html bilinear.png

\author Sergey Bakhurin www.dsplib.org
***************************************************************************** */
#endif
#ifdef DOXYGEN_RUSSIAN
/*! ****************************************************************************
\ingroup IIR_FILTER_DESIGN_GROUP
\fn int bilinear(double* bs, double* as, int ord, double* bz, double* az)

\brief
Билинейное преобразование передаточной характеристики аналогового 
фильтра \f$H(s)\f$, в передаточную характеристику цифрового фильтра \f$H(z)\f$.

Функция рассчитывает коэффициенты передаточной характеристики \f$H(z)\f$
цифрового фильтра путем дробно-рациональной подстановки вида 

\f[
s \leftarrow \frac{1 - z^{-1}}{1 - z^{-1}}.
\f]

Порядок цифрового фильтра при этом остается равным порядку аналогового фильтра,
а ось частот \f$\Omega\f$ аналогового фильтра связана c осью частот 
\f$\omega\f$ цифрового фильтра соотношением:

\f[
\Omega = \tan(\omega / 2).
\f]



\param[in] bs
Указатель на вектор коэффициентов числителя передаточной функции \f$H(s)\f$ 
исходного аналогового фильтра. \n
Размер вектора `[ord+1 x 1]`. \n
Память должна быть выделена. \n 
\n

\param[in] as
Указатель на вектор коэффициентов знаменателя передаточной функции \f$H(s)\f$ 
исходного аналогового фильтра. \n
Размер вектора `[ord+1 x 1]`. \n
Память должна быть выделена. \n 
\n

\param[in] ord
Порядок фильтра. \n
Количество коэффициентов числителя и знаменателя передаточных функций 
\f$H(s)\f$ и \f$H(z)\f$ аналогового и цифрового фильтров равно `ord+1`. \n 
\n

\param[out] bz
Указатель на вектор коэффициентов числителя передаточной функции \f$H(z)\f$ 
полученного цифрового фильтра. \n
Размер вектора `[ord+1 x 1]`. \n
Память должна быть выделена. \n 
\n

\param[out] az
Указатель на вектор коэффициентов знаменателя передаточной функции \f$H(z)\f$ 
полученного цифрового фильтра. \n
Размер вектора `[ord+1 x 1]`. \n
Память должна быть выделена. \n 
\n

\return
`RES_OK` --- фильтр рассчитан успешно. \n 
В противном случае \ref ERROR_CODE_GROUP "код ошибки". \n


Пример использования функции `bilinear`:

\include bilinear_test.c

Данная программа производит расчет передаточной характеристики аналогового 
фильтра Чебышева первого рода, с частотой среза равной 1 рад/с, и производит 
билинейное преобразование в цифровой, с частотой среза равной 0.5.

Результат работы программы:

\verbatim
bz[0] =     0.246        az[0] =     4.425
bz[1] =     0.983        az[1] =    -3.318
bz[2] =     1.474        az[2] =     4.746
bz[3] =     0.983        az[3] =    -2.477
bz[4] =     0.246        az[4] =     1.034
err = 0
\endverbatim

Кроме этого производится расчет АЧХ полученного цифрового фильтра и строится
график АЧХ пакетом GNUPLOT

\image html bilinear.png

\author Бахурин Сергей www.dsplib.org
***************************************************************************** */
#endif
int DSPL_API bilinear(double* bs, double* as, int ord, double* bz, double* az)
{
    double c[2] = {1.0, -1.0};
    double d[2] = {1.0,  1.0};
    return ratcompos(bs, as, ord, c, d, 1, bz, az);
}


