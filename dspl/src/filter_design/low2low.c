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




#ifdef DOXYGEN_ENGLISH
/*! ****************************************************************************
\ingroup IIR_FILTER_DESIGN_GROUP
\fn int  low2low(double* b, double* a, int ord, double w0, double w1,
                 double* beta, double* alpha)

Lowpass to lowpass filter frequency transform

Function transforms lowpass filter transfer function \f$ H(s) \f$ 
to the lowpass filter transfer function \f$ F(s) \f$ 
with other cutoff frequency.

Filter order, magnitude ripple in passband and stopband 
supression still the same.

\param[in]  b
Pointer to the input lowpass filter transfer function \f$H(s)\f$ numerator 
coefficients vector. \n
Vector size is `[ord+1 x 1]`. \n 
\n

\param[in]  a
Pointer to the input lowpass filter transfer function \f$H(s)\f$ denominator 
coefficients vector. \n
Vector size is `[ord+1 x 1]`. \n 
\n

\param[in]  ord
Filter order. \n 
\n

\param[in]  w0
Input lowpass filter cutoff frequency. \n 
\n

\param[in]  w1
Lowpass filter cutoff frequency after transformation. \n 
\n

\param[in,out]  beta
Pointer to the lowpass filter transfer function \f$F(s)\f$ numerator 
coefficients vector after transformation. \n
Vector size is `[ord+1 x 1]`. \n 
Memory must be allocated. \n
\n

\param[in,out]  alpha
Pointer to the lowpass filter transfer function \f$F(s)\f$ denominator 
coefficients vector after transformation. \n
Vector size is `[ord+1 x 1]`. \n 
Memory must be allocated. \n
\n

\return
`RES_OK` if filter coefficients is calculated successfully. \n 
Else \ref ERROR_CODE_GROUP "code error".

\author Sergey Bakhurin www.dsplib.org
***************************************************************************** */
#endif
#ifdef DOXYGEN_RUSSIAN
/*! ****************************************************************************
\ingroup IIR_FILTER_DESIGN_GROUP
\fn int  low2low(double* b, double* a, int ord, double w0, double w1,
                 double* beta, double* alpha)
\brief  Частотное преобразование ФНЧ-ФНЧ

Функция производит преобразование передаточной функции \f$ H(s) \f$ 
аналогового ФНЧ с частотой среза `w0` рад/c 
в передаточную функцию \f$ F(s) \f$ аналоговго ФНЧ с частотой среза `w1` рад/c.

Неравномерность АЧХ в полосе пропускания, уровень подавления в полосе 
заграждения и порядок фильтра остаются неизменными.

\param[in]  b
Указатель на вектор коэффициентов числителя передаточной функции \f$H(s)\f$
исходного аналогового ФНЧ. \n
Размер вектора `[ord+1 x 1]`. \n
Память должна быть выделена. \n 
\n

\param[in]  a
Указатель на вектор коэффициентов знаменателя передаточной функции \f$H(s)\f$
исходного аналогового ФНЧ. \n
Размер вектора `[ord+1 x 1]`. \n
Память должна быть выделена. \n 
\n

\param[in]  ord
Порядок исходного фильтра и фильтра после преобразования. \n 
\n

\param[in]  w0
Частота среза исходного ФНЧ. \n 
\n

\param[in]  w1
Требуемая частота среза ФНЧ после преобразования. \n 
\n

\param[in,out]  beta   Указатель на вектор коэффициентов числителя 
передаточной функции \f$F(s)\f$ ФНЧ после преобразования. \n
Размер вектора `[ord+1 x 1]`. \n
Память должна быть выделена. \n 
\n

\param[in,out]  alpha
Указатель на вектор коэффициентов знаменателя передаточной функции \f$F(s)\f$
аналогового ФНЧ после преобразования. \n
Размер вектора `[ord+1 x 1]`. \n
Память должна быть выделена. \n 
\n

\return
`RES_OK` --- Преоборазование расчитано успешно. \n
В противном случае \ref ERROR_CODE_GROUP "код ошибки". \n

\author Бахурин Сергей www.dsplib.org
***************************************************************************** */
#endif
int DSPL_API low2low(double* b, double* a, int ord, double w0, double w1,
                     double* beta, double* alpha)
{

    double num[2] = {0.0, 1.0};
    double den[2] = {0.0, 0.0};

    if(!b || !a || !beta || !alpha)
        return ERROR_PTR;
    if(ord < 1)
        return ERROR_FILTER_ORD;
    if(w0 <= 0.0 || w1 <= 0.0)
        return ERROR_FILTER_FT;

    den[0] = w1 / w0;

    return ratcompos(b, a, ord, num, den, 1, beta, alpha);
}

