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
\fn int butter_ap(double Rp, int ord, double* b, double* a)

\brief
Function calculates the transfer function \f$ H(s) \f$ coefficients of
analog normalized lowpass Butterworth filter.

Analog normalized lowpass filter magnitude ripple equals \f$ -R_p \f$ dB 
for angular frequency \f$ \omega \f$ from 0 to 1 rad/s.

\param[in]  Rp
Magnitude ripple in passband (dB). \n
This parameter sets maximum filter distortion from 0 to 1 rad/s frequency. \n
Parameter must be positive. \n 
\n

\param[in]  ord
Filter order. \n
Filter coefficients number equals `ord+1` for numerator and denominator
of transfer function \f$ H(s) \f$ \n 
\n

\param[out]  b
Pointer to the vector of transfer function \f$H(s)\f$ 
numerator coefficient. \n
Vector size is `[ord+1 x 1]`. \n
Memory must be allocated. \n 
\n

\param[out]  a
Pointer to the vector of transfer function \f$H(s)\f$ 
denominator coefficient. \n
Vector size is `[ord+1 x 1]`. \n
Memory must be allocated. \n 
\n

\return
`RES_OK` if filter coefficients is calculated successfully. \n 
Else \ref ERROR_CODE_GROUP "code error".
\n

Example:

\include butter_ap_test.c

Result:

\verbatim
b[ 0] =     1.965     a[ 0] =     1.965
b[ 1] =     0.000     a[ 1] =     3.138
b[ 2] =     0.000     a[ 2] =     2.505
b[ 3] =     0.000     a[ 3] =     1.000 
\endverbatim
\n

In `dat` folder will be created 3 files: \n

\verbatim
butter_ap_test_mag.txt    magnitude
butter_ap_test_phi.txt    phase response
butter_ap_test_tau.txt    group delay
\endverbatim

In addition, GNUPLOT will build the following graphs from data stored in files:

\image html butter_ap_test.png

\author Sergey Bakhurin www.dsplib.org
***************************************************************************** */
#endif
#ifdef DOXYGEN_RUSSIAN
/*! ****************************************************************************
\ingroup IIR_FILTER_DESIGN_GROUP
\fn int butter_ap(double Rp, int ord, double* b, double* a)

\brief
Расчет передаточной характеристики \f$ H(s) \f$ аналогового 
нормированного ФНЧ Баттерворта.

Функция рассчитывает коэффициенты передаточной характеристики \f$H(s)\f$
аналогового нормированного ФНЧ Баттерворта порядка `ord` с частотой среза 
1 рад/с по уровню \f$ -R_p \f$ дБ.

\param[in]  Rp
Неравномерность АЧХ в полосе пропускания (дБ). \n
Параметр задает уровень искажений в полосе от 0 до 1 рад/с. \n
Значение должно быть положительным. \n 
\n

\param[in]  ord
Порядок фильтра. \n
Количество коэффициентов числителя и знаменателя 
передаточной функции \f$H(s)\f$ равно `ord+1`. \n 
\n

\param[out]  b
Указатель на вектор коэффициентов числителя 
передаточной функции \f$H(s)\f$. \n
Размер вектора `[ord+1 x 1]`. \n
Память должна быть выделена. \n 
\n

\param[out]  a
Указатель на вектор коэффициентов знаменателя 
передаточной функции \f$H(s)\f$. \n
Размер вектора `[ord+1 x 1]`. \n
Память должна быть выделена. \n 
\n

\return
`RES_OK` --- фильтр рассчитан успешно. \n 
В противном случае \ref ERROR_CODE_GROUP "код ошибки". \n 
\n

Пример использования функции `butter_ap`:

\include butter_ap_test.c

Результат работы программы:

\verbatim
b[ 0] =     1.965     a[ 0] =     1.965
b[ 1] =     0.000     a[ 1] =     3.138
b[ 2] =     0.000     a[ 2] =     2.505
b[ 3] =     0.000     a[ 3] =     1.000 
\endverbatim
\n

В каталоге `dat` будут созданы три файла: \n

\verbatim
butter_ap_test_mag.txt    АЧХ фильтра   
butter_ap_test_phi.txt    ФЧХ фильтра
butter_ap_test_tau.txt    ГВЗ фильтра
\endverbatim

Кроме того программа GNUPLOT произведет построение следующих графиков 
по сохраненным в файлах данным:

\image html butter_ap_test.png

\author Бахурин Сергей www.dsplib.org
***************************************************************************** */
#endif
int DSPL_API butter_ap(double rp, int ord, double* b, double* a)
{
    int res;
    complex_t *z = NULL;
    complex_t *p = NULL;
    int nz, np;

    if(rp < 0.0)
        return ERROR_FILTER_RP;
    if(ord < 1)
        return ERROR_FILTER_ORD;
    if(!a || !b)
        return ERROR_PTR;

    z = (complex_t*) malloc(ord*sizeof(complex_t));
    p = (complex_t*) malloc(ord*sizeof(complex_t));


    res = butter_ap_zp(ord, rp, z, &nz, p, &np);
    if(res != RES_OK)
        goto exit_label;

    res = filter_zp2ab(z, nz, p, np, ord, b, a);
    if(res != RES_OK)
        goto exit_label;

    b[0] = a[0];


exit_label:
    if(z)
        free(z);
    if(p)
        free(p);
    return res;
}

