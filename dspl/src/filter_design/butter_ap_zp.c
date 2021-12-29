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
\fn int butter_ap_zp(int ord, double rp, complex_t* z, int* nz, 
                     complex_t* p, int* np)

\brief
Function calculates arrays of zeros and poles for analog normlized lowpass 
Batterworth filter transfer function \f$ H(s) \f$ order `ord` .

Analog normalized lowpass filter magnitude ripple equals \f$ -R_p \f$ dB 
for angular frequency \f$ \omega \f$ from 0 to 1 rad/s.


\param[in]  ord
Filter order. \n
Number of zeros and poles of filter can be less or equal `ord`. \n
\n

\param[in]  rp
Magnitude ripple in passband (dB). \n
This parameter sets maximum filter distortion from 0 to 1 rad/s frequency. \n
Parameter must be positive. \n 
\n

\param[out]  z
Pointer to the \f$ H(s) \f$ zeros array. \n
Maximum vector size is `[ord x 1]`. \n
Memory must be allocated for maximum vector size. \n 
\n

\param[out]  nz
Pointer to the variable which keep number of finite zeros \f$ H(s) \f$. \n
Number of finite zeros which was calculated and saved in vector `z`. \n
Pointer cannot be `NULL`. \n
\n

\param[out]  p
Pointer to the \f$ H(s) \f$ poles array. \n
Maximum vector size is `[ord x 1]`. \n
Memory must be allocated for maximum vector size. \n 
\n

\param[out]  np
Pointer to the variable which keep number of 
calculated poles of \f$ H(s) \f$. \n
Pointer cannot be `NULL`. \n
\n

\return
`RES_OK` if zeros and poles is calculated successfully. \n 
Else \ref ERROR_CODE_GROUP "code error".
\n

\note
Normalized Butterworth lowpass filter has no finite zeros. 
So `z` vector will not changed and in pointer `nz` will write 0 value. \n


Example of normalized Butterworth lowpass filter zeros and poles calculation:
\include butter_ap_zp_test.c

Result:

\verbatim
Butterworth filter zeros: 0
Butterworth filter poles: 7
p[ 0] =    -1.101    +0.000 j
p[ 1] =    -0.245    +1.074 j
p[ 2] =    -0.245    -1.074 j
p[ 3] =    -0.687    +0.861 j
p[ 4] =    -0.687    -0.861 j
p[ 5] =    -0.992    +0.478 j
p[ 6] =    -0.992    -0.478 j  
\endverbatim
\n

In `dat` folder will be created `butter_ap_zp.txt` file. \n

In addition, GNUPLOT will build the following graphs 
from data stored in `dat/butter_ap_zp.txt` file:

\image html butter_ap_zp_test.png

\author Sergey Bakhurin www.dsplib.org
***************************************************************************** */
#endif
#ifdef DOXYGEN_RUSSIAN
/*! ****************************************************************************
\ingroup IIR_FILTER_DESIGN_GROUP
\fn int butter_ap_zp(int ord, double rp, complex_t* z, int* nz, 
                     complex_t* p, int* np)

\brief
Расчет массивов нулей и полюсов передаточной функции 
\f$ H(s) \f$ аналогового нормированного ФНЧ Баттерворта.

Функция рассчитывает значения нулей и полюсов передаточной функции 
\f$ H(s)\f$ аналогового нормированного ФНЧ Баттерворта порядка `ord` 
с частотой среза 1 рад/с по уровню \f$-R_p\f$ дБ. \n 


\param[in]  ord
Порядок фильтра. \n
\n

\param[in]  rp
Неравномерность АЧХ в полосе пропускания (дБ). \n
Параметр задает уровень искажений в полосе от 0 до 1 рад/с. \n
Значение должно быть положительным. \n 
\n

\param[out]  z
Указатель на массив комплексных нулей 
передаточной характеристики \f$ H(s)\f$. \n
Максимальный размер вектора  вектора `[ord x 1]`. \n
Память должна быть выделена. \n 
\n

\param[out]  nz
Указатель на переменную количества нулей 
передаточной характеристики \f$ H(s)\f$. \n
По данному указателю будет записано количество 
нулей фильтра, которые были рассчитаны и 
помещены в вектор `z`. \n
Память должна быть выделена. \n
\n

\param[out]  p
Указатель на массив комплексных полюсов 
передаточной характеристики \f$ H(s)\f$. \n
Максимальный размер вектора  вектора `[ord x 1]`. \n
Память должна быть выделена. \n
\n

\param[out]  np
Указатель на переменную количества полюсов 
передаточной характеристики \f$ H(s)\f$. \n
По данному укащзателю будет записано количество нулей фильтра, которые 
были рассчитны и помещены в вектор `p`. \n
Память должна быть выделена. \n
\n

\return
`RES_OK` --- массивы нулей и полюсов рассчитаны успешно. \n
В противном случае \ref ERROR_CODE_GROUP "код ошибки". \n
\n

\note
Нормированный ФНЧ Баттерворта не имеет нулей, поэтому массив нулей `z` 
не будет изменен, а по указателю `nz` будет записан 0. \n


Пример программы рассчета нулей и полюсов нормированного ФНЧ Баттерворта:
\include butter_ap_zp_test.c

Результат выполнения программы:

\verbatim
Butterworth filter zeros: 0
Butterworth filter poles: 7
p[ 0] =    -1.101    +0.000 j
p[ 1] =    -0.245    +1.074 j
p[ 2] =    -0.245    -1.074 j
p[ 3] =    -0.687    +0.861 j
p[ 4] =    -0.687    -0.861 j
p[ 5] =    -0.992    +0.478 j
p[ 6] =    -0.992    -0.478 j  
\endverbatim
\n

В каталоге `dat` будет создан файл `butter_ap_zp.txt`. \n

Пакет GNUPLOT произведет построение карты полюсов по
сохранненным в  `dat/butter_ap_zp.txt` данным:

\image html butter_ap_zp_test.png

\author
Бахурин Сергей
www.dsplib.org 
***************************************************************************** */
#endif
int DSPL_API butter_ap_zp(int ord, double rp, complex_t* z, int* nz,
                          complex_t *p, int* np)
{
    double alpha;
    double theta;
    double ep;
    int r;
    int L;
    int ind = 0, k;

    if(rp < 0 || rp == 0)
        return ERROR_FILTER_RP;
    if(ord < 1)
        return ERROR_FILTER_ORD;
    if(!z || !p || !nz || !np)
        return ERROR_PTR;

    ep = sqrt(pow(10.0, rp*0.1) - 1.0);
    r = ord % 2;
    L = (int)((ord-r)/2);

    alpha = pow(ep, -1.0/(double)ord);
    if(r)
    {
        RE(p[ind]) = -alpha;
        IM(p[ind]) = 0.0;
        ind++;
    }
    for(k = 0; k < L; k++)
    {
        theta = M_PI*(double)(2*k + 1)/(double)(2*ord);
        RE(p[ind]) = RE(p[ind+1]) = -alpha *    sin(theta);
        IM(p[ind]) =     alpha *    cos(theta);
        IM(p[ind+1]) =    -alpha *    cos(theta);
        ind+=2;
    }
    *np = ord;
    *nz = 0;
    return RES_OK;
}


