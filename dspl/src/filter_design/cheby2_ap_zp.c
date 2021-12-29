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
\fn int cheby2_ap_zp(int ord, double rs, complex_t* z, int* nz, 
                     complex_t* p, int* np)

\brief
Function calculates arrays of zeros and poles for analog normlized lowpass 
Chebyshev type 2 filter transfer function \f$ H(s) \f$ order `ord` .

Analog normalized Chebyshev type 2 filter lowpass filter has \f$Rs\f$ dB 
suppression in stopband. 
Also analog normalized Chebyshev type 2 filter magnitude equals \f$-Rs\f$ dB 
for angular frequency \f$\omega = 1\f$ rad/s.


\param[in]  ord
Filter order. \n
Number of zeros and poles of filter can be less or equal `ord`. \n
\n

\param[in]  rs
Suppression level in stopband (dB). \n
This parameter sets filter supression for \f$\omega \geq 1\f$ rad/s frequency. \n
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

Example of normalized Chebyshev type 2 lowpass filter
 zeros and poles calculation:
\include cheby2_ap_zp_test.c

Result:

\verbatim
Chebyshev type 2 filter zeros: 6
z[ 0] =     0.000    +1.026 j
z[ 1] =     0.000    -1.026 j
z[ 2] =     0.000    +1.279 j
z[ 3] =     0.000    -1.279 j
z[ 4] =     0.000    +2.305 j
z[ 5] =     0.000    -2.305 j
Chebyshev type 2 filter poles: 7
p[ 0] =    -1.203    +0.000 j
p[ 1] =    -0.113    +0.772 j
p[ 2] =    -0.113    -0.772 j
p[ 3] =    -0.398    +0.781 j
p[ 4] =    -0.398    -0.781 j
p[ 5] =    -0.852    +0.642 j
p[ 6] =    -0.852    -0.642 j 
\endverbatim
\n

In `dat` folder will be created `cheby2_ap_z.txt` and 
`cheby2_ap_z.txt` files which keeps zeros and poles vectors. \n

In addition, GNUPLOT will build the following graphs 
from data stored in the files:

\image html cheby2_ap_zp_test.png


\author Sergey Bakhurin www.dsplib.org
***************************************************************************** */
#endif
#ifdef DOXYGEN_RUSSIAN
/*! ****************************************************************************
\ingroup IIR_FILTER_DESIGN_GROUP
\fn int cheby2_ap_zp(int ord, double rs, complex_t* z, int* nz, 
                     complex_t* p, int* np)

\brief
Расчет массивов нулей и полюсов передаточной функции \f$ H(s) \f$ 
аналогового  нормированного ФНЧ Чебышёва второго рода.

Функция рассчитывает значения нулей и полюсов передаточной функции 
\f$H(s)\f$ аналогового нормированного ФНЧ Чебышёва второго рода порядка `ord` с 
частотой   заграждения 1 рад/с по уровню \f$-R_s\f$ дБ. \n 


\param[in]  ord
Порядок фильтра. \n
\n

\param[in]  rs
Уровень подавления АЧХ в полосе загражения (дБ). \n
Параметр задает уровень подавления сигнала в полосе частот от 1 рад/с и выше. \n
Значение должно быть положительным. \n 
\n

\param[out]  z
Указатель на массив комплексных нулей передаточной функции \f$H(s)\f$. \n
Максимальный размер вектора  вектора `[ord x 1]`. \n
Память должна быть выделена. \n 
\n

\param[out]  nz
Указатель на переменную количества нулей передаточной функции \f$H(s)\f$. \n
По данному указателю будет записано количество нулей фильтра, которые были 
рассчитаны и помещены в вектор `z`. \n
Память должна быть выделена. \n 
\n

\param[out]  p
Указатель на массив комплексных полюсов передаточной функции \f$H(s)\f$. \n
Максимальный размер вектора  вектора `[ord x 1]`. \n
Память должна быть выделена. \n 
\n

\param[out]  np
Указатель на переменную количества полюсов передаточной функции \f$H(s)\f$. \n
По данному указателю будет записано количество нулей 
фильтра, которые были 
рассчитаны и помещены в вектор `p`. \n
Память должна быть выделена. \n 
\n

\return
`RES_OK` --- массивы нулей и полюсов рассчитаны успешно. \n
В противном случае \ref ERROR_CODE_GROUP "код ошибки". \n

Пример использования функции `cheby2_ap_zp`:

Пример программы рассчета нулей и полюсов нормированного 
ФНЧ Чебышева первого рода:
\include cheby2_ap_zp_test.c

Результат выполнения программы:

\verbatim
Chebyshev type 2 filter zeros: 6
z[ 0] =     0.000    +1.026 j
z[ 1] =     0.000    -1.026 j
z[ 2] =     0.000    +1.279 j
z[ 3] =     0.000    -1.279 j
z[ 4] =     0.000    +2.305 j
z[ 5] =     0.000    -2.305 j
Chebyshev type 2 filter poles: 7
p[ 0] =    -1.203    +0.000 j
p[ 1] =    -0.113    +0.772 j
p[ 2] =    -0.113    -0.772 j
p[ 3] =    -0.398    +0.781 j
p[ 4] =    -0.398    -0.781 j
p[ 5] =    -0.852    +0.642 j
p[ 6] =    -0.852    -0.642 j 
\endverbatim
\n

В каталоге `dat` будет создан файлы `cheby2_ap_z.txt` и `cheby2_ap_z.txt`,
хранящие наборы нулей и полюсов на комплексной плоскости. \n

Пакет GNUPLOT произведет построение карты полюсов по
сохранненным в `dat/cheby2_ap_z.txt` и `dat/cheby2_ap_p.txt` данным:

\image html cheby2_ap_zp_test.png


\author
Бахурин Сергей
www.dsplib.org 
***************************************************************************** */
#endif
int DSPL_API cheby2_ap_zp(int ord, double rs, complex_t* z, int* nz,
                          complex_t *p, int* np)
{
    double es;
    int L, r, k;
    double beta;
    int iz, ip;

    double alpha;
    double chb, shb, sa, ca;
    double ssh2, cch2;

    if(rs < 0 || rs == 0)
        return ERROR_FILTER_RS;
    if(ord < 1)
        return ERROR_FILTER_ORD;
    if(!z || !p || !nz || !np)
        return ERROR_PTR;

    es = sqrt(pow(10.0, rs*0.1) - 1.0);
    r = ord % 2;
    L = (int)((ord-r)/2);

    beta     = asinh(es)/(double)ord;

    chb = cosh(beta);
    shb = sinh(beta);

    iz = ip = 0;

    if(r)
    {
        RE(p[0]) = -1.0 / sinh(beta);
        IM(p[0]) =    0.0;
        ip = 1;
    }

    for(k = 0; k < L; k++)
    {
        alpha = M_PI*(double)(2*k + 1)/(double)(2*ord);
        sa = sin(alpha);
        ca = cos(alpha);
        ssh2    = sa*shb;
        ssh2 *= ssh2;

        cch2    = ca*chb;
        cch2 *= cch2;

        RE(z[iz]) = RE(z[iz+1]) = 0.0;
        IM(z[iz]) = 1.0 / ca;
        IM(z[iz+1]) = -IM(z[iz]);
        iz+=2;

        RE(p[ip]) = RE(p[ip+1]) = -sa*shb / (ssh2 + cch2);
        IM(p[ip]) = ca*chb / (ssh2 + cch2);
        IM(p[ip+1]) = -IM(p[ip]);
        ip+=2;
    }
    *nz = iz;
    *np = ip;

    return RES_OK;
}

