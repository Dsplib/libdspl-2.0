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
\fn int ellip_ap_zp(int ord, double rp, double rs, complex_t* z, int* nz, 
                     complex_t* p, int* np)

\brief
Function calculates arrays of zeros and poles for analog normlized lowpass 
elliptic filter transfer function \f$ H(s) \f$ order `ord` .

\param[in]  ord
Filter order. \n
Number of zeros and poles of filter can be less or equal `ord`. \n
\n

\param[in]  rp
Magnitude ripple in passband (dB). \n
This parameter sets maximum filter distortion from 0 to 1 rad/s frequency. \n
Parameter must be positive. \n 
\n

\param[in]  rs
Suppression level in stopband (dB). \n
This parameter sets filter suppression
for \f$\omega \geq 1\f$ rad/s frequency. \n
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

Example of normalized elliptic lowpass filter zeros and poles calculation:
\include ellip_ap_zp_test.c

Result:

\verbatim
Elliptic filter zeros: 6
z[ 0] =     0.000    +1.053 j
z[ 1] =     0.000    -1.053 j
z[ 2] =     0.000    +1.136 j
z[ 3] =     0.000    -1.136 j
z[ 4] =     0.000    +1.626 j
z[ 5] =     0.000    -1.626 j
Elliptic filter poles: 7
p[ 0] =    -0.358    +0.000 j
p[ 1] =    -0.011    +1.000 j
p[ 2] =    -0.011    -1.000 j
p[ 3] =    -0.060    +0.940 j
p[ 4] =    -0.060    -0.940 j
p[ 5] =    -0.206    +0.689 j
p[ 6] =    -0.206    -0.689 j
\endverbatim
\n

In `dat` folder will be created `ellip_ap_z.txt` and 
`ellip_ap_z.txt` files which keeps zeros and poles vectors. \n

In addition, GNUPLOT will build the following graphs 
from data stored in the files:

\image html ellip_ap_zp_test.png


\author Sergey Bakhurin www.dsplib.org
***************************************************************************** */
#endif
#ifdef DOXYGEN_RUSSIAN
/*! ****************************************************************************
\ingroup IIR_FILTER_DESIGN_GROUP
\fn int ellip_ap_zp(int ord, double rp, double rs, complex_t* z, int* nz, 
                     complex_t* p, int* np)

\brief
Расчет массивов нулей и полюсов передаточной функции \f$ H(s) \f$ 
аналогового нормированного эллиптического ФНЧ.

\param[in]  ord
Порядок фильтра. \n
\n


\param[in]  rp
Неравномерность АЧХ в полосе пропускания (дБ). \n
Параметр задает уровень искажений в полосе от 0 до 1 рад/с. \n
Значение должно быть положительным. \n
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
эллиптического ФНЧ :
\include ellip_ap_zp_test.c

Результат выполнения программы:

\verbatim
Elliptic filter zeros: 6
z[ 0] =     0.000    +1.053 j
z[ 1] =     0.000    -1.053 j
z[ 2] =     0.000    +1.136 j
z[ 3] =     0.000    -1.136 j
z[ 4] =     0.000    +1.626 j
z[ 5] =     0.000    -1.626 j
Elliptic filter poles: 7
p[ 0] =    -0.358    +0.000 j
p[ 1] =    -0.011    +1.000 j
p[ 2] =    -0.011    -1.000 j
p[ 3] =    -0.060    +0.940 j
p[ 4] =    -0.060    -0.940 j
p[ 5] =    -0.206    +0.689 j
p[ 6] =    -0.206    -0.689 j
\endverbatim
\n

В каталоге `dat` будет создан файлы `ellip_ap_z.txt` и `ellip_ap_z.txt`,
хранящие наборы нулей и полюсов на комплексной плоскости. \n

Пакет GNUPLOT произведет построение карты полюсов по
сохранненным в `dat/ellip_ap_z.txt` и `dat/ellip_ap_p.txt` данным:

\image html ellip_ap_zp_test.png


\author
Бахурин Сергей
www.dsplib.org 
***************************************************************************** */
#endif
int DSPL_API ellip_ap_zp(int ord, double rp, double rs,
                         complex_t* z, int* nz, complex_t* p, int* np)
{
    double es, ep;
    int L, r, n, res;
    int iz, ip;
    double ke, k, u, t;
    complex_t tc, v0, jv0;


    if(rp < 0 || rp == 0)
        return ERROR_FILTER_RP;
    if(rs < 0 || rs == 0)
        return ERROR_FILTER_RS;
    if(ord < 1)
        return ERROR_FILTER_ORD;
    if(!z || !p || !nz || !np)
        return ERROR_PTR;

    es = sqrt(pow(10.0, rs*0.1) - 1.0);
    ep = sqrt(pow(10.0, rp*0.1) - 1.0);
    ke = ep / es;

    r = ord % 2;
    L = (int)((ord-r)/2);

    res = ellip_modulareq(rp, rs, ord, &k);
    if(res != RES_OK)
        return res;
    // v0
    RE(tc) = 0.0;
    IM(tc) = 1.0 / ep;

    ellip_asn_cmplx(&tc, 1, ke, &v0);

    t = RE(v0);
    RE(v0) = IM(v0) / (double)ord;
    IM(v0) = -t / (double)ord;

    RE(jv0) = -IM(v0);
    IM(jv0) =    RE(v0);

    iz = ip = 0;

    if(r)
    {
        res = ellip_sn_cmplx(&jv0, 1, k, &tc);
        if(res != RES_OK)
            return res;
        RE(p[0]) = -IM(tc);
        IM(p[0]) =    RE(tc);
        ip = 1;
    }

    for(n = 0; n < L; n++)
    {
        u = (double)(2 * n + 1)/(double)ord;

        res = ellip_cd(& u, 1, k, &t);
        if(res != RES_OK)
            return res;

        RE(z[iz]) = RE(z[iz+1]) = 0.0;
        IM(z[iz])     =    1.0/(k*t);
        IM(z[iz+1]) = -1.0/(k*t);
        iz+=2;

        RE(tc) = u - RE(jv0);
        IM(tc) =     - IM(jv0);

        res = ellip_cd_cmplx(&tc, 1, k, p+ip+1);
        if(res != RES_OK)
            return res;

        RE(p[ip]) = -IM(p[ip+1]);
        IM(p[ip]) =    RE(p[ip+1]);

        RE(p[ip+1]) =     RE(p[ip]);
        IM(p[ip+1]) =    -IM(p[ip]);

        ip+=2;
    }
    *nz = iz;
    *np = ip;

    return RES_OK;
}

