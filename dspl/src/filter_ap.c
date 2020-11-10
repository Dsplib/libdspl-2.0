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






#ifdef DOXYGEN_ENGLISH
/*! ****************************************************************************
\ingroup IIR_FILTER_DESIGN_GROUP
\fn int cheby1_ap(double Rp, int ord, double* b, double* a)

\brief
Function calculates the transfer function \f$ H(s) \f$ coefficients of
analog normalized lowpass Chebyshev type 1 filter.

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

\include cheby1_ap_test.c

Result:

\verbatim
b[ 0] =     0.125     a[ 0] =     0.177
b[ 1] =     0.000     a[ 1] =     0.405
b[ 2] =     0.000     a[ 2] =     1.169
b[ 3] =     0.000     a[ 3] =     0.582
b[ 4] =     0.000     a[ 4] =     1.000
\endverbatim
\n

In `dat` folder will be created 3 files: \n

\verbatim
cheby1_ap_test_mag.txt    magnitude
cheby1_ap_test_phi.txt    phase response
cheby1_ap_test_tau.txt    group delay
\endverbatim

In addition, GNUPLOT will build the following graphs from data stored in files:

\image html cheby1_ap_test.png

\author Sergey Bakhurin www.dsplib.org
***************************************************************************** */
#endif
#ifdef DOXYGEN_RUSSIAN
/*! ****************************************************************************
\ingroup IIR_FILTER_DESIGN_GROUP
\fn int cheby1_ap(double Rp, int ord, double* b, double* a)

\brief
Расчет передаточной характеристики \f$ H(s) \f$ аналогового 
нормированного ФНЧ Чебышёва первого рода.

Функция рассчитывает коэффициенты передаточной характеристики 
\f$ H(s)\f$ аналогового нормированного ФНЧ Чебышёва первого рода 
порядка `ord` с частотой среза 1 рад/с по уровню \f$-R_p\f$ дБ. \n 

Особенностью фильтра Чебышёва первого рода являются 
равноволновые пульсации АЧХ в полосе пропускания.

\param[in]  Rp
Неравномерность АЧХ в полосе пропускания (дБ). \n
Параметр задает уровень искажений в полосе от 0 до 1 рад/с. \n
Значение должно быть положительным. \n
\n

\param[in]  ord
Порядок фильтра. \n
Количество коэффициентов числителя и знаменателя
передаточной функции \f$ H(s)\f$ равно `ord+1`. \n
\n

\param[out]  b
Указатель на вектор коэффициентов числителя передаточной функции \f$H(s)\f$. \n
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

Пример использования функции `cheby1_ap`:

\include cheby1_ap_test.c

Результат работы программы:

\verbatim
b[ 0] =     0.125     a[ 0] =     0.177
b[ 1] =     0.000     a[ 1] =     0.405
b[ 2] =     0.000     a[ 2] =     1.169
b[ 3] =     0.000     a[ 3] =     0.582
b[ 4] =     0.000     a[ 4] =     1.000
\endverbatim
\n

В каталоге `dat` будут созданы три файла: \n

\verbatim
cheby1_ap_test_mag.txt    АЧХ фильтра
cheby1_ap_test_phi.txt    ФЧХ фильтра
cheby1_ap_test_tau.txt    ГВЗ фильтра
\endverbatim
\n

Кроме того программа GNUPLOT произведет построение следующих графиков 
по сохраненным в файлах данным:

\image html cheby1_ap_test.png

\author
Бахурин Сергей
www.dsplib.org
***************************************************************************** */
#endif
int DSPL_API cheby1_ap(double rp, int ord, double* b, double* a)
{
    int res;
    complex_t *z = NULL;
    complex_t *p = NULL;
    int nz, np, k;
    complex_t h0 = {1.0, 0.0};
    double tmp;


    if(rp < 0.0)
        return ERROR_FILTER_RP;
    if(ord < 1)
        return ERROR_FILTER_ORD;
    if(!a || !b)
        return ERROR_PTR;

    z = (complex_t*) malloc(ord*sizeof(complex_t));
    p = (complex_t*) malloc(ord*sizeof(complex_t));


    res = cheby1_ap_zp(ord, rp, z, &nz, p, &np);
    if(res != RES_OK)
        goto exit_label;

    res = filter_zp2ab(z, nz, p, np, ord, b, a);
    if(res != RES_OK)
        goto exit_label;


    if(!(ord % 2))
        RE(h0) = 1.0 / pow(10.0, rp*0.05);

    for(k = 0; k < np; k++)
    {
        tmp    = CMRE(h0, p[k]);
        IM(h0) = CMIM(h0, p[k]);
        RE(h0) = tmp;
    }

    b[0] = fabs(RE(h0));

exit_label:
    if(z)
        free(z);
    if(p)
        free(p);
    return res;
}





#ifdef DOXYGEN_ENGLISH
/*! ****************************************************************************
\ingroup IIR_FILTER_DESIGN_GROUP
\fn int cheby1_ap_zp( int ord, double rp, complex_t* z, int* nz, 
                      complex_t* p, int* np)
\brief
Function calculates arrays of zeros and poles for analog normlized lowpass 
Chebyshev type 1 filter transfer function \f$ H(s) \f$ order `ord` .

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
Normalized Chebyshev type 1 lowpass filter has no finite zeros. 
So `z` vector will not changed and in pointer `nz` will write 0 value. \n

Example of normalized Chebyshev type 1 lowpass filter
 zeros and poles calculation:
\include cheby1_ap_zp_test.c

Result:

\verbatim
Chebyshev type 1 filter zeros: 0
Chebyshev type 1 filter poles: 7
p[ 0] =    -0.256    +0.000 j
p[ 1] =    -0.057    +1.006 j
p[ 2] =    -0.057    -1.006 j
p[ 3] =    -0.160    +0.807 j
p[ 4] =    -0.160    -0.807 j
p[ 5] =    -0.231    +0.448 j
p[ 6] =    -0.231    -0.448 j 
\endverbatim
\n

In `dat` folder will be created `cheby1_ap_zp.txt` file. \n

In addition, GNUPLOT will build the following graphs 
from data stored in `dat/cheby1_ap_zp.txt` file:

\image html cheby1_ap_zp_test.png


\author Sergey Bakhurin www.dsplib.org
***************************************************************************** */
#endif
#ifdef DOXYGEN_RUSSIAN
/*! ****************************************************************************
\ingroup IIR_FILTER_DESIGN_GROUP
\fn int cheby1_ap_zp(int ord, double rp, complex_t* z, int* nz, complex_t* p, int* np)
\brief
Расчет массивов нулей и полюсов передаточной функции   \f$ H(s) \f$ 
аналогового нормированного ФНЧ Чебышёва первого рода.

Функция рассчитывает значения нулей и полюсов передаточной функции
\f$ H(s)\f$  аналогового нормированного ФНЧ Чебышёва первого рода 
порядка `ord` с частотой среза 1 рад/с по уровню \f$-R_p\f$ дБ, с 
неравномерностью в полосе пропускания \f$ R_p \f$ дБ. \n 

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
Максимальный размер вектора `[ord x 1]`. \n
Память должна быть выделена. \n
\n

\param[out]  nz
Указатель на переменную количества нулей 
передаточной функции \f$H(s)\f$. \n
По данному указателю будет записано количество нулей фильтра, 
которые были рассчитаны и помещены в вектор `z`. \n
Память должна быть выделена. \n
\n

\param[out]  p
Указатель на массив комплексных полюсов 
передаточной характеристики \f$H(s)\f$. \n
Максимальный размер вектора  вектора `[ord x 1]`. \n
Память должна быть выделена. \n
\n

\param[out]  np
Указатель на переменную количества полюсов передаточной функции \f$ H(s)\f$. \n
По данному указателю будет записано количество нулей фильтра, которые были 
рассчитаны и помещены в вектор `p`. \n
Память должна быть выделена. \n
\n

\return
`RES_OK` --- массивы нулей и полюсов рассчитаны успешно. \n
В противном случае \ref ERROR_CODE_GROUP "код ошибки". \n

\note
Нормированный ФНЧ Чебышёва первого рода не имеет нулей, поэтому массив 
нулей `z` не будет изменен, а по указателю `nz`  будет записан 0. \n


Пример программы рассчета нулей и полюсов нормированного 
ФНЧ Чебышева первого рода:
\include cheby1_ap_zp_test.c

Результат выполнения программы:

\verbatim
Chebyshev type 1 filter zeros: 0
Chebyshev type 1 filter poles: 7
p[ 0] =    -0.256    +0.000 j
p[ 1] =    -0.057    +1.006 j
p[ 2] =    -0.057    -1.006 j
p[ 3] =    -0.160    +0.807 j
p[ 4] =    -0.160    -0.807 j
p[ 5] =    -0.231    +0.448 j
p[ 6] =    -0.231    -0.448 j
\endverbatim
\n

В каталоге `dat` будет создан файл `cheby1_ap_zp.txt`. \n

Пакет GNUPLOT произведет построение карты полюсов по
сохранненным в  `dat/cheby1_ap_zp.txt` данным:

\image html cheby1_ap_zp_test.png

\author Бахурин Сергей www.dsplib.org 
***************************************************************************** */
#endif
int DSPL_API cheby1_ap_zp(int ord, double rp, complex_t* z, int* nz,
                          complex_t* p, int* np)
{
    double theta;
    double ep;
    double beta;
    double shbeta;
    double chbeta;
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


    beta     = asinh(1.0/ep)/(double)ord;
    chbeta = cosh(beta);
    shbeta = sinh(beta);

    if(r)
    {
        RE(p[ind]) = -shbeta;
        IM(p[ind]) =    0.0;
        ind++;
    }
    for(k = 0; k < L; k++)
    {
        theta = M_PI*(double)(2*k + 1)/(double)(2*ord);
        RE(p[ind]) = RE(p[ind+1]) = -shbeta *    sin(theta);
        IM(p[ind]) =    chbeta *    cos(theta);
        IM(p[ind+1]) = -IM(p[ind]);
        ind+=2;
    }
    *np = ord;
    *nz = 0;
    return RES_OK;
}



#ifdef DOXYGEN_ENGLISH
/*! ****************************************************************************
\ingroup IIR_FILTER_DESIGN_GROUP
\fn int cheby2_ap(double Rs, int ord, double *b, double *a)

\brief
Function calculates the transfer function \f$ H(s) \f$ coefficients of
analog normalized lowpass Chebyshev type 2 filter.

Analog normalized Chebyshev type 2 filter lowpass filter has \f$Rs\f$ dB 
suppression in stopband. 
Also analog normalized Chebyshev type 2 filter magnitude equals \f$-Rs\f$ dB 
for angular frequency \f$\omega = 1\f$ rad/s.

\param[in]  Rs
Suppression level in stopband (dB). \n
This parameter sets filter supression for \f$\omega \geq 1\f$ rad/s frequency. \n
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

\include cheby2_ap_test.c

Result:

\verbatim
b[ 0] =     0.008     a[ 0] =     0.008
b[ 1] =     0.000     a[ 1] =     0.068
b[ 2] =     0.008     a[ 2] =     0.300
b[ 3] =     0.000     a[ 3] =     0.774
b[ 4] =     0.001     a[ 4] =     1.000 
\endverbatim
\n

In `dat` folder will be created 3 files: \n

\verbatim
cheby2_ap_test_mag.txt    magnitude
cheby2_ap_test_phi.txt    phase response
cheby2_ap_test_tau.txt    group delay
\endverbatim

In addition, GNUPLOT will build the following graphs from data stored in files:

\image html cheby2_ap_test.png

\author Sergey Bakhurin www.dsplib.org
***************************************************************************** */
#endif
#ifdef DOXYGEN_RUSSIAN
/*! ****************************************************************************
\ingroup IIR_FILTER_DESIGN_GROUP
\fn int cheby2_ap(double Rs, int ord, double *b, double *a)

\brief
Расчет передаточной характеристики \f$ H(s) \f$ аналогового 
нормированного ФНЧ Чебышёва второго рода.

Функция рассчитывает коэффициенты передаточной характеристики \f$H(s)\f$
аналогового нормированного ФНЧ Чебышёва второго рода порядка `ord` 
с частотой заграждения 1 рад/с по уровню \f$-R_s\f$ дБ. \n

Особенностью фильтра Чебышёва второго рода являются:  \n
1) равноволновые пульсации  АЧХ в полосе заграждения. \n
2) уровень АЧХ \f$H(j\cdot 1) = -R_s\f$ дБ. \n

\param[in]  Rs
Уровень подавления в полосе пропускания (дБ). \n
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

Пример использования функции `cheby1_ap`:

\include cheby2_ap_test.c

Результат работы программы:

\verbatim
b[ 0] =     0.008     a[ 0] =     0.008
b[ 1] =     0.000     a[ 1] =     0.068
b[ 2] =     0.008     a[ 2] =     0.300
b[ 3] =     0.000     a[ 3] =     0.774
b[ 4] =     0.001     a[ 4] =     1.000 
\endverbatim
\n

В каталоге `dat` будут созданы три файла: \n

\verbatim
cheby2_ap_test_mag.txt    АЧХ фильтра   
cheby2_ap_test_phi.txt    ФЧХ фильтра
cheby2_ap_test_tau.txt    ГВЗ фильтра
\endverbatim
\n

Кроме того программа GNUPLOT произведет построение следующих графиков 
по сохраненным в файлах данным:

\image html cheby2_ap_test.png


\return
`RES_OK` --- фильтр рассчитан успешно. \n
В противном случае \ref ERROR_CODE_GROUP "код ошибки". \n

\author
Бахурин Сергей
www.dsplib.org
***************************************************************************** */
#endif
int DSPL_API cheby2_ap(double rs, int ord, double* b, double* a)
{
    int res;
    complex_t *z = NULL;
    complex_t *p = NULL;
    int nz, np;
    double norm;


    if(rs < 0.0)
        return ERROR_FILTER_RP;
    if(ord < 1)
        return ERROR_FILTER_ORD;
    if(!a || !b)
        return ERROR_PTR;

    z = (complex_t*) malloc(ord*sizeof(complex_t));
    p = (complex_t*) malloc(ord*sizeof(complex_t));


    res = cheby2_ap_zp(ord, rs, z, &nz, p, &np);
    if(res != RES_OK)
        goto exit_label;

    res = filter_zp2ab(z, nz, p, np, ord, b, a);
    if(res != RES_OK)
        goto exit_label;

    norm = a[0] / b[0];

    for(nz = 0; nz < ord+1; nz++)
        b[nz]*=norm;

exit_label:
    if(z)
        free(z);
    if(p)
        free(p);
    return res;
}



#ifdef DOXYGEN_ENGLISH

#endif
#ifdef DOXYGEN_RUSSIAN

#endif
int DSPL_API cheby2_ap_wp1(double rp, double rs, int ord, double* b, double* a)
{
    int err;
    double es, gp, alpha, beta, y, wp;

    if(rp <= 0)
        return    ERROR_FILTER_RP;

    err = cheby2_ap(rs, ord, b, a);
    if(err!=RES_OK)
        goto exit_label;

    es = sqrt(pow(10.0, rs*0.1) - 1.0);
    gp = pow(10.0, -rp*0.05);
    alpha = gp * es / sqrt(1.0 - gp*gp);
    beta = alpha + sqrt(alpha * alpha - 1.0);
    y = log(beta)/ (double)ord;
    wp = 2.0 / (exp(y) + exp(-y));
    
    err = low2low(b, a, ord, wp, 1.0, b, a);

exit_label:
    return err;
}




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







#ifdef DOXYGEN_ENGLISH
/*! ****************************************************************************
\ingroup IIR_FILTER_DESIGN_GROUP
\fn int ellip_ap(double rp, double rs, int ord, double* b, double* a)

\brief
Function calculates the transfer function \f$ H(s) \f$ coefficients of
analog normalized lowpass elliptic filter order `ord` with passband ripple
`rp` dB and stopband suppression equals `rs` dB.

\param[in]  rp
Magnitude ripple in passband (dB). \n
This parameter sets maximum filter distortion from 0 to 1 rad/s frequency. \n
Parameter must be positive. \n 
\n


\param[in]  rs
Suppression level in stopband (dB). \n
This parameter sets filter supression for \f$\omega \geq 1\f$ rad/s frequency. \n
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

\include ellip_ap_test.c

Result:

\verbatim
b[ 0] =     0.268     a[ 0] =     0.301
b[ 1] =     0.000     a[ 1] =     0.764
b[ 2] =     0.045     a[ 2] =     1.472
b[ 3] =     0.000     a[ 3] =     0.948
b[ 4] =     0.001     a[ 4] =     1.000
\endverbatim
\n

In `dat` folder will be created 3 files: \n

\verbatim
ellip_ap_test_mag.txt    magnitude
ellip_ap_test_phi.txt    phase response
ellip_ap_test_tau.txt    group delay
\endverbatim

In addition, GNUPLOT will build the following graphs from data stored in files:

\image html ellip_ap_test.png

\author Sergey Bakhurin www.dsplib.org
***************************************************************************** */
#endif
#ifdef DOXYGEN_RUSSIAN
/*! ****************************************************************************
\ingroup IIR_FILTER_DESIGN_GROUP
\fn int ellip_ap(double rp, double rs, int ord, double* b, double* a)

\brief
Расчет передаточной характеристики \f$ H(s) \f$ аналогового 
нормированного эллиптического ФНЧ.

Функция рассчитывает коэффициенты передаточной характеристики \f$H(s)\f$
аналогового нормированного эллиптического ФНЧ порядка `ord` 
с частотой среза 1 рад/с по уровню \f$-R_p\f$ дБ. \n

Особенностью эллиптического фильтра являются равноволновые пульсации 
АЧХ как в полосе пропускания, так и в полосе заграждения, в результате
чего обеспечиваеся минимальная переходная полоса фильтра. \n

\param[in]  rp
Уровень пульсаций в полосе пропускания (дБ). \n
Значение должно быть положительным. \n 
\n

\param[in]  rs
Уровень подавления в полосе заграждения (дБ). \n
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

Пример использования функции `ellip_ap`:

\include ellip_ap_test.c

Результат работы программы:

\verbatim
b[ 0] =     0.268     a[ 0] =     0.301
b[ 1] =     0.000     a[ 1] =     0.764
b[ 2] =     0.045     a[ 2] =     1.472
b[ 3] =     0.000     a[ 3] =     0.948
b[ 4] =     0.001     a[ 4] =     1.000
\endverbatim
\n

В каталоге `dat` будут созданы три файла: \n

\verbatim
ellip_ap_test_mag.txt    АЧХ фильтра   
ellip_ap_test_phi.txt    ФЧХ фильтра
ellip_ap_test_tau.txt    ГВЗ фильтра
\endverbatim
\n

Кроме того программа GNUPLOT произведет построение следующих графиков 
по сохраненным в файлах данным:

\image html ellip_ap_test.png


\return
`RES_OK` --- фильтр рассчитан успешно. \n
В противном случае \ref ERROR_CODE_GROUP "код ошибки". \n

\author
Бахурин Сергей
www.dsplib.org
***************************************************************************** */
#endif
int DSPL_API ellip_ap(double rp, double rs, int ord, double* b, double* a)
{
    int res;
    complex_t *z = NULL;
    complex_t *p = NULL;
    int nz, np;
    double norm, g0;


    if(rp < 0.0)
        return ERROR_FILTER_RP;
    if(rs < 0.0)
        return ERROR_FILTER_RS;
    if(ord < 1)
        return ERROR_FILTER_ORD;
    if(!a || !b)
        return ERROR_PTR;

    z = (complex_t*) malloc(ord*sizeof(complex_t));
    p = (complex_t*) malloc(ord*sizeof(complex_t));


    res = ellip_ap_zp(ord, rp, rs, z, &nz, p, &np);
    if(res != RES_OK)
        goto exit_label;

    res = filter_zp2ab(z, nz, p, np, ord, b, a);
    if(res != RES_OK)
        goto exit_label;


    g0 = 1.0;
    if(!(ord % 2))
    {
        g0 = 1.0 / pow(10.0, rp*0.05);
    }


    norm = g0 * a[0] / b[0];

    for(nz = 0; nz < ord+1; nz++)
        b[nz]*=norm;

    exit_label:
    if(z)
        free(z);
    if(p)
        free(p);
    return res;
}






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




#ifdef DOXYGEN_ENGLISH
/*! ****************************************************************************
\ingroup IIR_FILTER_DESIGN_GROUP
\fn int filter_zp2ab(complex_t *z, int nz, complex_t *p, int np, int ord, 
 double* b, double* a)  
\brief  
Function recalculates complex zeros and poles of transfer function \f$ H(s) \f$
to the coefficients of \f$ H(s) \f$ numerator and denominator polynomials. 

Transfer function can we described as:
\f[
H(s) = 
\frac{\sum\limits_{n = 0}^{N_z} b_n  s^n}{\sum\limits_{m = 0}^{N_p} a_m  s^m} = 
\frac{\prod\limits_{n = 0}^{N_z}(s-z_n)}{\prod\limits_{m = 0}^{N_p} (s-p_m)}
\f]

\param[in]  z
Pointer to the vector of transfer function zeros. \n
Vector size is `[nz x 1]`. \n
Pointer can be `NULL` if filter has no finite zeros (`nz=0`). \n 
\n
    
\param[in]  nz
Number of fitite zeros (can be zero). \n
\n 

\param[in]  p
Pointer to the vector of transfer function poles. \n
Vector size is `[np x 1]`. \n
This pointer cannot be `NULL`. \n
\n

\param[in]  np
Size of vector of transfer function poles (`p` vector size). \n 
\n

\param[in]  ord
Filter order. \n
Number of \f$H(s)\f$ numerator and denominator coefficients equals `ord+1`. \n 
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

\note
Function calculates real `b` and  `a` coefficients of \f$H(s)\f$.
It means that zeros and poles vectors must have real values or conjugate pairs 
to get zeros image part of `b` and  `a` coefficients. This function ignores 
image part of `b` and  `a` coeeffitients if the requirements for zeros 
and poles are not fulfilled.

\author Sergey Bakhurin www.dsplib.org  
***************************************************************************** */
#endif
#ifdef DOXYGEN_RUSSIAN
/*! ****************************************************************************
\ingroup IIR_FILTER_DESIGN_GROUP
\fn int filter_zp2ab(complex_t *z, int nz, complex_t *p, int np, int ord, 
 double* b, double* a)  
\brief  Функция пересчета нулей и полюсов аналогового фильтра в коэффициенты 
        передаточной характеристики \f$ H(s) \f$
  
\f[
H(s) = 
\frac{\sum_{n = 0}^{N_z} b_n \cdot s^n}{\sum_{m = 0}^{N_p} a_m \cdot s^m} = 
\frac{\prod_{n = 0}^{N_z}(s-z_n)}{\prod_{m = 0}^{N_p} (s-p_m)}
\f]

\param[in]  z
Указатель на массив нулей передаточной характеристики. \n
Размер вектора `[nz x 1]`. \n
Указатель может быть `NULL` если фильтр не имеет конечных нулей (`nz=0`). \n 
\n
    
\param[in]  nz
Размер вектора нулей передаточной характеристики (может быть равен 0). \n
\n 

\param[in]  p
Указатель на массив полюсов передаточной характеристики. \n
Размер вектора `[np x 1]`. \n
Указатель не может быть `NULL`. \n
Память должна быть выделена. \n 
\n

\param[in]  np
Размер вектора полюсов передаточной характеристики (не может быть равен 0). \n 
\n

\param[in]  ord
Порядок фильтра для которого рассчитаны нули и полюса. \n
Количество коэффициентов числителя и знаменателя 
передаточной функции \f$H(s)\f$ равно `ord+1`. \n \n

\param[out]  b
Указатель на вектор коэффициентов числителя передаточной функции \f$H(s)\f$. \n
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
`RES_OK` --- пересчет произведен успешно. \n
В противном случае \ref ERROR_CODE_GROUP "код ошибки". \n

\note
Функция возвращает вещественные значения коэффициентов `b` и `a` 
передаточной функции. Это означает, что вектора нулей и полюсов 
должны хранить вещественные значения или комплексно-сопряженные пары 
нулей и полюсов, потому что мнимая часть коэффициентов `b` и `a`
игнорируется и не сохраняется.

\author
Бахурин Сергей
www.dsplib.org  
***************************************************************************** */
#endif
int DSPL_API filter_zp2ab(complex_t* z, int nz, complex_t* p, int np,
                          int ord, double* b, double* a)
{
    complex_t *acc = NULL;
    int res;

    if(!z || !p || !b || !a)
        return ERROR_PTR;
    if(nz < 0 || np < 0)
        return ERROR_SIZE;
    if(nz > ord || np > ord)
        return ERROR_POLY_ORD;

    acc = (complex_t*) malloc((ord+1) * sizeof(complex_t));
    res = poly_z2a_cmplx(z, nz, ord, acc);
    if(res != RES_OK)
        goto exit_label;

    res = cmplx2re(acc, ord+1, b, NULL);
    if(res != RES_OK)
        goto exit_label;

    res = poly_z2a_cmplx(p, np, ord, acc);
    if(res != RES_OK)
        goto exit_label;

    res = cmplx2re(acc, ord+1, a, NULL);
    if(res != RES_OK)
        goto exit_label;

exit_label:
    if(acc)
        free(acc);
    return res;

}

