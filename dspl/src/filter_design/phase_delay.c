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
#include <string.h>
#include <math.h>
#include "dspl.h"




#ifdef DOXYGEN_ENGLISH
/*! ****************************************************************************
\ingroup FILTER_ANALYSIS_GROUP
\fn int DSPL_API phase_delay(double* b, double* a, int ord, int flag,
                         double* w, int n, double* tau)

\brief
Phase delay calculation for digital or analog filter corresponds to 
\f$H(s)\f$, or \f$H(z)\f$ transfer function.

Group delay is describes as:
\f[
\tau_{\varphi}(\omega) = - \frac{\Phi(\omega)}{\omega},
\f]
here \f$\Phi(\omega)\f$ -- filter phase response, \f$\omega\f$ is angular
frequency for analog filter, or normalized frequency for digital filter.

\param[in]  b
Pointer to the \f$ H(s) \f$ or \f$H(z)\f$  transfer function 
numerator coefficients vector. \n 
Vector size is `[ord+1 x 1]`. \n \n 

\param[in]  a
Pointer to the \f$ H(s) \f$ or \f$H(z)\f$ transfer function 
denominator coefficients vector. \n 
Vector size is `[ord+1 x 1]`. \n \n 

\param[in]  ord
Filter order. \n
Transfer function \f$ H(s) \f$ or \f$H(z)\f$ numerator 
and denominator coefficients number equals `ord+1`. \n \n 

\param[in]  flag
Binary flags to set calculation rules: \n 
\verbatim
DSPL_FLAG_ANALOG  Coefficients corresponds to analog filter
DSPL_FLAG_DIGITAL Coefficients corresponds to digital filter
\endverbatim
\n \n

\param[in]  w
Pointer to the angular frequency  \f$ \omega \f$ (rad/s), 
which used for analog filter characteristics calculation
(flag sets as `DSPL_FLAG_ANALOG`). \n
For digital filter (flag sets as `DSPL_FLAG_DIGITAL`),
 parameter `w` describes normalized frequency of
frequency response \f$ H \left(\mathrm{e}^{j\omega} \right) \f$.
Digital filter frequency response is \f$ 2\pi \f$-periodic function,
and vector `w` advisable to set from 0 to \f$ \pi \f$,
or from 0 to \f$ 2\pi \f$, or from \f$ -\pi \f$ to \f$ \pi \f$.
Vector size is `[n x 1]`. \n \n

\param[in]  n
Size of frequency vector `w`. \n  \n 

\param[out]  tau
Pointer to the phase delay vector. \n 
Vector size is  `[n x 1]`. \n
Memory must be allocated. \n \n

\return
\return `RES_OK` if function is calculated successfully. \n
Else \ref ERROR_CODE_GROUP "code error".

\author Sergey Bakhurin www.dsplib.org
***************************************************************************** */
#endif
#ifdef DOXYGEN_RUSSIAN
/*! ****************************************************************************
\ingroup FILTER_ANALYSIS_GROUP
\fn int DSPL_API phase_delay(double* b, double* a, int ord, int flag,
                         double* w, int n, double* tau)
\brief
Расчет фазовой задержки цифрового или аналогового фильтра.

Фазовая задержка определяется как:
\f[
\tau_{\varphi}(\omega) = - \frac{\Phi(\omega)}{\omega},
\f]
где \f$\Phi(\omega)\f$ -- ФЧХ фильтра, \f$\omega\f$ циктическая частот в случае 
аналогового фильтра, или нормированная частота цифрового фильтра.

\param[in]  b
Указатель на вектор коэффициентов числителя передаточной функции 
аналогового фильтра \f$ H(s) \f$  или цифрового фильтра \f$ H(z) \f$. \n 
Размер вектора `[ord+1 x 1]`. \n \n 

\param[in]  a
Указатель на вектор коэффициентов числителя передаточной функции 
аналогового фильтра \f$ H(s) \f$  или цифрового фильтра \f$ H(z) \f$. \n 
Размер вектора `[ord+1 x 1]`. \n  
Параметр может быть `NULL`. В этом случае расчет производится для цифрового 
КИХ-фильтра с коэффициентами, заданными вектором `b`. \n\n

\param[in]  ord
Порядок фильтра. Количество коэффициентов 
числителя и знаменателя передаточной
функции \f$ H(s) \f$ или \f$ H(z) \f$ равно `ord+1`. \n  \n 

\param[in]  flag
Флаг который задает тип фильтра: \n 
\verbatim
DSPL_FLAG_ANALOG  Коэффициенты относятся к аналоговому фильтру
DSPL_FLAG_DIGITAL  Коэффициенты относятся к цифровому фильтру
\endverbatim

\param[in]  w
Указатель на вектор значений циклической частоты  \f$ \omega \f$ (рад/с), 
для которого будет рассчитаны  АЧХ, ФЧХ и ГВЗ аналогового фильтра, 
если установлен  флаг `DSPL_FLAG_ANALOG`. \n
В случае если флаг `DSPL_FLAG_ANALOG` не установлен, то  вектор частоты `w` 
используется как нормированная частота комплексного коэффициента передачи 
\f$ H \left(\mathrm{e}^{j\omega} \right) \f$  цифрового фильтра. \n
В этом случае характеристика цифрового фильтра является 
\f$ 2\pi \f$-периодической, и вектор частоты может содержать 
произвольные значения, однако целесообразно задавать 
его от 0 до \f$ \pi \f$, а такжет от 0 до \f$ 2\pi \f$, или 
от \f$ -\pi \f$ до \f$ \pi \f$. \n
Размер вектора `[n x 1]`. \n \n 

\param[in]  n
Размер вектора циклической частоты `w`. \n  \n 

\param[out]  tau
Указатель на вектор фазовой задержки. \n 
Размер вектора `[n x 1]`. \n
Память должна быть выделена. \n

\return
`RES_OK` фазовая задержка фильтра рассчитана успешно. \n
В противном случае \ref ERROR_CODE_GROUP "код ошибки". \n

\author Бахурин Сергей www.dsplib.org
***************************************************************************** */
#endif
int DSPL_API phase_delay(double* b, double* a, int ord, int flag,
                         double* w, int n, double* tau)
{
    int err, i;
    double *phi =  NULL;

    if(n > 0)
        phi = (double*)malloc(n*sizeof(double));
    else
        return ERROR_SIZE;

    err = filter_freq_resp(b, a, ord, w, n, flag | DSPL_FLAG_UNWRAP, 
	                                                           NULL, phi, NULL);
    if(err!=RES_OK)
        goto exit_label;
    for(i = 0; i < n; i++)
    {
        tau[i] = w[i] ? ( - phi[i] / w[i]) : ( - phi[i] / (w[i] + 1E-9) );
    }
exit_label:
    if(phi)
        free(phi);
    return err;
}

