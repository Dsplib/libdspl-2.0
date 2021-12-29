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
\fn int DSPL_API group_delay(double* b, double* a, int ord, int flag,
                         double* w, int n, double* tau)

\brief
Group delay calculation for digital or analog filter corresponds to 
\f$H(s)\f$, or \f$H(z)\f$ transfer function.

Group delay is describes as:
\f[
\tau_g(\omega) = - \frac{d\Phi(\omega)}{d\omega},
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
Pointer to the group delay vector. \n 
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
\fn int DSPL_API group_delay(double* b, double* a, int ord, int flag,
                         double* w, int n, double* tau)
\brief
Расчет группового времени запаздывания цифрового или аналогового фильтра.

Групповое время запаздывания определяется как:
\f[
\tau_g(\omega) = - \frac{d\Phi(\omega)}{d\omega},
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
Указатель на вектор групповой задержки. \n 
Размер вектора `[n x 1]`. \n
Память должна быть выделена. \n

\return
`RES_OK`  групповая задержка фильтра рассчитана успешно. \n
В противном случае \ref ERROR_CODE_GROUP "код ошибки". \n

\author Бахурин Сергей www.dsplib.org
***************************************************************************** */
#endif
int DSPL_API group_delay(double* pb, double* pa, int ord, int flag,
                         double* w, int n, double* tau)
{
    double a, b, c, d, da, db, dc, dd, f, e;
    int t, m;
    
    double *qa = NULL;
    
    if(!pb || !w || !tau || (!pa && (flag & DSPL_FLAG_ANALOG)))
        return ERROR_PTR;
    if(ord < 1)
        return ERROR_FILTER_ORD;
    if(n < 1)
        return ERROR_SIZE;


    if(pa)
      qa = pa;
    else
    {
        qa = (double*)malloc((ord+1) * sizeof(double));
        memset(qa, 0, (ord+1) * sizeof(double));
        qa[0] = 1.0;
    }
    
    for(t = 0; t < n; t++)
    {
        a = b = c = d = da = db = dc = dd = 0.0;
        if(flag & DSPL_FLAG_ANALOG)
        { 
            for(m = 0; m < ord+1; m+=4)
            {
                a  += pb[m] * pow(w[t], (double)m);
                c  += qa[m] * pow(w[t], (double)m);
                da += pb[m] * (double) m * pow(w[t], (double)(m-1));
                dc += qa[m] * (double) m * pow(w[t], (double)(m-1));
            }
            for(m = 2; m < ord+1; m+=4)
            {
                a  -= pb[m] * pow(w[t], (double)m);
                c  -= qa[m] * pow(w[t], (double)m);
                da -= pb[m] * (double) m * pow(w[t], (double)(m-1));
                dc -= qa[m] * (double) m * pow(w[t], (double)(m-1));
            }
            
            for(m = 1; m < ord+1; m+=4)
            {
                b  += pb[m] * pow(w[t], (double)m) ;
                d  += qa[m] * pow(w[t], (double)m) ;
                db += pb[m] * (double) m * pow(w[t], (double)(m-1)) ;
                dd += qa[m] * (double) m * pow(w[t], (double)(m-1)) ;
            }
            
            for(m = 3; m < ord+1; m+=4)
            {
                b  -= pb[m] * pow(w[t], (double)m) ;
                d  -= qa[m] * pow(w[t], (double)m) ;
                db -= pb[m] * (double) m * pow(w[t], (double)(m-1)) ;
                dd -= qa[m] * (double) m * pow(w[t], (double)(m-1)) ;
            }
            
        }
        else
        {
            for(m = 0; m < ord+1; m++)
            {
                a += pb[m] * cos(w[t]*(double)m);
                b -= pb[m] * sin(w[t]*(double)m);
                c += qa[m] * cos(w[t]*(double)m);
                d -= qa[m] * sin(w[t]*(double)m);
                  
                da -= pb[m] *(double)m * sin(w[t]*(double)m);
                db -= pb[m] *(double)m * cos(w[t]*(double)m);
                dc -= qa[m] *(double)m * sin(w[t]*(double)m);
                dd -= qa[m] *(double)m * cos(w[t]*(double)m);
            }
        }
        
        f = da * c + a * dc + db * d + b * dd;
        e = db * c + b * dc - da * d - a * dd;
        tau[t] = (f * (b * c - a * d) - e * (a * c + b * d)) / 
                 ((a * a + b * b) * (c * c + d * d));
    }
    
    if(qa != pa)
      free(qa);
    
    return RES_OK;
}
