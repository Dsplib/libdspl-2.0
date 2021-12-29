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
\fn int freqz(double* b, double* a, int ord, double* w, int n, complex_t *h)

\brief Function calculates the digital filter frequency response 
\f$ H \left(e^{j \omega} \right)\f$ corresponds to transfer function \f$H(z)\f$.

Digital filter transfer function:
\f[
H(z) = \frac{\sum\limits_{k = 0}^{N} b_k  z^{-k}}
    {\sum\limits_{m = 0}^{N} a_m z^{-m}},
\f]
here \f$N\f$ --- filter order (parameter `ord`). \n

Frequency response \f$ H \left(e^{j \omega} \right)\f$ we can get
if substitute \f$z = e^{j \omega} \f$. \n

\param[in]  b
Pointer to the \f$ H(z) \f$ transfer function 
numerator coefficients vector. \n 
Vector size is `[ord+1 x 1]`. \n \n 

\param[in]  a
Pointer to the \f$H(z)\f$ transfer function 
denominator coefficients vector. \n 
Vector size is `[ord+1 x 1]`. \n \n 

\param[in]  ord
Filter order. \n
Transfer function \f$H(z)\f$ numerator 
and denominator coefficients number equals `ord+1`. \n \n 

\param[in]  w
Pointer to the normalized frequency of digital filter 
frequency response \f$ H \left(\mathrm{e}^{j\omega} \right) \f$. \n
Digital filter frequency response is \f$ 2\pi \f$-periodic function,
and vector `w` advisable to set from 0 to \f$ \pi \f$,
or from 0 to \f$ 2\pi \f$, or from \f$ -\pi \f$ to \f$ \pi \f$.
Vector size is `[n x 1]`. \n \n

\param[in]  n
Size of frequency vector `w`. \n  \n 

\param[out]  h
Pointer to the frequency response vector 
\f$ H \left(\mathrm{e}^{j\omega} \right) \f$, 
corresponds to normalized frequency `w`. \n 
Vector size is `[n x 1]`. \n 
Memory must be allocated. \n  \n

\return `RES_OK` if frequency response vector is calculated successfully. \n
Else \ref ERROR_CODE_GROUP "code error".

\author Sergey Bakhurin www.dsplib.org
***************************************************************************** */
#endif
#ifdef DOXYGEN_RUSSIAN
/*! ****************************************************************************
\ingroup FILTER_ANALYSIS_GROUP
\fn int freqz(double* b, double* a, int ord, double* w, int n, complex_t *h)

\brief  Расчет комплексного коэффициента передачи 
  \f$ H \left(e^{j \omega} \right)\f$ цифрового фильтра.

Функция рассчитывает значения комплексного коэффициента передачи 
\f$ H \left(e^{j \omega} \right)\f$  цифрового фильтра, заданного 
коэффициентами передаточной функции \f$H(z)\f$:

\f[
H(z) = \frac  {\sum_{k = 0}^{N} b_k  z^{-k}}
    {\sum_{m = 0}^{N} a_m z^{-m}},
\f]

где \f$N\f$ --- порядок фильтра (параметр `ord`). \n

Комплексный коэффициент передачи рассчитывается путем 
подстановки \f$z = e^{j \omega} \f$. \n

\param[in]  b
Указатель на вектор коэффициентов числителя 
передаточной функции \f$H(z)\f$. \n
Размер вектора `[ord+1 x 1]`. \n \n

\param[in]  a
Указатель на вектор коэффициентов знаменателя 
передаточной функции \f$H(z)\f$. \n
Размер вектора `[ord+1 x 1]`. \n \n

\param[in]  ord
Порядок фильтра. Количество коэффициентов числителя и знаменателя 
передаточной функции \f$H(z)\f$ равно `ord+1`. \n \n

\param[in]  w
Указатель на вектор значений  нормированной циклической частоты \f$\omega\f$, 
для которого будет рассчитан комплексный коэффициент передачи 
\f$ H \left(e^{j \omega} \right)\f$. \n
Размер вектора `[n x 1]`. \n \n

\param[in]  n
Размер вектора нормированной циклической частоты `w`. \n \n

\param[out]  h
Указатель на вектор комплексного коэффициента передачи 
\f$ H \left(e^{j \omega} \right)\f$, рассчитанного для 
циклической частоты `w`. \n
Размер вектора `[n x 1]`. \n
Память должна быть выделена. \n \n

\return
`RES_OK` Комплексный коэффициент передачи рассчитан успешно. \n
В противном случае \ref ERROR_CODE_GROUP "код ошибки". \n

\note
Комплексный коэффициент передачи \f$ H \left(e^{j \omega} \right)\f$ 
цифрового фильтра представляет собой \f$ 2 \pi-\f$периодическую функцию 
нормированной циклической частоты \f$\omega\f$.
Поэтому анализ цифровых фильтров целесообразно вести на одном периоде 
повторения \f$ H \left(e^{j \omega} \right)\f$, т.е. в интервале 
\f$\omega\f$ от 0 до \f$2 \pi\f$, или от \f$-\pi\f$ до \f$ \pi\f$.  \n
Кроме того известно, что для фильтра с вещественными коэффициентами 
\f$ H \left(e^{j \omega} \right) = H^* \left(e^{-j \omega} \right)\f$,
а значит, анализ цифрового фильтра с вещественными коэффициентами 
достаточно вести для нормированной частоты \f$\omega\f$ от 0 до \f$\pi\f$.

\author Бахурин Сергей www.dsplib.org
***************************************************************************** */
#endif
int DSPL_API freqz(double* b, double* a, int ord, double* w,
                   int n, complex_t *h)
{
    complex_t jw;
    complex_t *bc = NULL;
    complex_t *ac = NULL;
    complex_t num, den;
    double mag;
    int k;
    int res;

    if(!b || !w || !h)
        return ERROR_PTR;
    if(ord<0)
        return ERROR_FILTER_ORD;
    if(n<1)
        return ERROR_SIZE;


    bc = (complex_t*) malloc((ord+1) * sizeof(complex_t));
    res = re2cmplx(b, ord+1, bc);
    if( res!=RES_OK )
        goto exit_label;

    if(a)
    {
        /* IIR filter if a != NULL */
        ac = (complex_t*) malloc((ord+1) * sizeof(complex_t));
        res = re2cmplx(a, ord+1, ac);
        if( res!=RES_OK )
            goto exit_label;
        for(k = 0; k < n; k++)
        {
            RE(jw) =  cos(w[k]);
            IM(jw) = -sin(w[k]);
            res = polyval_cmplx(bc, ord, &jw, 1, &num);
            if(res != RES_OK)
                goto exit_label;
            res = polyval_cmplx(ac, ord, &jw, 1, &den);
            if(res != RES_OK)
                goto exit_label;
            mag = ABSSQR(den);
            if(mag == 0.0)
            {
                res = ERROR_DIV_ZERO;
                goto exit_label;
            }
            mag = 1.0 / mag;
            RE(h[k]) = CMCONJRE(num, den) * mag;
            IM(h[k]) = CMCONJIM(num, den) * mag;
        }
    }
    else
    {
        /* FIR filter if a == NULL */
        for(k = 0; k < n; k++)
        {
            RE(jw) =    cos(w[k]);
            IM(jw) = -sin(w[k]);
            res = polyval_cmplx(bc, ord, &jw, 1, h+k);
            if(res != RES_OK)
                goto exit_label;
        }
    }
    res = RES_OK;
exit_label:
    if(bc)
        free(bc);
    if(ac)
        free(ac);
    return res;
}

