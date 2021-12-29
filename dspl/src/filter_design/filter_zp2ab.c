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

