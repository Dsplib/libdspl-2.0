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

#endif
#ifdef DOXYGEN_RUSSIAN

#endif
double DSPL_API filter_ws1(int ord, double rp, double rs, int type)
{
    double es2, ep2, gs2, x, ws;

    if(ord<1 || rp < 0.0 || rs < 0.0)
        return -1.0;

    es2 = pow(10.0, rs*0.1) - 1.0;
    ep2 = pow(10.0, rp*0.1) - 1.0;
    gs2 = 1.0 / (1.0 + es2);

    x = (1.0 - gs2) / (gs2 * ep2);

    switch( type & DSPL_FILTER_APPROX_MASK)
    {
        case DSPL_FILTER_BUTTER:
            ws = pow(x, 0.5 / (double)ord);
            break;
        case DSPL_FILTER_CHEBY1:
        case DSPL_FILTER_CHEBY2:
            x = sqrt(x) + sqrt(x - 1.0);
            x = log(x) / (double)ord;
            ws    = 0.5 * (exp(-x) + exp(x));
            break;
        case DSPL_FILTER_ELLIP:
        {
            double k, k1;
            complex_t y, z;
            int res;
            k = sqrt(ep2 / es2);
            res = ellip_modulareq(rp, rs, ord, &k1);
            if(res != RES_OK)
            {
                ws = -1.0;
                break;
            }
            RE(z) = sqrt(x);
            IM(z) = 0.0;

            res = ellip_acd_cmplx(&z, 1, k, &y);
            if(res != RES_OK)
            {
                ws = -1.0;
                break;
            }
            RE(y) /= (double)ord;
            IM(y) /= (double)ord;
            res = ellip_cd_cmplx(&y, 1, k1, &z);
            if(res != RES_OK)
            {
                ws = -1.0;
                break;
            }
            ws = RE(z);
            break;
        }
        default:
            ws    = -1.0;
            break;
    }
    return ws;
}


#ifdef DOXYGEN_ENGLISH

#endif
#ifdef DOXYGEN_RUSSIAN

#endif
int DSPL_API low2bp(double* b, double* a, int ord,
                    double w0, double wpl, double wph,
                    double* beta, double* alpha)
{

    double num[3] = {0.0, 0.0, 1.0};
    double den[3] = {0.0, 0.0, 0.0};

    if(!b || !a || !beta || !alpha)
        return ERROR_PTR;
    if(ord < 1)
        return ERROR_FILTER_ORD;
    if(w0 <= 0.0 || wpl <= 0.0 || wph <= 0.0 || wph <= wpl)
        return ERROR_FILTER_FT;

    num[0] = (wph * wpl) / (w0 * w0);
    den[1] = (wph - wpl) / w0;

    return ratcompos(b, a, ord, num, den, 2, beta, alpha);
}






#ifdef DOXYGEN_ENGLISH

#endif
#ifdef DOXYGEN_RUSSIAN

#endif
int DSPL_API low2bs(double* b, double* a, int ord,
                    double w0, double wsl, double wsh,
                    double* beta, double* alpha)
{

    double den[3] = {0.0, 0.0, 1.0};
    double num[3] = {0.0, 0.0, 0.0};

    if(!b || !a || !beta || !alpha)
        return ERROR_PTR;
    if(ord < 1)
        return ERROR_FILTER_ORD;
    if(w0 <= 0.0 || wsl <= 0.0 || wsh <= 0.0 || wsh <= wsl)
        return ERROR_FILTER_FT;

    den[0] = (wsh * wsl) / (w0 * w0);
    num[1] = (wsh - wsl) / w0;

    return ratcompos(b, a, ord, num, den, 2, beta, alpha);
}





#ifdef DOXYGEN_ENGLISH
/*! ****************************************************************************
\ingroup IIR_FILTER_DESIGN_GROUP
\fn int  low2high (double* b, double* a, int ord, double w0, double w1,
                   double* beta, double* alpha)
\brief Lowpass to highpass filter frequency transform

Function transforms lowpass filter transfer function \f$ H(s) \f$ 
to the highpass filter transfer function \f$ F(s) \f$.

Filter order, magnitude ripple in passband and stopband 
supression still the same.

\param[in]  b
Pointer to the lowpass filter transfer function \f$H(s)\f$ numerator 
coefficients vector. \n
Vector size is `[ord+1 x 1]`. \n 
\n

\param[in]  a
Pointer to the lowpass filter transfer function \f$H(s)\f$ denominator 
coefficients vector. \n
Vector size is `[ord+1 x 1]`. \n 
\n

\param[in]  ord
Filter order. \n 
\n

\param[in]  w0
Lowpass filter cutoff frequency. \n 
\n

\param[in]  w1
Highpass filter cutoff frequency after transformation. \n 
\n

\param[in,out]  beta
Pointer to the highwpass filter transfer function \f$F(s)\f$ numerator 
coefficients vector after transformation. \n
Vector size is `[ord+1 x 1]`. \n 
Memory must be allocated. \n
\n

\param[in,out]  alpha
Pointer to the highwpass filter transfer function \f$F(s)\f$ denominator 
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
\fn int  low2high (double* b, double* a, int ord, double w0, double w1,
                   double* beta, double* alpha)
\brief  Частотное преобразование ФНЧ-ФВЧ

Функция производит перобразование передаточной функции \f$ H(s) \f$ 
аналогового ФНЧ с частотой среза `w0` рад/c 
в передаточную функцию \f$ F(s) \f$ аналоговго ФВЧ с частотой среза `w1` рад/c.

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
Порядок исходного фильтра и фильтра после переобразования. \n 
\n

\param[in]  w0
Частота среза исходного ФНЧ. \n 
\n

\param[in]  w1
Требуемая частота среза  ФВЧ после преобразования. \n 
\n

\param[in,out]  beta
Указатель на вектор коэффициентов числителя передаточной функции \f$F(s)\f$
ФВЧ после преобразования. \n
Размер вектора `[ord+1 x 1]`. \n
Память должна быть выделена. \n 
\n

\param[in,out]  alpha
Указатель на вектор коэффициентов знаменателя передаточной функции \f$F(s)\f$
аналогового ФВЧ после преобразования. \n
Размер вектора `[ord+1 x 1]`. \n
Память должна быть выделена. \n 
\n

\return
`RES_OK` --- преобразование рассчитано успешно. \n
В противном случае \ref ERROR_CODE_GROUP "код ошибки". \n

\author Бахурин Сергей www.dsplib.org
***************************************************************************** */
#endif
int DSPL_API low2high(double* b, double* a, int ord, double w0, double w1,
                      double* beta, double* alpha)
{

    double num[2] = {0.0, 0.0};
    double den[2] = {0.0, 1.0};

    if(!b || !a || !beta || !alpha)
        return ERROR_PTR;
    if(ord < 1)
        return ERROR_FILTER_ORD;
    if(w0 <= 0.0 || w1 <= 0.0)
        return ERROR_FILTER_FT;

    num[0] = w1 / w0;

    return ratcompos(b, a, ord, num, den, 1, beta, alpha);
}








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




#ifdef DOXYGEN_ENGLISH
/*! ****************************************************************************
\ingroup IIR_FILTER_DESIGN_GROUP
\fn int ratcompos(  double* b, double* a, int n, 
                    double* c, double* d, int p, 
                    double* beta, double* alpha)  
\brief  Rational composition

Function calcultes composition \f$Y(s) = (H \circ F)(s) = H(F(s))\f$, here

\f[
H(s) = \frac{\sum\limits_{m = 0}^{n} b_m s^m}
{\sum\limits_{k = 0}^{n} a_k s^k}, \quad 
F(s) = \frac{\sum\limits_{m = 0}^{p} d_m s^m}
{\sum\limits_{k = 0}^{p} c_k s^k}, \quad 
Y(s) = \frac{\sum\limits_{m = 0}^{n p} \beta_m s^m}
{\sum\limits_{k = 0}^{n  p} \alpha_k s^k}
\f]

This function is using for filter frequency transform.

\param[in]  b
Pointer to the \f$H(s)\f$ polynomial function 
numerator coefficients vector. \n
Vector size is `[n+1 x 1]`. \n
\n

\param[in]  a
Pointer to the \f$H(s)\f$ polynomial function 
denominator coefficients vector. \n
Vector size is `[n+1 x 1]`. \n
\n

\param[in]  n
Order of \f$H(s)\f$ numerator and denominator polynomials. \n 
\n

\param[in]  c
Pointer to the \f$F(s)\f$ polynomial function 
numerator coefficients vector. \n
Vector size is `[p+1 x 1]`. \n
\n

\param[in]  d
Pointer to the \f$F(s)\f$ polynomial function 
denominator coefficients vector. \n
Vector size is `[p+1 x 1]`. \n
\n

\param[in]  p
Order of \f$F(s)\f$ numerator and denominator polynomials. \n 
\n

\param[in,out]  beta
Pointer to the numerator coefficients vector of 
\f$Y(s) = (H \circ F)(s)\f$. \n
Vector size is `[n*p+1 x 1]`. \n
Memory must be allocated. \n 
\n

\param[in,out]  alpha
Pointer to the denominator coefficients vector of 
\f$Y(s) = (H \circ F)(s)\f$. \n
Vector size is `[n*p+1 x 1]`. \n
Memory must be allocated. \n 
\n


\return
`RES_OK` if rational composition is calculated successfully. \n 
Else \ref ERROR_CODE_GROUP "code error".

\author Sergey Bakhurin www.dsplib.org
***************************************************************************** */
#endif
#ifdef DOXYGEN_RUSSIAN
/*! ****************************************************************************
\ingroup IIR_FILTER_DESIGN_GROUP
\fn int ratcompos(  double* b, double* a, int n, 
                    double* c, double* d, int p, 
                    double* beta, double* alpha)  
\brief  Рациональная композиця 

Функция рассчитывает композицию вида \f$Y(s) = (H \circ F)(s) = H(F(s))\f$, где

\f[
H(s) = \frac{\sum\limits_{m = 0}^{n} b_m s^m}
{\sum\limits_{k = 0}^{n} a_k s^k}, \quad 
F(s) = \frac{\sum\limits_{m = 0}^{p} d_m s^m}
{\sum\limits_{k = 0}^{p} c_k s^k}, \quad 
Y(s) = \frac{\sum\limits_{m = 0}^{n p} \beta_m s^m}
{\sum\limits_{k = 0}^{n  p} \alpha_k s^k}
\f]

Функция рациональной композиции необходима для произведения частотных
преобразований передаточных характеристик аналоговых и цифровых фильтров,
а также для билинейного преобразования передаточных характеристик аналоговых 
фильтров в соответствующие передаточные характеристики цифровых фильтров.

\param[in]  b
Указатель на вектор коэффициентов числителя функции \f$H(s)\f$. \n
Размер вектора `[n+1 x 1]`. \n
Память должна быть выделена. \n 
\n

\param[in]  a
Указатель на вектор коэффициентов знаменателя функции \f$H(s)\f$. \n
Размер вектора `[n+1 x 1]`. \n
Память должна быть выделена. \n 
\n

\param[in]  n
Порядок полиномов рациональной функции \f$H(s)\f$. \n 
\n

\param[in]  c
Указатель на вектор коэффициентов числителя функции \f$F(s)\f$. \n
Размер вектора `[p+1 x 1]`. \n
Память должна быть выделена. \n 
\n

\param[in]  d
Указатель на вектор коэффициентов знаменателя функции \f$F(s)\f$. \n
Размер вектора `[p+1 x 1]`. \n
Память должна быть выделена. \n 
\n

\param[in]  p
Порядок полиномов рациональной 
функции \f$F(s)\f$. \n 
\n

\param[in,out]  beta
Указатель на вектор коэффициентов 
числителя функции \f$Y(s) = (H \circ F)(s)\f$. \n
Размер вектора `[n*p+1 x 1]`. \n
Память должна быть выделена. \n 
\n

\param[in,out]  alpha
Указатель на вектор коэффициентов знаменателя 
функции \f$Y(s) = (H \circ F)(s)\f$. \n
Размер вектора `[n*p+1 x 1]`. \n
Память должна быть выделена. \n 
\n
      

\return
`RES_OK` --- Рациональная композиция рассчитана успешно. \n
В противном случае \ref ERROR_CODE_GROUP "код ошибки". \n

\author Бахурин Сергей www.dsplib.org
***************************************************************************** */
#endif
int DSPL_API ratcompos(double* b, double* a, int n,
                       double* c, double* d, int p,
                       double* beta, double* alpha)
{

    int k2, i, k,    pn, pd, ln, ld, k2s, nk2s;
    double *num = NULL, *den = NULL, *ndn = NULL, *ndd = NULL;
    int res;

    if (!a || !b || !c || !d || !beta || !alpha)
    {
        res = ERROR_PTR;
        goto exit_label;
    }
    if(n < 1 || p < 1)
    {
        res =    ERROR_SIZE;
        goto exit_label;
    }

    k2   = (n*p)+1;
    k2s  = k2*sizeof(double);     /* alpha and beta size */
    nk2s = (n+1)*k2*sizeof(double); /* num, den, ndn and ndd size */

    num = (double*)malloc(nk2s);
    den = (double*)malloc(nk2s);
    ndn = (double*)malloc(nk2s);
    ndd = (double*)malloc(nk2s);

    memset(num, 0, nk2s);
    memset(den, 0, nk2s);
    memset(ndn, 0, nk2s);
    memset(ndd, 0, nk2s);


    num[0] = den[0] = 1.0;
    pn = 0;
    ln = 1;
    for(i = 1; i < n+1; i++)
    {
        res = conv(num+pn, ln, c, p+1, num+pn+k2);
        if(res!=RES_OK)
            goto exit_label;
        res = conv(den+pn, ln, d, p+1, den+pn+k2);
        if(res!=RES_OK)
            goto exit_label;
        pn += k2;
        ln += p;
    }

    pn = 0;
    pd = n*k2;
    ln = 1;
    ld = k2;

    for (i = 0; i < n+1; i++)
    {
        res = conv(num + pn, ln, den + pd, ld, ndn + i*k2);
        if(res!=RES_OK)
            goto exit_label;
        ln += p;
        ld -= p;
        pn += k2;
        pd -= k2;
    }

    for (i = 0; i < n+1; i++)
    {
        for (k = 0; k < k2; k++)
        {
            ndd[i*k2 + k] = ndn[i*k2 + k] * a[i];
            ndn[i*k2 + k] *= b[i];
        }
    }


    memset(alpha, 0, k2s);
    memset(beta,  0, k2s);

    for (k = 0; k < k2; k++)
    {
        for (i = 0; i < n+1; i++)
        {
            beta[k]  += ndn[i*k2 + k];
            alpha[k] += ndd[i*k2 + k];
        }
    }

    res = RES_OK;

exit_label:
    if(num)
        free(num);
    if(den)
        free(den);
    if(ndn)
        free(ndn);
    if(ndd)
        free(ndd);

    return res;
}

