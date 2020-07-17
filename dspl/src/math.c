/*
* Copyright (c) 2015-2020 Sergey Bakhurin
* Digital Signal Processing Library [http://dsplib.org]
*
* This file is part of DSPL.
*
* is free software: you can redistribute it and/or modify
* it under the terms of the GNU General Public License as published by
* the Free Software Foundation, either version 3 of the License, or
* (at your option) any later version.
*
* DSPL is distributed in the hope that it will be useful,
* but WITHOUT ANY WARRANTY; without even the implied warranty of
* MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.    See the
* GNU General Public License for more details.
*
* You should have received a copy of the GNU General Public License
* along with Foobar.    If not, see <http://www.gnu.org/licenses/>.
*/


#include <stdio.h>
#include <stdlib.h>
#include "dspl.h"




#if DOXYGEN_ENGLISH
/*! ****************************************************************************
\ingroup SPEC_MATH_TRIG_GROUP
\fn int acos_cmplx(complex_t* x, int n, complex_t *y)
\brief    The inverse of the cosine function the complex vector argument `x`

Function calculates the inverse of the cosine function as: \n

\f[
\textrm{Arccos}(x) = \frac{\pi}{2} - \textrm{Arcsin}(x) = 
\frac{\pi}{2} -j \textrm{Ln}\left( j x + \sqrt{1 - x^2} \right)
\f]


\param[in]    x
Pointer to the argument vector `x`. \n
Vector size is `[n x 1]`. \n
\n

\param[in]    n
Input vector `x` and the inverse cosine vector `y` size. \n
\n


\param[out] y
Pointer to the output complex vector `y`, 
corresponds to the input vector `x`. \n 
Vector size is `[n x 1]`. \n
Memory must be allocated. \n
\n

\return
`RES_OK` if function calculated successfully. \n
Else \ref ERROR_CODE_GROUP "code error". \n

Example: \n
\code{.cpp}
    complex_t x[3] = {{1.0, 2.0}, {3.0, 4.0}, {5.0, 6.0}};
    complex_t y[3];
    int k;
    
    acos_cmplx(x, 3, y);
    
    for(k = 0; k < 3; k++)
        printf("acos_cmplx(%.1f%+.1fj) = %.3f%+.3fj\n", 
                RE(x[k]), IM(x[k]), RE(y[k]), IM(y[k]));
\endcode 
\n

Output is: \n
\verbatim
acos_cmplx(1.0+2.0j) = 1.144-1.529j
acos_cmplx(3.0+4.0j) = 0.937-2.306j
acos_cmplx(5.0+6.0j) = 0.880-2.749j
\endverbatim

\author
Sergey Bakhurin www.dsplib.org
***************************************************************************** */
#endif

#if DOXYGEN_RUSSIAN
/*! ****************************************************************************
\ingroup SPEC_MATH_TRIG_GROUP
\fn int acos_cmplx(complex_t* x, int n, complex_t *y)
\brief  Арккосинус комплексного аргумента `x`

Функция рассчитывает значения арккосинуса комплексного аргумента, 
заданного вектором `x` длины `n`:  \n
\f[
\textrm{Arccos}(x) = \frac{\pi}{2} - \textrm{Arcsin}(x) = 
\frac{\pi}{2} -j \textrm{Ln}\left( j x + \sqrt{1 - x^2} \right)
\f]  


\param[in]  x
Указатель на вектор аргумента комплексного арккосинуса. \n
Размер вектора `[n x 1]`.  \n \n

\param[in]  n
Размер входного и выходного векторов `x` и `y`. \n \n


\param[out] y
Указатель на вектор значений комплексного арккосинуса,
соответствующего входному вектору `x`. \n
Размер массива `[n x 1]`.  \n
Память должна быть выделена.  \n \n

\return
`RES_OK` если значение функции рассчитано успешно   .  \n
В противном случае \ref ERROR_CODE_GROUP "код ошибки": \n

\note
Функция может использоваться для расчета арккосинуса аргумента 
большего единицы, когда вещественная функция `acos` не определена.

Например при выполнении следующего кода 
\code{.cpp}
    complex_t x[3] = {{1.0, 2.0}, {3.0, 4.0}, {5.0, 6.0}};
    complex_t y[3];
    int k;
  
    acos_cmplx(x, 3, y);
  
    for(k = 0; k < 3; k++)
        printf("acos_cmplx(%.1f%+.1fj) = %.3f%+.3fj\n", 
               RE(x[k]), IM(x[k]), RE(y[k]), IM(y[k]));
\endcode 
\n

Результатом работы будет

\verbatim
acos_cmplx(1.0+2.0j) = 1.144-1.529j
acos_cmplx(3.0+4.0j) = 0.937-2.306j
acos_cmplx(5.0+6.0j) = 0.880-2.749j
\endverbatim

\author
Бахурин Сергей www.dsplib.org 
***************************************************************************** */
#endif
int DSPL_API acos_cmplx(complex_t* x, int n, complex_t *y)
{
    int k, res;
    double pi2 = 0.5 * M_PI;

    res = asin_cmplx(x, n, y);
    if(res != RES_OK)
        return res;

    for(k = 0; k < n; k++)
    {
        RE(y[k]) = pi2 - RE(y[k]);
        IM(y[k]) =     - IM(y[k]);
    }
    return RES_OK;
}




/******************************************************************************
The inverse of the sine function the complex vector argument `x`
*******************************************************************************/
int DSPL_API asin_cmplx(complex_t* x, int n, complex_t *y)
{
    int k;
    complex_t tmp;
    if(!x || !y)
        return ERROR_PTR;
    if(n < 1)
        return ERROR_SIZE;

    for(k = 0; k < n; k++)
    {
        RE(tmp) = 1.0 - CMRE(x[k], x[k]);     /* 1-x[k]^2    */
        IM(tmp) =         - CMIM(x[k], x[k]); /* 1-x[k]^2    */
        sqrt_cmplx(&tmp, 1, y+k);             /* sqrt(1 - x[k]^2) */
        RE(y[k]) -= IM(x[k]);      /* j * x[k] + sqrt(1 - x[k]^2) */
        IM(y[k]) += RE(x[k]);      /* j * x[k] + sqrt(1 - x[k]^2) */
        log_cmplx(y+k, 1, &tmp);   /* log( j * x[k] + sqrt(1 - x[k]^2) ) */
        RE(y[k]) =    IM(tmp);     /* -j * log( j * x[k] + sqrt(1 - x[k]^2) ) */
        IM(y[k]) = -RE(tmp);       /* -j * log( j * x[k] + sqrt(1 - x[k]^2) ) */
    }
    return RES_OK;
}




/*******************************************************************************
Modified Bessel Function of the First Kind – I0(x) [1]

[1] Rational Approximations for the Modified Bessel Function
        of the First Kind – I0(x) for Computations with Double Precision
        by PAVEL HOLOBORODKO on NOVEMBER 11, 2015

        https://www.advanpix.com/2015/11/11/
*******************************************************************************/
int DSPL_API bessel_i0(double* x, int n, double* y)
{
    double P16[17] = { 1.0000000000000000000000801e+00,
                       2.4999999999999999999629693e-01,
                       2.7777777777777777805664954e-02,
                       1.7361111111111110294015271e-03,
                       6.9444444444444568581891535e-05,
                       1.9290123456788994104574754e-06,
                       3.9367598891475388547279760e-08,
                       6.1511873265092916275099070e-10,
                       7.5940584360755226536109511e-12,
                       7.5940582595094190098755663e-14,
                       6.2760839879536225394314453e-16,
                       4.3583591008893599099577755e-18,
                       2.5791926805873898803749321e-20,
                       1.3141332422663039834197910e-22,
                       5.9203280572170548134753422e-25,
                       2.0732014503197852176921968e-27,
                       1.1497640034400735733456400e-29};

    double P22[23] = { 3.9894228040143265335649948e-01,
                       4.9867785050353992900698488e-02,
                       2.8050628884163787533196746e-02,
                       2.9219501690198775910219311e-02,
                       4.4718622769244715693031735e-02,
                       9.4085204199017869159183831e-02,
                      -1.0699095472110916094973951e-01,
                       2.2725199603010833194037016e+01,
                      -1.0026890180180668595066918e+03,
                       3.1275740782277570164423916e+04,
                      -5.9355022509673600842060002e+05,
                       2.6092888649549172879282592e+06,
                       2.3518420447411254516178388e+08,
                      -8.9270060370015930749184222e+09,
                       1.8592340458074104721496236e+11,
                      -2.6632742974569782078420204e+12,
                       2.7752144774934763122129261e+13,
                      -2.1323049786724612220362154e+14,
                       1.1989242681178569338129044e+15,
                      -4.8049082153027457378879746e+15,
                       1.3012646806421079076251950e+16,
                      -2.1363029690365351606041265e+16,
                       1.6069467093441596329340754e+16};

    double x2;
    int k;

    if(!x || !y)
        return ERROR_PTR;
    if(n < 1)
        return ERROR_SIZE;

    for(k =0; k < n; k++)
    {
        if(x[k] < 0.0)
            return ERROR_NEGATIVE;

        if(x[k] < 7.75)
        {
            x2 = x[k] * x[k] * 0.25;
            polyval(P16, 16, &x2, 1, y+k);
            y[k] = x2 * y[k] + 1.0;
        }
        else
        {
            x2 = 1.0 / x[k];
            polyval(P22, 22, &x2, 1, y+k);
            y[k] *= exp(x[k]) / sqrt(x[k]);
        }
    }
    return RES_OK;
}


/*******************************************************************************
module operator for double
*******************************************************************************/
double DSPL_API dmod (double x, double y)
{
    if(y == 0.0)
        return x;
    return x - floor(x/y) * y;
}




/******************************************************************************
The cosine function the complex vector argument `x`
*******************************************************************************/
int DSPL_API cos_cmplx(complex_t* x, int n, complex_t *y)
{
    int k;
    double ep, em, sx, cx;
    if(!x || !y)
        return ERROR_PTR;
    if(n < 1)
        return ERROR_SIZE;

    for(k = 0; k < n; k++)
    {
        ep = exp( IM(x[k]));
        em = exp(-IM(x[k]));
        sx = 0.5 * sin(RE(x[k]));
        cx = 0.5 * cos(RE(x[k]));
        RE(y[k]) = cx * (em + ep);
        IM(y[k]) = sx * (em - ep);
    }
    return RES_OK;
}




/******************************************************************************
The logarithm function the complex vector argument `x`
*******************************************************************************/
int DSPL_API log_cmplx(complex_t* x, int n, complex_t *y)
{
    int k;
    if(!x || !y)
        return ERROR_PTR;
    if(n < 1)
        return ERROR_SIZE;

    for(k = 0; k < n; k++)
    {
        RE(y[k]) = 0.5 * log(ABSSQR(x[k]));
        IM(y[k]) = atan2(IM(x[k]), RE(x[k]));
    }
    return RES_OK;
}







/******************************************************************************
\brief    The sine function the complex vector argument `x`
*******************************************************************************/
int DSPL_API sin_cmplx(complex_t* x, int n, complex_t *y)
{
    int k;
    double ep, em, sx, cx;
    if(!x || !y)
        return ERROR_PTR;
    if(n < 1)
        return ERROR_SIZE;

    for(k = 0; k < n; k++)
    {
        ep = exp( IM(x[k]));
        em = exp(-IM(x[k]));
        sx = 0.5 * sin(RE(x[k]));
        cx = 0.5 * cos(RE(x[k]));
        RE(y[k]) = sx * (em + ep);
        IM(y[k]) = cx * (ep - em);
    }
    return RES_OK;
}






/*******************************************************************************
sinc(x) = sin(pi*x)/(pi*x)
*******************************************************************************/
int DSPL_API sinc(double* x, int n, double a, double* y)
{
    int k;

    if(!x || !y)
        return ERROR_PTR;
    if(n<1)
        return ERROR_SIZE;

    for(k = 0; k < n; k++)
        y[k] = (x[k]==0.0) ? 1.0 : sin(a*x[k])/(a*x[k]);

    return RES_OK;


}



/*******************************************************************************
Sine integral
--------------------------------------------------------------------------------
This function uses Padé approximants of the convergent Taylor series [1]


[1]
https://www.sciencedirect.com/science/article/pii/S221313371500013X?via%3Dihub

*******************************************************************************/
int DSPL_API sine_int(double* x, int n, double* si)
{
    int k, sgn, p;
    double num, den, y, x2, x22, z, f, g;

    double A[8] =     {+1.00000000000000000E0,
                       -4.54393409816329991E-2,
                       +1.15457225751016682E-3,
                       -1.41018536821330254E-5,
                       +9.43280809438713025E-8,
                       -3.53201978997168357E-10,
                       +7.08240282274875911E-13,
                       -6.05338212010422477E-16};



    double B[7]    =    {+1.0,
                         +1.01162145739225565E-2,
                         +4.99175116169755106E-5,
                         +1.55654986308745614E-7,
                         +3.28067571055789734E-10,
                         +4.50490975753865810E-13,
                         +3.21107051193712168E-16};



    double FA[11] = {+1.000000000000000000000E0,
                     +7.444370681619367006180E2,
                     +1.963963728951468698010E5,
                     +2.377503101254318340340E7,
                     +1.430734038212746368880E9,
                     +4.33736238870432522765E10,
                     +6.40533830574022022911E11,
                     +4.20968180571076940208E12,
                     +1.00795182980368574617E13,
                     +4.94816688199951963482E12,
                     -4.94701168645415959931E11};

    double FB[10] = {+1.000000000000000000000E0,
                     +7.464370681619276780310E2,
                     +1.978652470315839514500E5,
                     +2.415356701651268451440E7,
                     +1.474789521929854649580E9,
                     +4.58595115847765779830E10,
                     +7.08501308149515401563E11,
                     +5.06084464593475076774E12,
                     +1.43468549171581016479E13,
                     +1.11535493509914254097E13};



    double GA[11] = {+1.000000000000000000E0,
                     +8.135952011516861500E2,
                     +2.352391816264782000E5,
                     +3.125575707957787310E7,
                     +2.062975951467633540E9,
                     +6.83052205423625007E10,
                     +1.09049528450362786E12,
                     +7.57664583257834349E12,
                     +1.81004487464664575E13,
                     +6.43291613143049485E12,
                     -1.36517137670871689E12};


    double GB[10] = {+1.000000000000000000E0,
                     +8.195952011514515640E2,
                     +2.400367528355787770E5,
                     +3.260266616470908220E7,
                     +2.233555432780993600E9,
                     +7.87465017341829930E10,
                     +1.39866710696414565E12,
                     +1.17164723371736605E13,
                     +4.01839087307656620E13,
                     +3.99653257887490811E13};

    if(!x || !si)
        return ERROR_PTR;
    if(n<1)
        return ERROR_SIZE;


    for(p = 0; p < n; p++)
    {
        sgn = x[p] > 0.0 ?    0 : 1;
        y     = x[p] < 0.0 ? -x[p] : x[p];

        if(y < 4)
        {
            x2 = y * y;
            z = 1.0;
            num = 0.0;
            for(k = 0; k < 8; k++)
            {
                num += A[k] * z;
                z*=x2;
            }
            z = 1.0;
            den = 0.0;
            for(k = 0; k < 7; k++)
            {
                den += B[k]*z;
                z*=x2;
            }
            si[p] = x[p] * num/den;
        }
        else
        {

            x2 = 1.0/y;
            x22 = x2*x2;
            z = 1.0;
            num = 0.0;
            for(k = 0; k < 11; k++)
            {
                num += FA[k] * z;
                z*=x22;
            }
            z = 1.0;
            den = 0.0;
            for(k = 0; k < 10; k++)
            {
                den += FB[k]*z;
                z*=x22;
            }

            f = x2 * num / den;

            z = 1.0;
            num = 0.0;
            for(k = 0; k < 11; k++)
            {
                num += GA[k] * z;
                z*=x22;
            }
            z = 1.0;
            den = 0.0;
            for(k = 0; k < 10; k++)
            {
                den += GB[k]*z;
                z*=x22;
            }

            g = x22 * num / den;

            si[p] = sgn ? f * cos(y) + g * sin(y) - M_PI * 0.5 :
                          M_PI * 0.5 - f * cos(y) - g * sin(y);
        }
    }
    return RES_OK;
}








/******************************************************************************
\ingroup SPEC_MATH_COMMON_GROUP
\fn int sqrt_cmplx(complex_t* x, int n, complex_t *y)
\brief Square root of the complex vector argguument `x`.

Function calculates square root value of vector `x` length `n`:    \n
\f[
y(k) = \sqrt{x(k)}, \qquad k = 0 \ldots n-1.
\f]


\param[in]    x     Pointer to the input complex vector `x`. \n
                                Vector size is `[n x 1]`.    \n \n

\param[in]    n     Size of input and output vectors `x` and `y`. \n \n


\param[out] y     Pointer to the square root vector `y`. \n
                                Vector size is `[n x 1]`.    \n
                                Memory must be allocated.    \n \n

\return `RES_OK` if function is calculated successfully.    \n
Else \ref ERROR_CODE_GROUP "code error". \n

Example
\code{.cpp}
    complex_t x[3] = {{1.0, 2.0}, {3.0, 4.0}, {5.0, 6.0}};
    complex_t y[3]
    int k;

    sqrt_cmplx(x, 3, y);

    for(k = 0; k < 3; k++)
        printf("sqrt_cmplx(%.1f%+.1fj) = %.3f%+.3fj\n",
                        RE(x[k]), IM(x[k]), RE(y[k]), IM(y[k]));

 \endcode
 \n

Результатом работы будет

\verbatim
sqrt_cmplx(1.0+2.0j) = 1.272+0.786j
sqrt_cmplx(3.0+4.0j) = 2.000+1.000j
sqrt_cmplx(5.0+6.0j) = 2.531+1.185j
\endverbatim

\author Sergey Bakhurin www.dsplib.org
*******************************************************************************/
int DSPL_API sqrt_cmplx(complex_t* x, int n, complex_t *y)
{
    int k;
    double r, zr, at;
    complex_t t;
    if(!x || !y)
        return ERROR_PTR;
    if(n < 1)
        return ERROR_SIZE;

    for(k = 0; k < n; k++)
    {
        r = ABS(x[k]);
        if(r == 0.0)
        {
            RE(y[k]) = 0.0;
            IM(y[k]) = 0.0;
        }
        else
        {
            RE(t) = RE(x[k]) + r;
            IM(t) = IM(x[k]);
            at = ABS(t);
            if(at == 0.0)
            {
                RE(y[k]) = 0.0;
                IM(y[k]) = sqrt(r);
            }
            else
            {
                zr = 1.0 / ABS(t);
                r = sqrt(r);
                RE(y[k]) = RE(t) * zr * r;
                IM(y[k]) = IM(t) * zr * r;
            }
        }
    }
    return RES_OK;
}