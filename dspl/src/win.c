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
#include <math.h>
#include "dspl.h"
#include "dspl_internal.h"



#ifdef DOXYGEN_ENGLISH
/*! ****************************************************************************
\ingroup WIN_GROUP
\fn int window(double* w, int n, int win_type, double param)
\brief Window function calculation

The function calculates a periodic or symmetric window function
according to parameter `win_type`. \n

A periodic window function is used for spectral analysis, 
and a symmetric window function can be used to design FIR filters.

\param [in,out] w
Pointer to the window. \n
Vector size is `[n x 1]`. \n
Memory must be allocated. \n
The calculated window function will be placed at the given address. \n \n

\param [in]  n
Size of window function `w` vector. \n \n

\param [in]  win_type
Combination of flags for specifying the type of window function. \n
Combination of `DSPL_WIN_MASK | DSPL_WIN_SYM_MASK` bit masks 
is used to set the window type.\n
Bit mask `DSPL_WIN_MASK` sets the window type.
Can be one of follow: \n
\verbatim
-------------------------------------------------------------------------
 win_type                    |  Description
-----------------------------|-------------------------------------------
 DSPL_WIN_BARTLETT           | Nonparametric Bartlett window
-----------------------------|-------------------------------------------
 DSPL_WIN_BARTLETT_HANN      | Nonparametric Bartlett-Hann window
-----------------------------|-------------------------------------------
 DSPL_WIN_BLACKMAN           | Nonparametric  Blackman window 
-----------------------------|-------------------------------------------
 DSPL_WIN_BLACKMAN_HARRIS    | Nonparametric Blackman-Harris window
-----------------------------|-------------------------------------------
 DSPL_WIN_BLACKMAN_NUTTALL   | Nonparametric Blackman-Nuttall
-----------------------------|-------------------------------------------
 DSPL_WIN_CHEBY              | Parametric Dolph-Chebyshev window.
                             | Parametr  `win_param` sets sidelobe attenuation 
                             | level in dB.
-----------------------------|-------------------------------------------
 DSPL_WIN_COS                | Nonparametric Cosine window
-----------------------------|-------------------------------------------
 DSPL_WIN_FLAT_TOP           | Nonparametric maxflat window
-----------------------------|-------------------------------------------
 DSPL_WIN_GAUSSIAN           | Nonparametric Gauss window
-----------------------------|-------------------------------------------
 DSPL_WIN_HAMMING            | Nonparametric Hamming window
-----------------------------|-------------------------------------------
 DSPL_WIN_HANN               | Nonparametric Hann window
-----------------------------|-------------------------------------------
 DSPL_WIN_KAISER             | Parametric Kaiser window
-----------------------------|-------------------------------------------
 DSPL_WIN_LANCZOS            | Nonparametric Lanczos window
-----------------------------|-------------------------------------------
 DSPL_WIN_NUTTALL            | Nonparametric Nuttall window
-----------------------------|-------------------------------------------
 DSPL_WIN_RECT               | Nonparametric rectangular window
-------------------------------------------------------------------------
\endverbatim
\n 
Bit mask `DSPL_WIN_SYM_MASK` sets window function symmetry: \n
\verbatim
-------------------------------------------------------------------------
 DSPL_WIN_SYM_MASK           | Description
-----------------------------|-------------------------------------------
 DSPL_WIN_SYMMETRIC          | Symmetry window (default value)
 DSPL_WIN_PERIODIC           | Periodic window
-------------------------------------------------------------------------
\endverbatim
 \n \n

\param [in]  param
Window function parameter. \n
This parameter is using only to parametric window functions, 
and ignored for nonparametric windows. \n
\n

\return  
`RES_OK` if window function is calculated successfully. \n
Else \ref ERROR_CODE_GROUP "error code".

The following program calculates 64 samples window functions, 
draws their spectrum when using the bin indices of 
the discrete Fourier transform along the frequency axis.

\include windows_test.c

A personal graph is displayed for each type of window function.

Rectangular window
\image html win_rect.png
\n
\n

Nonparametric windows
\image html win_bartlett.png
\image html win_flattop.png
\image html win_bartletthann.png
\image html win_hann.png
\image html win_hamming.png
\image html win_lanczos.png
\image html win_blackman.png
\image html win_blackmanharris.png
\image html win_blackmannuttall.png
\image html win_cos.png
\image html win_nuttall.png
\n
\n

Parametric Dolph-Chebyshev windows
\image html win_cheby50.png
\image html win_cheby80.png
\image html win_cheby120.png
\n
\n

Parametric Gaussian windows
\image html win_gaussian0p5.png
\image html win_gaussian0p3.png

\n
\n
Parametric Kaiser windows
\image html win_kaiser4p0.png
\image html win_kaiser8p0.png
\image html win_kaiser12p0.png
\n
\n

\author Sergey Bakhurin. www.dsplib.org
***************************************************************************** */
#endif
#ifdef DOXYGEN_RUSSIAN
/*! ****************************************************************************
\ingroup WIN_GROUP
\fn int window(double* w, int n, int win_type, double param)
\brief Расчет функции оконного взвешивания

Функция рассчитывает периодическую или симметричную оконную функцию
в соответствии с параметром `win_type`. \n

Периодическая оконная функция используется для спектрального анализа,
а симметричная оконная функция может быть использована для синтеза 
КИХ-фильтров.

\param [in,out] w
Указатель на вектор оконной функции. \n
Размер вектора `[n x 1]`. \n
Память должна быть выделена. \n
Рассчитанная оконная функция будет помещена по данному адресу. \n \n

\param [in]  n
Размер вектора `w` оконной функции. \n \n

\param [in]  win_type
Комбинация флагов для задания типа оконной функции. \n
Для задания типа окна используется комбинация битовых масок
`DSPL_WIN_MASK | DSPL_WIN_SYM_MASK`. \n
Маска `DSPL_WIN_MASK` задает тип оконной функции.
Может принимать следующие значения: \n
\verbatim
-------------------------------------------------------------------------
Значение  DSPL_WIN_MASK      |  Описание
-----------------------------|-------------------------------------------
 DSPL_WIN_BARTLETT           | Непараметрическое окно Бартлетта
-----------------------------|-------------------------------------------
 DSPL_WIN_BARTLETT_HANN      | Непараметрическое окно Бартлетта-Ханна
-----------------------------|-------------------------------------------
 DSPL_WIN_BLACKMAN           | Непараметрическое окно Блэкмана 
-----------------------------|-------------------------------------------
 DSPL_WIN_BLACKMAN_HARRIS    | Непараметрическое окно Блэкмана-Харриса
-----------------------------|-------------------------------------------
 DSPL_WIN_BLACKMAN_NUTTALL   | Непараметрическое окно Блэкмана-Натталла
-----------------------------|-------------------------------------------
 DSPL_WIN_CHEBY              | Параметрическое окно Дольф-Чебышева. 
                             | Данное окно всегда является симметричным и
                             | игнорирует параметр  DSPL_WIN_SYM_MASK . 
                             | Параметр  param  задает уровень боковых 
                             | лепестков в дБ.
-----------------------------|-------------------------------------------
 DSPL_WIN_COS                | Непараметрическое косинус-окно
-----------------------------|-------------------------------------------
 DSPL_WIN_FLAT_TOP           | Непараметрическое окно с максимально 
                             | плоской вершиной
-----------------------------|-------------------------------------------
 DSPL_WIN_GAUSSIAN           | Параметрическое окно Гаусса
-----------------------------|-------------------------------------------
 DSPL_WIN_HAMMING            | Непараметрическое окно Хемминга
-----------------------------|-------------------------------------------
 DSPL_WIN_HANN               | Непараметрическое окно Ханна
-----------------------------|-------------------------------------------
 DSPL_WIN_KAISER             | Параметрическое окно Кайзера
-----------------------------|-------------------------------------------
 DSPL_WIN_LANCZOS            | Непараметрическое окно Ланкзоса
-----------------------------|-------------------------------------------
 DSPL_WIN_NUTTALL            | Непараметрическое окно Натталла
-----------------------------|-------------------------------------------
 DSPL_WIN_RECT               | Непараметрическое прямоугольное окно
-------------------------------------------------------------------------
\endverbatim
\n 
Маска `DSPL_WIN_SYM_MASK` задает симметричное 
или периодическое окно: \n
\verbatim
-------------------------------------------------------------------------
Значение  DSPL_WIN_SYM_MASK  | Описание
-----------------------------|-------------------------------------------
 DSPL_WIN_SYMMETRIC          | Симметричное окно (по умолчанию)
 DSPL_WIN_PERIODIC           | Периодическое окно
-------------------------------------------------------------------------
\endverbatim
 \n \n

\param [in]  param
Параметр окна. \n
Данный параметр применяется только для параметрических оконных функций. \n
Для непараметрических окон игнорируется. \n\n

\return  
`RES_OK` если оконная функция рассчитана успешно. \n
В противном случае \ref ERROR_CODE_GROUP "код ошибки".

Следующая программа производит расчет оконных функций длительности 64 отсчета, 
строит их спектральную плотность при использовании по оси частот индексы бинов 
дискретного преобразования Фурье.

\include windows_test.c

Для каждого вида оконной функция выводится персональный график.

Прямоугольное окно
\image html win_rect.png
\n
\n

Непраметрические окна
\image html win_bartlett.png
\image html win_flattop.png
\image html win_bartletthann.png
\image html win_hann.png
\image html win_hamming.png
\image html win_lanczos.png
\image html win_blackman.png
\image html win_blackmanharris.png
\image html win_blackmannuttall.png
\image html win_cos.png
\image html win_nuttall.png
\n
\n

Параметрические окна Дольф-Чебышева
\image html win_cheby50.png
\image html win_cheby80.png
\image html win_cheby120.png
\n
\n

Параметрические окна Гаусса
\image html win_gaussian0p5.png
\image html win_gaussian0p3.png

\n
\n
Параметрические окна Кайзера
\image html win_kaiser4p0.png
\image html win_kaiser8p0.png
\image html win_kaiser12p0.png
\n
\n

\author Бахурин Сергей. www.dsplib.org
***************************************************************************** */
#endif
int DSPL_API window(double* w, int n, int win_type, double param)
{
    switch(win_type & DSPL_WIN_MASK)
    {
        case DSPL_WIN_BARTLETT:
            return win_bartlett(w, n, win_type);
        case DSPL_WIN_BARTLETT_HANN:
            return win_bartlett_hann(w, n, win_type);
        case DSPL_WIN_BLACKMAN:
            return win_blackman(w, n, win_type);
        case DSPL_WIN_BLACKMAN_HARRIS:
            return win_blackman_harris(w, n, win_type);
        case DSPL_WIN_BLACKMAN_NUTTALL:
            return win_blackman_nuttall(w, n, win_type);
        case DSPL_WIN_CHEBY:
            return win_cheby(w, n, param);
        case DSPL_WIN_FLAT_TOP:
            return win_flat_top(w, n, win_type);
        case DSPL_WIN_GAUSSIAN:
            return win_gaussian(w, n, win_type, param);
        case DSPL_WIN_HAMMING:
            return win_hamming(w, n, win_type);
        case DSPL_WIN_HANN:
            return win_hann(w, n, win_type);
        case DSPL_WIN_KAISER:
            return win_kaiser(w, n, win_type, param);
        case DSPL_WIN_LANCZOS:
            return win_lanczos(w, n, win_type);
        case DSPL_WIN_NUTTALL:
            return win_nuttall(w, n, win_type);
        case DSPL_WIN_RECT:
            return win_rect(w, n);
        case DSPL_WIN_COS:
            return win_cos(w, n, win_type);
        default:
            return ERROR_WIN_TYPE;
    }
    return RES_OK;
}



/******************************************************************************
Barlett window function
*******************************************************************************/
int win_bartlett(double *w, int n, int win_type)
{
    double x = 0.0;
    int i;

    if(!w)
        return ERROR_PTR;
    if(n<2)
        return ERROR_SIZE;

    switch(win_type & DSPL_WIN_SYM_MASK)
    {
        case DSPL_WIN_SYMMETRIC: x = (double)(n-1); break;
        case DSPL_WIN_PERIODIC : x = (double)n;     break;
        default: return ERROR_WIN_SYM;
    }

    for(i = 0; i < n; i++)
    {
        w[i] = 2.0 / x * (x * 0.5-fabs((double)i - x * 0.5));
    }
    return RES_OK;
}





/******************************************************************************
Barlett - Hann    window function
******************************************************************************/
int win_bartlett_hann(double *w, int n, int win_type)
{
    double y;
    double x = 0.0;
    int i;

    if(!w)
    return ERROR_PTR;
    if(n<2)
        return ERROR_SIZE;

    switch(win_type & DSPL_WIN_SYM_MASK)
    {
        case DSPL_WIN_SYMMETRIC: x = 1.0/(double)(n-1); break;
        case DSPL_WIN_PERIODIC : x = 1.0/(double)n; break;
        default: return ERROR_WIN_SYM;
    }

    y = 0.0;
    for(i = 0; i<n; i++)
    {
        w[i] = 0.62 - 0.48 * fabs(y-0.5)-0.38*cos(M_2PI*y);
        y += x;
    }
    return RES_OK;
}





/******************************************************************************
Blackman    window function
******************************************************************************/
int win_blackman(double *w, int n, int win_type)
{
    double y;
    double x = 0.0;
    int i;

    if(!w)
    return ERROR_PTR;
    if(n<2)
        return ERROR_SIZE;

    switch(win_type & DSPL_WIN_SYM_MASK)
    {
        case DSPL_WIN_SYMMETRIC: x = M_2PI/(double)(n-1); break;
        case DSPL_WIN_PERIODIC : x = M_2PI/(double)n; break;
        default: return ERROR_WIN_SYM;
    }

    y = 0.0;
    for(i = 0; i<n; i++)
    {
        w[i] = 0.42 - 0.5* cos(y)+0.08*cos(2.0*y);
        y += x;
    }
    return RES_OK;
}




/******************************************************************************
Blackman - Harris window function
******************************************************************************/
int win_blackman_harris(double *w, int n, int win_type)
{
    double y;
    double x  = 0.0;
    double a0 = 0.35875;
    double a1 = 0.48829;
    double a2 = 0.14128;
    double a3 = 0.01168;
    int i;

    if(!w)
        return ERROR_PTR;
    if(n<2)
        return ERROR_SIZE;

    switch(win_type & DSPL_WIN_SYM_MASK)
    {
        case DSPL_WIN_SYMMETRIC: x = M_2PI/(double)(n-1); break;
        case DSPL_WIN_PERIODIC : x = M_2PI/(double)n; break;
        default: return ERROR_WIN_SYM;
    }

    y = 0.0;
    for(i = 0; i<n; i++)
    {
        w[i] = a0 - a1* cos(y)+a2*cos(2.0*y)-a3*cos(3.0*y);
        y += x;
    }
    return RES_OK;
}




/******************************************************************************
Blackman - Nuttull     window function
******************************************************************************/
int win_blackman_nuttall(double *w, int n, int win_type)
{
    double y;
    double x    = 0.0;
    double a0 = 0.3635819;
    double a1 = 0.4891775;
    double a2 = 0.1365995;
    double a3 = 0.0106411;
    int i;


    if(!w)
        return ERROR_PTR;
    if(n<2)
        return ERROR_SIZE;

    switch(win_type & DSPL_WIN_SYM_MASK)
    {
        case DSPL_WIN_SYMMETRIC:        x = M_2PI/(double)(n-1);        break;
        case DSPL_WIN_PERIODIC :        x = M_2PI/(double)n;                break;
        default: return ERROR_WIN_SYM;
    }

    y = 0.0;
    for(i = 0; i<n; i++)
    {
        w[i] = a0 - a1* cos(y)+a2*cos(2.0*y)-a3*cos(3.0*y);
        y += x;
    }
    return RES_OK;
}






/******************************************************************************
Chebyshev parametric window function
param sets spectrum sidelobes level in dB
ATTENTION! ONLY SYMMETRIC WINDOW
*******************************************************************************/
int win_cheby(double *w, int n, double param)
{
    int k, i, m;
    double z, dz, sum = 0, wmax=0, r1, x0, chx, chy, in;

    if(!w)
        return ERROR_PTR;

    if(n<2)
        return ERROR_SIZE;

    if(param <= 0.0)
        return ERROR_WIN_PARAM;

    r1 = pow(10, param/20);
    x0 = cosh((1.0/(double)(n-1)) * acosh(r1));

    /* check window length even or odd */
    if(n%2==0)
    {
        dz = 0.5;
        m = n/2-1;
    }
    else
    {
        m = (n-1)/2;
        dz = 0.0;
    }

    for(k = 0; k < m+2; k++)
    {
        z = (double)(k - m) - dz;
        sum = 0;

        for(i = 1; i <= m; i++)
        {
            in = (double)i / (double)n;
            chx = x0 * cos(M_PI * in);
            cheby_poly1(&chx, 1, n-1, &chy);
            sum += chy * cos(2.0 * z * M_PI * in);
        }

        w[k] = r1 + 2.0 * sum;
        w[n-1-k] = w[k];

        /* max value calculation */
        if(w[k]>wmax)
            wmax=w[k];
    }

    /* normalization */
    for(k=0; k < n; k++)
        w[k] /= wmax;

    return RES_OK;
}



/******************************************************************************
Cosine window function
******************************************************************************/
int win_cos(double *w, int n, int win_type)
{
    double y;
    double x = 0.0;
    int i;

    if(!w)
        return ERROR_PTR;
    if(n<2)
        return ERROR_SIZE;

    switch(win_type & DSPL_WIN_SYM_MASK)
    {
        case DSPL_WIN_SYMMETRIC: x = M_PI/(double)(n-1); break;
        case DSPL_WIN_PERIODIC : x = M_PI/(double)n; break;
        default: return ERROR_WIN_SYM;
    }

    y = 0.0;
    for(i = 0; i<n; i++)
    {
        w[i] = sin(y);
        y += x;
    }
    return RES_OK;
}






/******************************************************************************
Flat - Top     window function
******************************************************************************/
int win_flat_top(double *w, int n, int win_type)
{
    double y;
    double x  = 0.0;
    double a0 = 1.0;
    double a1 = 1.93;
    double a2 = 1.29;
    double a3 = 0.388;
    double a4 = 0.032;
    int i;

    if(!w)
        return ERROR_PTR;
    if(n<2)
        return ERROR_SIZE;

    switch(win_type & DSPL_WIN_SYM_MASK)
    {
        case DSPL_WIN_SYMMETRIC: x = M_2PI/(double)(n-1); break;
        case DSPL_WIN_PERIODIC : x = M_2PI/(double)n; break;
        default: return ERROR_WIN_SYM;
    }

    y = 0.0;
    for(i = 0; i<n; i++)
    {
        w[i] = a0 - a1* cos(y)+a2*cos(2.0*y)-a3*cos(3.0*y)+a4*cos(4.0*y);
        y += x;
    }
    return RES_OK;
}






/******************************************************************************
Gaussian window function
******************************************************************************/
int win_gaussian(double *w, int n, int win_type, double alpha)
{
    double a = 0.0;
    double y;
    double sigma;
    int i;

    if(!w)
        return ERROR_PTR;
    if(n<2)
        return ERROR_SIZE;

    switch(win_type & DSPL_WIN_SYM_MASK)
    {
        case DSPL_WIN_SYMMETRIC: a = (double)(n-1)*0.5; break;
        case DSPL_WIN_PERIODIC : a = (double)(n)*0.5; break;
        default: return ERROR_WIN_SYM;
    }


    sigma = 1.0 / (alpha * a);
    for(i = 0; i<n; i++)
    {
        y = ((double)i - a)*sigma;
        w[i] = exp(-0.5*y*y);
    }
        return RES_OK;
}






/******************************************************************************
Hamming window function
******************************************************************************/
int win_hamming(double *w, int n, int win_type)
{
    double x = 0.0;
    double y;
    int i;

    if(!w)
        return ERROR_PTR;
    if(n<2)
        return ERROR_SIZE;

    switch(win_type & DSPL_WIN_SYM_MASK)
    {
        case DSPL_WIN_SYMMETRIC: x = M_2PI/(double)(n-1); break;
        case DSPL_WIN_PERIODIC : x = M_2PI/(double)n; break;
        default: return ERROR_WIN_SYM;
    }

    y = 0.0;
    for(i = 0; i < n; i++)
    {
        w[i] = 0.54-0.46*cos(y);
        y += x;
    }
    return RES_OK;
}




/******************************************************************************
Hann window function
******************************************************************************/
int win_hann(double *w, int n, int win_type)
{
    double x = 0.0;
    double y;
    int i;

    if(!w)
        return ERROR_PTR;
    if(n<2)
        return ERROR_SIZE;

    switch(win_type & DSPL_WIN_SYM_MASK)
    {
        case DSPL_WIN_SYMMETRIC: x = M_2PI/(double)(n-1); break;
        case DSPL_WIN_PERIODIC : x = M_2PI/(double)n; break;
        default: return ERROR_WIN_SYM;
    }

    y = 0.0;
    for(i = 0; i < n; i++)
    {
        w[i] = 0.5*(1-cos(y));
        y += x;
    }
    return RES_OK;
}


/******************************************************************************
Kaiser window function
******************************************************************************/
int win_kaiser(double* w, int n, int win_type, double param)
{
    double num, den, x, y, L;
    int i, err;
    if(!w)
        return ERROR_PTR;
    if(n<2)
        return ERROR_SIZE;

    switch(win_type & DSPL_WIN_SYM_MASK)
    {
        case DSPL_WIN_SYMMETRIC: L = (double)(n-1) / 2.0; break;
        case DSPL_WIN_PERIODIC : L = (double)n / 2.0; break;
        default: return ERROR_WIN_SYM;
    }

    err = bessel_i0(&param, 1, &den);
    if(err != RES_OK)
        return err;
    for(i = 0; i < n; i++)
    { 
        x = 2.0*((double)i - L) / (double)n;
        y = param * sqrt(1.0 - x*x);
        err = bessel_i0(&y, 1, &num);
        if(err != RES_OK)
            return err;
        w[i] = num / den;
    }
    return err;
}



/******************************************************************************
Lanczos window function
******************************************************************************/
int win_lanczos(double *w, int n, int win_type)
{
    double y;
    double x = 0.0;
    int i;

    if(!w)
        return ERROR_PTR;
    if(n<2)
        return ERROR_SIZE;

    switch(win_type & DSPL_WIN_SYM_MASK)
    {
        case DSPL_WIN_SYMMETRIC: x = M_2PI/(double)(n-1); break;
        case DSPL_WIN_PERIODIC : x = M_2PI/(double)n; break;
        default: return ERROR_WIN_SYM;
    }

    y = 0.0;
    for(i = 0; i < n; i++)
    {
        if((y - M_PI)==0.0)
            w[i] = 1.0;
        else
            w[i] = sin(y - M_PI)/(y - M_PI);
        y += x;
    }
    return RES_OK;

}



/******************************************************************************
Nuttall window function
******************************************************************************/
int win_nuttall(double *w, int n, int win_type)
{
    double y;
    double x  = 0.0;
    double a0 = 0.355768;
    double a1 = 0.487396;
    double a2 = 0.144232;
    double a3 = 0.012604;
    int i;

    if(!w)
        return ERROR_PTR;
    if(n<2)
        return ERROR_SIZE;

    switch(win_type & DSPL_WIN_SYM_MASK)
    {
        case DSPL_WIN_SYMMETRIC: x = M_2PI/(double)(n-1); break;
        case DSPL_WIN_PERIODIC : x = M_2PI/(double)n; break;
        default: return ERROR_WIN_SYM;
    }

    y = 0.0;
    for(i = 0; i < n; i++)
    {
        w[i] = a0 - a1* cos(y)+a2*cos(2.0*y)-a3*cos(3.0*y);
        y += x;
    }
    return RES_OK;
}





/******************************************************************************
Rectangle window function
******************************************************************************/
int win_rect(double *w, int n)
{
    int i;

    if(!w)
        return ERROR_PTR;
    if(n<2)
        return ERROR_SIZE;

    for(i = 0; i < n; i++)
        w[i] = 1.0;
    return RES_OK;
}

