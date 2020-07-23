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
#include "dspl_internal.h"




#ifdef DOXYGEN_ENGLISH

#endif
#ifdef DOXYGEN_RUSSIAN
/*!*****************************************************************************
\ingroup RESAMPLING_GROUP
\fn  int farrow_lagrange(double *s, int n, double p, double q, 
                         double frd, double **y, int *ny)
\brief Передискретизация вещественного сигнала на основе 
полиномиальной Лагранжевой интерполяции.

Данная функция осуществляет передискретизацию входного сигнала `s` в `p/q` раз 
со смещением дробной задержки `frd`. \n

Для передискретизации используется 
<a href = "http://ru.dsplib.org/content/resampling_lagrange/resampling_lagrange.html">
полиномиальная Лагранжева интерполяция
</a> (структура Фарроу для полиномиальной интерполяции). \n

\param [in] s
Указатель на вектор входного вещественного сигнала. \n
Размер вектора `[n x 1]`. \n \n

\param [in] n
Размер вектора входного сигнала. \n \n

\param [in] p
Числитель коэффициента передискретизации. \n \n

\param [in] q
Знаменатель коэффициента передискретизации. \n\n

\param [in] frd 
Значение смещения дробной задержки в пределах одного отсчета. \n
Значение должно быть от 0 до 1. \n 
\n
        
\param [out] y 
Указатель на адрес результата передискретизации. \n
По данному адресу будет произведено динамическое выделение памяти 
для результата передискретизации. \n
Будет выделено памяти под `n*q/p` отсчетов выходного сигнала. \n
Данный указатель не может быть `NULL`. \n\n

\param [in] ny
Указатель на переменную, в которую будет записан 
размер вектора `(*y)` после выделения памяти. \n \n

\return
`RES_OK` --- передискретизация рассчитана успешно.  \n
В противном случае \ref ERROR_CODE_GROUP "код ошибки".  \n

\author Бахурин Сергей. www.dsplib.org
***************************************************************************** */
#endif
int DSPL_API farrow_lagrange(double *s, int n, double p, double q,
                             double frd, double **y, int *ny)
{
    double a[4];
    double t, x, dt;
    int ind, k, res;
    double g[4];
    double *z;

    if(!s || !y)
        return ERROR_PTR;

    if(n<1)
        return ERROR_SIZE;

    if(p <= 0.0 || q <= 0.0)
        return ERROR_RESAMPLE_RATIO;

    if(frd <= -1.0 || frd >= 1.0)
        return ERROR_RESAMPLE_FRAC_DELAY;

    dt = q/p;

    if((*ny) != (int)((double)(n-1)/dt)+1 || !(*y))
    {

        *ny = (int)((double)(n-1)/dt)+1;
        (*y) = (double*)realloc((*y), (*ny)*sizeof(double));
    }

    t = -frd;
    k = 0;
    while(k < (*ny))
    {
        ind = (int)floor(t)+1;
        x = t - (double)ind;
        ind-=2;
        if(ind < 0)
        {
            memset(g, 0, 4*sizeof(double));
            if(ind > (-3))
                memcpy(g-ind, s, (4+ind)*sizeof(double));
            z = g;
        }
        else
        {
            if(ind < n-3)
                z = s+ind;
            else
            {
                memset(g, 0, 4*sizeof(double));
                if((n-ind)>0)
                    memcpy(g, s+ind, (n-ind)*sizeof(double));
                z = g;
            }
        }
        a[0] = z[2];
        a[3] = DSPL_FARROW_LAGRANGE_COEFF*(z[3] -z[0]) + 0.5*(z[1] - z[2]);
        a[1] = 0.5*(z[3] - z[1])-a[3];
        a[2] = z[3] - z[2] -a[3]-a[1];

        res = polyval(a, 3, &x, 1, (*y)+k);

        if(res != RES_OK)
            goto exit_label;
        t+=dt;
        k++;
    }

exit_label:
    return res;
}





#ifdef DOXYGEN_ENGLISH

#endif
#ifdef DOXYGEN_RUSSIAN
/*! ****************************************************************************
\ingroup RESAMPLING_GROUP
\fn  int int farrow_spline(double *s, int n, double p, double q, double frd, 
                           double **y, int *ny)
\brief Передискретизация вещественного сигнала на основе сплайн интерполяции.

Данная функция осуществляет передискретизацию 
входного сигнала `s` в `p/q` раз со смещением дробной задержки `frd`.  \n
Для передискретизации используются 
<a href = "http://ru.dsplib.org/content/resampling_spline/resampling_spline.html">
кубические сплайны Эрмита 
</a> 
(структура Фарроу для для сплайн-интерполяции). \n


\param [in] s
Указатель на вектор входного вещественного сигнала. \n
Размер вектора `[n x 1]`. \n \n
        
\param [in] n
Размер вектора входного сигнала. \n \n

\param [in] p
Числитель коэффициента передискретизации. \n \n

\param [in] q
Знаменатель коэффициента передискретизации. \n \n

\param [in] frd
Значение смещения дробной задержки в пределах одного отсчета. \n
Значение должно быть от 0 до 1. \n\n

\param [out] y
Указатель на адрес результата передискретизации. \n
По данному адресу будет произведено динамическое выделение памяти 
для результата передискретизации. \n
Будет выделено памяти под `n*q/p` отсчетов выходного сигнала. \n
Данный указатель не может быть `NULL`. \n \n

\param [in] ny
Указатель на переменную, в которую будет записан 
размер вектора `(*y)` после выделения памяти. \n \n

\return
`RES_OK` --- передискретизация рассчитана успешно.  \n
В противном случае \ref ERROR_CODE_GROUP "код ошибки".  \n        

\author Бахурин Сергей. www.dsplib.org
***************************************************************************** */
#endif
int DSPL_API farrow_spline(double *s, int n, double p, double q,
                           double frd, double **y, int *ny)
{
    double a[4];
    double t, x, dt;
    int ind, k, res;
    double g[4];
    double *z;

    if(!s || !y)
        return ERROR_PTR;

    if(n<1)
        return ERROR_SIZE;

    if(p <= 0.0 || q <= 0.0)
        return ERROR_RESAMPLE_RATIO;

    if(frd <= -1.0 || frd >= 1.0)
        return ERROR_RESAMPLE_FRAC_DELAY;

    dt = q/p;

    if((*ny) != (int)((double)(n-1)/dt)+1 || !(*y))
    {

        *ny = (int)((double)(n-1)/dt)+1;
        (*y) = (double*)realloc((*y), (*ny)*sizeof(double));
    }

    t = -frd;
    k = 0;
    while(k < (*ny))
    {
        ind = (int)floor(t)+1;
        x = t - (double)ind;
        ind-=2;
        if(ind < 0)
        {
            memset(g, 0, 4*sizeof(double));
            if(ind > (-3))
                memcpy(g-ind, s, (4+ind)*sizeof(double));
            z = g;
        }
        else
        {
            if(ind < n-3)
                z = s+ind;
            else
            {
                memset(g, 0, 4*sizeof(double));
                if((n-ind)>0)
                    memcpy(g, s+ind, (n-ind)*sizeof(double));
                z = g;
            }
        }
        a[0] = z[2];
        a[1] = 0.5*(z[3] - z[1]);
        a[3] = 2.0*(z[1] - z[2]) + a[1] + 0.5*(z[2] - z[0]);
        a[2] = z[1] - z[2] +a[3] + a[1];

        res = polyval(a, 3, &x, 1, (*y)+k);

        if(res != RES_OK)
            goto exit_label;
        t+=dt;
        k++;
    }

exit_label:
    return res;
}

