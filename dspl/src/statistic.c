/*
* Copyright (c) 2015-2019 Sergey Bakhurin
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
#include <string.h>
#include "dspl.h"





#ifdef DOXYGEN_ENGLISH
/*! ****************************************************************************
\ingroup SPEC_MATH_STAT_GROUP
\fn int find_max_abs(double* a, int n, double* m, int* ind)
\brief Find maximum absolute value from the real vector `a`

Function searches maximum absolute value in the real vector `a`.
This value writes to the address `m` and index keeps to te address `ind`.

\param[in] a
Pointer to the real vector `a`. \n
Vector size is `[n x 1]`. \n \n

\param[in] n
Size of the input vector `a`. \n \n


\param[out] m
Pointer to the variable which keeps vector `a`
maximum absolute value. \n
Pointer can be `NULL`, maximum value will not return
in this case. \n \n

\param[out] ind
Pointer to the variable which keeps index of a 
maximum absolute value inside vector `a`. \n
Pointer can be `NULL`, index will not return
in this case. \n \n

\return
`RES_OK` if function calculates successfully,
 else \ref ERROR_CODE_GROUP "code error".

Example:
\code{.cpp}
    double a[5] = {0.0, 2.0, -5.0, 4.0, 2.0};
    double m;
    int ind;
    find_max_abs(a, 5, &m, &ind);
    printf("\n\nmax absolute value:        %8.1f    (index %d)", m, ind);
\endcode 
As result the variable `m` will keep value `5`, 
and variable `ind` will keep `2`.

\author Sergey Bakhurin www.dsplib.org
***************************************************************************** */
#endif
#ifdef DOXYGEN_RUSSIAN
/*! ****************************************************************************
\ingroup SPEC_MATH_STAT_GROUP
\fn int find_max_abs(double* a, int n, double* m, int* ind)
\brief Поиск максимального по модулю элемента вещественного вектора `a`

Функция производит поиск максимального по модулю значения вектора `a`. \n
Максимальное по модулю значение `max|a[k]|` сохраняется по адресу `m`, а индекс
данного значения в векторе `a` сохраняется по адресу `ind`. \n

\param[in]  a
Указатель на вещественный вектор `a`. \n
Размер вектора `[n x 1]`. \n \n

\param[in]  n
Размер входного вектора `a`. \n \n

\param[out] m
Указатель на адрес памяти, в который сохранить
максимальное по модулю значение вектора `a`. \n
Указатель может быть `NULL`, в этом случае максимальное по модулю значение 
не сохраняется. \n \n

\param[out] ind Указатель на переменную, в которую будет сохранен
индекс максимального по модулю значению вектора `a`. \n
Указатель может быть `NULL`, в этом случае индекс не возвращается. \n \n

\return
`RES_OK` если функция выполнена успешно. \n
 В противном случае \ref ERROR_CODE_GROUP "код ошибки".

Пример:
\code{.cpp}
    double a[5] = {0.0, 2.0, -5.0, 4.0, 2.0};
    double m;
    int ind;
    find_max_abs(a, 5, &m, &ind);
    printf("\n\nmax absolute value:    %8.1f  (index %d)", m, ind);
\endcode 
В результате в переменную `m` будет записано значение `5`, 
а в переменную `ind` значение `2`.

\author Бахурин Сергей. www.dsplib.org
***************************************************************************** */
#endif
int DSPL_API find_max_abs(double* a, int n, double* m, int* ind)
{
    int k, i;
    double t;
    if(!a)
        return ERROR_PTR;
    if(n < 1)
        return ERROR_SIZE;
    t = fabs(a[0]);
    i = 0;
    for(k = 1; k < n; k++)
    {
        if(fabs(a[k]) > t)
        {
            t = fabs(a[k]);
            i = k;
        }
    }
    if(m)
        *m = t;
    if(ind)
        *ind = i;
    return RES_OK;
}


#ifdef DOXYGEN_ENGLISH

#endif
#ifdef DOXYGEN_RUSSIAN

#endif
int DSPL_API histogram(double* x, int n, int nh, double* pedges, double* ph)
{
    double xmin, xmax;
    int k, ind;
    int res;
    
    if(!x || !pedges || !ph)
        return ERROR_PTR;

    if(n<1 || nh<1)
        return ERROR_SIZE;

    res = minmax(x, n, &xmin, &xmax);
    if(res != RES_OK)
        return res;
    
    res = linspace(xmin, xmax, nh+1, DSPL_SYMMETRIC, pedges);
    if(res != RES_OK)
        return res;
    
    memset(ph, 0, nh*sizeof(double));
    for(k = 0; k < n; k++)
    {
        ind = 0;
        while(ind<nh && x[k]>=pedges[ind])
            ind++;
        ph[ind-1]+=1.0;
    }
    return RES_OK;
}




#ifdef DOXYGEN_ENGLISH

#endif
#ifdef DOXYGEN_RUSSIAN

#endif
int DSPL_API histogram_norm(double* y, int n, int nh, double* x, double* w)
{
    double *pedges = NULL;
    int k, res;
    
    if(!y || !x || !w)
        return ERROR_PTR;

    if(n<1 || nh<1)
        return ERROR_SIZE;
    
    pedges = (double*)malloc((nh+1)*sizeof(double));
    
    res = histogram(y, n, nh, pedges, w);
    if(res != RES_OK)
        goto exit_label;

    for(k = 1; k < nh+1; k++)
    {
        x[k-1] = 0.5*(pedges[k] + pedges[k-1]);
        w[k-1] /= ((double)n * (pedges[k] - pedges[k-1]));
    }
    
    res = RES_OK;    
exit_label: 

    if(pedges)
        free(pedges);
    return res;

}




#ifdef DOXYGEN_ENGLISH

#endif
#ifdef DOXYGEN_RUSSIAN

#endif
int DSPL_API minmax(double* x, int n, double* xmin, double* xmax)
{
    int k; 
    double min, max;
    
    if(!x)
        return ERROR_PTR;

    if(n<1)
        return ERROR_SIZE;

    min = max = x[0];
    for(k = 1; k < n; k++)
    {
        min = x[k] < min ? x[k] : min;
        max = x[k] > max ? x[k] : max;
    }
    
    if(xmin)
        *xmin = min;
    if(xmax)
        *xmax = max;
    
    return RES_OK;
}