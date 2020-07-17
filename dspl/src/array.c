/*
* Copyright (c) 2015-2020 Sergey Bakhurin
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
#include "blas.h"


#ifdef DOXYGEN_ENGLISH
/*! ****************************************************************************
\ingroup ARRAY_GROUP
\fn int array_scale_lin(double* x, int n,
                        double xmin, double xmax, double dx,
                        double h, double* y)
\brief Vector `x` linear transformation

Function transforms values \f$x(i)\f$, \f$i = 0,1,\ldots n\f$
to the \f$y(i)\f$, accordint to equation:

\f[
y(i) = k_x x(i) + d_x, \qquad k_x =
\frac{h}{x_{\textrm{max}} - x_{\textrm{min}}}.
\f]

All values of the vector `x` between
\f$x_{\textrm{min}}\f$ and \f$x_{\textrm{max}}\f$, transforms to
the vector `y` between \f$d_x\f$ and \f$h + d_x\f$.
Parameter \f$d_x\f$ sets mean shift of the vector `y`.

This function is convenient for translating values ​​
of different dimensions. For example it can be used
to transfer the values ​​of the vector `x`
to the graph of the height of` h`, where the height can
be set in the number of pixels, in centimeters, etc.

\param[in] x
Pointer to the input vector `x`. \n
Vector size is `[n x 1]`. \n
\n

\param[in] n
Size of vector `x`. \n
\n

\param[in] xmin
Parameter \f$x_{\textrm{min}}\f$. \n
\n

\param[in] xmax
Parameter \f$x_{\textrm{min}}\f$. \n
Value `xmax` must be more than `xmin`. \n
\n

\param[in] dx
Displacement after transformation. \n
This parameter must have output vector `y`
dimensions (pixels, centimeters). \n
\n

\param[in] h
Height of vector `y` after transforming between `dx` and `h+dx`. \n
\n

\param[out] y
Pointer to the output vector `y`. \n
Vector size is `[n x 1]`. \n
Memory must be allocated. \n
\note
Pointer    `y` can be equal to `x`.
Velues of vector `x` will be rewritten in this case. \n
\n

\return
`RES_OK` if function returns successfully. \n
Else \ref ERROR_CODE_GROUP "code error".

\author Sergey Bakhurin www.dsplib.org
***************************************************************************** */
#endif

#ifdef DOXYGEN_RUSSIAN
/*! ****************************************************************************
\ingroup ARRAY_GROUP
\fn int array_scale_lin(double* x, int n,
                        double xmin, double xmax, double dx,
                        double h, double* y)
\brief Линейное растяжение вектора данных `x`
Функция производит преобразование значений \f$x(i)\f$, \f$i = 0,1,\ldots n\f$
в значения \f$y(i)\f$, в соответствии с формулой:

\f[
y(i) = k_x x(i) + d_x, \qquad k_x =
\frac{h}{x_{\textrm{max}} - x_{\textrm{min}}}.
\f]

Таким образом, все значения входного вектора `x` в диапазоне от
\f$x_{\textrm{min}}\f$ до \f$x_{\textrm{max}}\f$, линейно растягиваются в
значения вектора `y` в диапазоне от \f$d_x\f$ до \f$h + d_x\f$.
Заметим, что \f$d_x\f$ задает линейное смещение значений вектора `y`.

Данная функция удобна для перевода величин разных размерностей, в частности,
для переноса значений вектора `x` на график высоты `h`, где высота может
быть задана в количестве пикселей, в сантиметрах и т.д.

\param[in] x
Указатель на вектор входных значений `x`. \n
Размер вектора `[n x 1]`. \n
\n

\param[in] n
Размер вектора `x`. \n
\n

\param[in] xmin
Нижняя граница диапазона трансформации. \n
\n

\param[in] xmax
Верхняя граница диапазона трансформации. \n
Значение `xmax` должно быть строго больше значения `xmin`. \n
\n

\param[in] dx
Смещение после трансформации. \n
Данный параметр должен иметь размерность выходного вектора `y`. \n
\n

\param[in] h
Диапазон значений вектора `y` после трансформации от `dx` до `h+dx`. \n
\n

\param[out] y
Указатель на вектора данных после трансформации. \n
Размер вектора `[n x 1]`. \n
Память должна быть выделена. \n
\note
Указатель `y` может совпадать с `x`, в этом случае,
данные вектора `x` будут перезаписаны линейно измененными в соответствии
с формулой выше. \n
\n

\return
`RES_OK` если функция выполнена успешно. \n
 В противном случае \ref ERROR_CODE_GROUP "код ошибки".

\author
Бахурин Сергей
www.dsplib.org
***************************************************************************** */
#endif
int DSPL_API array_scale_lin(double* x,   int n,
                             double xmin, double xmax, double dx,
                             double h,    double* y)
{
    double kx;
    int k;
    if(!x)
        return ERROR_PTR;
    if(n < 1)
        return ERROR_SIZE;
    if(h<0.0)
        return ERROR_NEGATIVE;

    if(xmin >= xmax)
        return ERROR_MIN_MAX;

    kx = h / (xmax - xmin);

    for(k = 0; k < n; k++)
        y[k] = (x[k] - xmin) * kx + dx;

    return RES_OK;
}



#ifdef DOXYGEN_ENGLISH
/*! ****************************************************************************
\ingroup ARRAY_GROUP
\fn int concat(void* a, size_t na, void* b, size_t nb, void* c)
\brief
Concatenate arrays `a` and `b`

Let's arrays `a` and `b` are vectors: \n
`a = [a(0), a(1), ... a(na-1)]`, \n
`b = [b(0), b(1), ... b(nb-1)]`, \n
concatenation of these arrays will be array `c` size `na+nb`: \n
`c = [a(0), a(1), ... a(na-1), b(0), b(1), ... b(nb-1)]`.


\param[in] a
Pointer to the first array `a`. \n
Array `a` size is `na` bytes. \n
\n

\param[in] na
Array `a` size (bytes). \n
\n

\param[in] b
Pointer to the second array `b`. \n
Array `b` size is `nb` bytes. \n
\n

\param[in] nb
Array `a` size (bytes). \n
\n

\param[out] c
Pointer to the concatenation result array `c`. \n
Array `c` size is `na + nb` bytes. \n
Memory must be allocated. \n
\n

\return
`RES_OK` if function returns successfully. \n
Else \ref ERROR_CODE_GROUP "code error".

Function uses pointer type `void*` and can be useful for an arrays
concatenation with different types. \n
For example two `double` arrays concatenation:
\code{.cpp}
double a[3] = {1.0, 2.0, 3.0};
double b[2] = {4.0, 5.0};
double c[5];

concat((void*)a, 3*sizeof(double), (void*)b, 2*sizeof(double), (void*)c);
\endcode
Vector `c` keeps follow data:
\verbatim
c = [1.0, 2.0, 3.0, 4.0, 5.0]
\endverbatim

\author Sergey Bakhurin www.dsplib.org
***************************************************************************** */
#endif

#ifdef DOXYGEN_RUSSIAN
/*! ****************************************************************************
\ingroup ARRAY_GROUP
\fn int concat(void* a, size_t na, void* b, size_t nb, void* c)
\brief Конкатенация двух массивов данных

Функция производит конкатенацию двух массивов. Пусть массивы `a` и `b`
заданы как векторы: \n
`a = [a(0), a(1), ... a(na-1)]`,    \n
`b = [b(0), b(1), ... b(nb-1)]`,    \n
тогда результатом конкатенации будет вектор размера `na+nb` вида: \n
`c = [a(0), a(1), ... a(na-1), b(0), b(1), ... b(nb-1)]`.


\param[in] a
Указатель на первый вектор `a`. \n
Размер вектора `na` байт. \n \n

\param[in] na
Размер первого вектора `a` в байт. \n \n

\param[in] b
Указатель на второй вектор `b`. \n
Размер памяти вектора `nb` байт. \n \n

\param[in] nb
Размер второго вектора `b` в байт. \n \n

\param[out] c
Указатель на вектор конкатенации `c`. \n
Размер памяти вектора `na + nb` байт. \n
Память должна быть выделена. \n \n

\return
`RES_OK` если функция выполнена успешно. \n
 В противном случае \ref ERROR_CODE_GROUP "код ошибки".

\note
Функция использует указатели типа `void*` и может быть использована для
конкатенации данных различного типа. \n
Например конкатенация массивов типа `double`:
\code{.cpp}
double a[3] = {1.0, 2.0, 3.0};
double b[2] = {4.0, 5.0};
double c[5];
concat((void*)a, 3*sizeof(double), (void*)b, 2*sizeof(double), (void*)c);
\endcode
в результате вектор `c` будет хранить массив данных:
\verbatim
c = [1.0, 2.0, 3.0, 4.0, 5.0]
\endverbatim

\author
Бахурин Сергей
www.dsplib.org
***************************************************************************** */
#endif
int DSPL_API concat(void* a, size_t na, void* b, size_t nb, void* c)
{
    if(!a || !b || !c || c == b)
        return ERROR_PTR;
    if(na < 1 || nb < 1)
        return ERROR_SIZE;

    if(c != a)
        memcpy(c, a, na);

    memcpy((char*)c+na, b, nb);
    return RES_OK;
}





#ifdef DOXYGEN_ENGLISH
/*! ****************************************************************************
\ingroup ARRAY_GROUP
\fn int decimate(double* x, int n, int d, double* y, int* cnt)
\brief
Real vector decimation

Function `d` times decimates real vector `x`. \n
Output vector `y` keeps values corresponds to:
`y(k) = x(k*d), k = 0...n/d-1` \n

\param[in] x
Pointer to the input real vector `x`. \n
Vector `x` size is `[n x 1]`. \n \n

\param[in] n
Size of input vector `x`. \n \n

\param[in] d
Decimation coefficient. \n
Each d-th vector will be copy from vector `x` to the
output vector `y`. \n \n

\param[out] y
Pointer to the output decimated vector `y`. \n
Output vector size is `[n/d x 1]` will be copy
to the address `cnt`. \n

\param[out] cnt
Address which will keep decimated vector `y` size. \n
Pointer can be `NULL`, vector `y` will not return
in this case. \n \n

\return
`RES_OK` if function calculated successfully. \n
Else \ref ERROR_CODE_GROUP "code error".

Two-times decimation example:
\code{.cpp}
double x[10] = {0.0, 1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0, 9.0};
double y[5];
int d = 2;
int cnt;

decimate(x, 10, d, y, &cnt);
\endcode
As result variable `cnt` will be written value 5 and
vector `y` will keep    array:
\verbatim
c = [0.0, 2.0, 4.0, 6.0, 8.0]
\endverbatim

\author Sergey Bakhurin www.dsplib.org
***************************************************************************** */
#endif
#ifdef DOXYGEN_RUSSIAN
/*! ****************************************************************************
\ingroup ARRAY_GROUP
\fn int decimate(double* x, int n, int d, double* y, int* cnt)
\brief Децимация вещественного вектора данных

Функция производит децимацию вещественного вектора `x` в `d` раз. \n
В результате выходной вектор `y` содержит значения:
`y(k) = x(k*d), k = 0...n/d-1` \n

\param[in] x
Указатель на вектор входных данных `x`. \n
Размер вектора `[n x 1]`. \n \n

\param[in] n
Размер входного вектора `x`. \n \n

\param[in] d
Коэффициент децимации. \n
В результате децимации из вектора `x` будет взять каждый
d-й элемент. \n \n

\param[out] y
Указатель на децимированный вектор `y`. \n
Размер выходного вектора равен `[n/d x 1]`
будет сохранен по адресу `cnt`. \n
Память должна быть выделена. \n \n

\param[out] cnt
Указатель переменную, в которую будет сохранен
размер выходного вектора после децимации. \n
Указатель может быть `NULL`, в этом случае
размер вектора `y` не возвращается. \n \n

\return
`RES_OK` если функция выполнена успешно. \n
 В противном случае \ref ERROR_CODE_GROUP "код ошибки".

Пример децимации вещественного массива данных в 2 раза:
\code{.cpp}
double x[10] = {0.0, 1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0, 9.0};
double y[5];
int d = 2;
int cnt;
decimate(x, 10, d, y, &cnt);
\endcode
В результате в переменную `cnt` будет записан размер 5,
а вектор `y` будет хранить массив данных:
\verbatim
c = [0.0, 2.0, 4.0, 6.0, 8.0]
\endverbatim

\author
Бахурин Сергей
www.dsplib.org
**************************************************************************** */
#endif
int DSPL_API decimate(double* x, int n, int d, double* y, int* cnt)
{
    int k = 0, i = 0;
    if(!x || !y)
        return ERROR_PTR;
    if(n < 1)
        return ERROR_SIZE;
    if(d < 1)
        return ERROR_NEGATIVE;

    k = i = 0;
    while(k + d <= n)
    {
        y[i] = x[k];
        k+=d;
        i++;
    }
    if(cnt)
        *cnt = i;

    return RES_OK;
}




#ifdef DOXYGEN_ENGLISH
/*! ****************************************************************************
\ingroup ARRAY_GROUP
\fn int decimate_cmplx(complex_t* x, int n, int d, complex_t* y, int* cnt)
\brief
Complex vector decimation

Function `d` times decimates a complex vector `x`. \n
Output vector `y` keeps values corresponds to:
`y(k) = x(k*d), k = 0...n/d-1` \n

\param[in]    x
Pointer to the input complex vector `x`. \n
Vector `x` size is `[n x 1]`. \n \n

\param[in]    n
Size of input vector `x`. \n \n

\param[in]    d
Decimation coefficient. \n
Each d-th vector will be copy from vector `x` to the
output vector `y`. \n \n

\param[out] y
Pointer to the output decimated vector `y`. \n
Output vector size is `[n/d x 1]` will be copy
to the address `cnt`. \n
Memory must be allocated. \n \n

\param[out] cnt
Address which will keep decimated vector `y` size. \n
Pointer can be `NULL`, vector `y` will not return
in this case. \n \n

\return
`RES_OK` if function calculated successfully. \n
Else \ref ERROR_CODE_GROUP "code error".

Two-times complex vector decimation example:

\code{.cpp}
compex_t x[10] = {{0.0, 0.0}, {1.0, 1.0}, {2.0, 2.0}, {3.0, 3.0}, {4.0, 4.0},
{5.0, 5.0}, {6.0, 6.0}, {7.0, 7.0}, {8.0, 8.0}, {9.0, 9.0}};
compex_t y[5];
int d = 2;
int cnt;

decimate_cmplx(x, 10, d, y, &cnt);
\endcode

As result variable `cnt` will be written value 5 and
vector `y` will keep    array:

\verbatim
c = [0.0+0.0j, 2.0+2.0j, 4.0+4.0j, 6.0+6.0j, 8.0+8.0j]
\endverbatim

\author Sergey Bakhurin www.dsplib.org
***************************************************************************** */
#endif
#ifdef DOXYGEN_RUSSIAN
/*! ****************************************************************************
\ingroup ARRAY_GROUP
\fn int decimate_cmplx(complex_t* x, int n, int d, complex_t* y, int* cnt)
\brief Децимация комплексного вектора данных

Функция производит децимацию комплексного вектора `x` в `d` раз. \n
В результате выходной вектор `y` содержит значения:
`y(k) = x(k*d), k = 0...n/d-1` \n

\param[in] x
Указатель на вектор входных данных `x`. \n
Размер вектора `[n x 1]`. \n \n

\param[in] n
Размер входного вектора `x`. \n \n

\param[in] d
Коэффициент децимации. \n
В результате децимации из вектора `x` будет взять каждый d-й элемент. \n \n

\param[out] y
Указатель на децимированный вектор `y`. \n
Размер выходного вектора равен `[n/d x 1]` будет сохранен по адресу `cnt`. \n
Память должна быть выделена. \n \n

\param[out] cnt
Указатель переменную, в которую будет сохранен
размер выходного вектора после децимации. \n
Указатель может быть `NULL`, в этом случае
размер вектора `y` не возвращается. \n \n

\return
`RES_OK` если функция выполнена успешно. \n
 В противном случае \ref ERROR_CODE_GROUP "код ошибки".

Пример децимации комплексного массива данных в 2 раза:
\code{.cpp}
compex_t x[10] = {{0.0, 0.0}, {1.0, 1.0}, {2.0, 2.0}, {3.0, 3.0}, {4.0, 4.0},
                  {5.0, 5.0}, {6.0, 6.0}, {7.0, 7.0}, {8.0, 8.0}, {9.0, 9.0}};
compex_t y[5];
int d = 2;
int cnt;
decimate_cmplx(x, 10, d, y, &cnt);
\endcode
В результате в переменную `cnt` будет записан размер 5, а вектор `y` будет
хранить массив данных:
\verbatim
c = [0.0+0.0j, 2.0+2.0j, 4.0+4.0j, 6.0+6.0j, 8.0+8.0j]
\endverbatim

\author
Бахурин Сергей
www.dsplib.org
***************************************************************************** */
#endif
int DSPL_API decimate_cmplx(complex_t* x, int n, int d, complex_t* y, int* cnt)
{
    int k = 0, i = 0;
    if(!x || !y)
        return ERROR_PTR;
    if(n < 1)
        return ERROR_SIZE;
    if(d < 1)
        return ERROR_NEGATIVE;

    k = i = 0;
    while(k + d < n)
    {
        RE(y[i]) = RE(x[k]);
        IM(y[i]) = IM(x[k]);
        k+=d;
        i++;
    }
    if(cnt)
        *cnt = i;

    return RES_OK;
}





#ifdef DOXYGEN_ENGLISH
/*! ****************************************************************************
\ingroup ARRAY_GROUP
\fn int flipip(double* x, int n)
\brief
Flip real vector `x` in place

Function flips real vector `x` length `n` in the memory. \n
For example real vector `x` length 6:\n
\verbatim
x = [0, 1, 2, 3, 4, 5]
\endverbatim
After flipping it will be as follow:
\verbatim
x = [5, 4, 3, 2, 1, 0]
\endverbatim

\param[in, out] x
Pointer to the real vector `x`. \n
Vector size is `[n x 1]`. \n
Flipped vector will be on the same address. \n
\n

\param[in] n
Length of the vector `x`. \n
\n

\return
`RES_OK` if function returns successfully. \n
 Else \ref ERROR_CODE_GROUP "error code".

Example:
\code{.cpp}
double x[5] = {0.0, 1.0, 2.0, 3.0, 4.0};
int i;
for(i = 0; i < 5; i++)
    printf("%6.1f    ", x[i]);
flipip(x, 5);
printf("\n");
for(i = 0; i < 5; i++)
    printf("%6.1f    ", x[i]);
\endcode
\n
Program result:
\verbatim
     0.0         1.0         2.0         3.0         4.0
     4.0         3.0         2.0         1.0         0.0
\endverbatim

\author Sergey Bakhurin www.dsplib.org
***************************************************************************** */
#endif
#ifdef DOXYGEN_RUSSIAN
/*! ****************************************************************************
\ingroup ARRAY_GROUP
\fn int flipip(double* x, int n)
\brief Функция отражения вещественного вектора `x`

Функция производит отражение вещественного вектора длины `n`
в памяти данных. \n
Например исходный вектор `x`    длины 6: \n
\verbatim
x = [0, 1, 2, 3, 4, 5]
\endverbatim
После отражения вектор `x` будет иметь вид:
\verbatim
x = [5, 4, 3, 2, 1, 0]
\endverbatim

\param[in, out] x
Указатель на вещественный вектор `x`. \n
Размер вектора `[n x 1]`. \n
Результат отражения будет помещен по этому же адресу. \n

\param[in] n
Размер вектора `x`. \n \n

\return
`RES_OK` если функция выполнена успешно. \n
В противном случае \ref ERROR_CODE_GROUP "код ошибки".

Пример:
\code{.cpp}
double x[5] = {0.0, 1.0, 2.0, 3.0, 4.0};
int i;
for(i = 0; i < 5; i++)
    printf("%6.1f    ", x[i]);
flipip(x, 5);
printf("\n");
for(i = 0; i < 5; i++)
    printf("%6.1f    ", x[i]);
\endcode
\n
Результат выполнения:
\verbatim
     0.0         1.0         2.0         3.0         4.0
     4.0         3.0         2.0         1.0         0.0
\endverbatim

\author
Бахурин Сергей
www.dsplib.org
***************************************************************************** */
#endif
int DSPL_API flipip(double* x, int n)
{
    int k;
    double tmp;
    if(!x)
        return ERROR_PTR;
    if(n<1)
        return ERROR_SIZE;

    for(k = 0; k < n/2; k++)
    {
        tmp = x[k];
        x[k] = x[n-1-k];
        x[n-1-k] = tmp;
    }
    return RES_OK;

}



#ifdef DOXYGEN_ENGLISH
/*! ****************************************************************************
\ingroup ARRAY_GROUP
\fn int flipip_cmplx(complex_t* x, int n)
\brief Flip complex vector `x` in place

Function flips complex vector `x` length `n` in the memory
 \n
For example complex vector `x`    length 6: \n
\verbatim
x = [0+0j, 1+1j, 2+2j, 3+3j, 4+4j, 5+5j]
\endverbatim
After flipping it will be as follow:
\verbatim
x = [5+5j, 4+4j, 3+3j, 2+2j, 1+1j, 0+0j]
\endverbatim

\param[in, out] x
Pointer to the complex vector `x`. \n
Vector size is `[n x 1]`. \n
Flipped vector will be on the same address. \n

\param[in] n
Length of the vector `x`. \n \n

\return
`RES_OK` if function returns successfully. \n
Else \ref ERROR_CODE_GROUP "error code".

Example:
\code{.cpp}
complex_t y[5] = {{0.0, 0.0}, {1.0, 1.0}, {2.0, 2.0}, {3.0, 3.0}, {4.0, 4.0}};
for(i = 0; i < 5; i++)
    printf("%6.1f%+.1fj    ", RE(y[i]), IM(y[i]));
flipip_cmplx(y, 5);
printf("\n");
for(i = 0; i < 5; i++)
    printf("%6.1f%+.1fj    ", RE(y[i]), IM(y[i]));
\endcode
 \n
Program result:
\verbatim
0.0+0.0j         1.0+1.0j         2.0+2.0j         3.0+3.0j         4.0+4.0j
4.0+4.0j         3.0+3.0j         2.0+2.0j         1.0+1.0j         0.0+0.0j
\endverbatim

\author Sergey Bakhurin www.dsplib.org
***************************************************************************** */
#endif
#ifdef DOXYGEN_RUSSIAN
/*! ****************************************************************************
\ingroup ARRAY_GROUP
\fn int flipip_cmplx(complex_t* x, int n)
\brief Функция отражения комплексного вектора `x`

Функция производит отражение комплексного вектора длины `n`
в памяти данных. \n
Например исходный вектор `x`    длины 6: \n
\verbatim
x = [0+0j, 1+1j, 2+2j, 3+3j, 4+4j, 5+5j]
\endverbatim
После отражения вектор `x` будет иметь вид:
\verbatim
x = [5+5j, 4+4j, 3+3j, 2+2j, 1+1j, 0+0j]
\endverbatim

\param[in, out] x
Указатель на комплексный вектор `x`. \n
Размер вектора `[n x 1]`. \n
Результат отражения будет помещен по этому же адресу. \n
\n

\param[in] n
Размер вектора `x`. \n
\n

\return
`RES_OK` если функция выполнена успешно. \n
 В противном случае \ref ERROR_CODE_GROUP "код ошибки".

Пример:
\code{.cpp}
complex_t y[5] = {{0.0, 0.0}, {1.0, 1.0}, {2.0, 2.0}, {3.0, 3.0}, {4.0, 4.0}};
for(i = 0; i < 5; i++)
    printf("%6.1f%+.1fj    ", RE(y[i]), IM(y[i]));
flipip_cmplx(y, 5);
printf("\n");
for(i = 0; i < 5; i++)
    printf("%6.1f%+.1fj    ", RE(y[i]), IM(y[i]));
\endcode
 \n
Результат выполнения:
\verbatim
0.0+0.0j         1.0+1.0j         2.0+2.0j         3.0+3.0j         4.0+4.0j
4.0+4.0j         3.0+3.0j         2.0+2.0j         1.0+1.0j         0.0+0.0j
\endverbatim

\author
Бахурин Сергей
www.dsplib.org
***************************************************************************** */
#endif
int DSPL_API flipip_cmplx(complex_t* x, int n)
{
    int k;
    complex_t tmp;
    if(!x)
        return ERROR_PTR;
    if(n<1)
        return ERROR_SIZE;

    for(k = 0; k < n/2; k++)
    {
        RE(tmp) = RE(x[k]);
        RE(x[k]) = RE(x[n-1-k]);
        RE(x[n-1-k]) = RE(tmp);

        IM(tmp) = IM(x[k]);
        IM(x[k]) = IM(x[n-1-k]);
        IM(x[n-1-k]) = IM(tmp);
    }
    return RES_OK;
}






#ifdef DOXYGEN_ENGLISH
/*! ****************************************************************************
\ingroup ARRAY_GROUP
\fn int linspace(double x0, double x1, int n, int type, double* x)
\brief Function fills a vector with `n` linearly spaced elements
between `x0` and `x1`.

Function supports two kinds of filling according to `type` parameter: \n

Symmetric fill (parameter `type=DSPL_SYMMETRIC`): \n

\f$x(k) = x_0 + k \cdot dx\f$,
\f$dx = \frac{x_1 - x_0}{n-1}\f$, \f$k = 0 \ldots n-1.\f$

Periodic fill (parameter `type=DSPL_PERIODIC`): \n

\f$x(k) = x_0 + k \cdot dx\f$,
\f$dx = \frac{x_1 - x_0}{n}\f$, \f$k = 0 \ldots n-1.\f$

\param[in] x0
Start point \f$x_0\f$. \n \n

\param[in] x1
End point \f$x_1\f$. \n \n

\param[in] n
Number of points `x` (size of vector `x`). \n \n

\param[in] type
Fill type: \n
`DSPL_SYMMETRIC` --- symmetric, \n
`DSPL_PERIODIC` --- periodic. \n \n

\param[in,out] x
Pointer to the output linearly spaced vector `x`. \n
Vector size is `[n x 1]`. \n
Memory must be allocated. \n \n

\return
`RES_OK` if function returns successfully. \n
 Else \ref ERROR_CODE_GROUP "error code".

\note
Difference between symmetric and periodic filling we can
understand from the follow examples.    \n
Example 1. Periodic fill.
    double x[5];
    linspace(0, 5, 5, DSPL_PERIODIC, x);
\endcode
Values in the vector `x` are:
\verbatim
0,    1,    2,    3,    4
\endverbatim
 \n \n
Example 2. Symmetric fill.
\code{.cpp}
    double x[5];
    linspace(0, 5, 5, DSPL_SYMMETRIC, x);
\endcode
Values in the vector `x` are:
\verbatim
0,    1.25,    2.5,    3.75,    5
\endverbatim

\author Sergey Bakhurin www.dsplib.org
***************************************************************************** */
#endif
#ifdef DOXYGEN_RUSSIAN
/*! ****************************************************************************
\ingroup ARRAY_GROUP
\fn int linspace(double x0, double x1, int n, int type, double* x)
\brief Функция заполняет массив линейно-нарастающими, равноотстоящими
значениями от `x0` до `x1`

Заполняет массив `x` длиной `n` значениями в диапазоне
от \f$x_0\f$ до \f$x_1\f$. Функция поддерживает два типа заполнения
в соответствии с параметром `type`: \n

Симметричное заполнение согласно выражению (параметр `type=DSPL_SYMMETRIC`): \n

\f$x(k) = x_0 + k \cdot dx\f$,
\f$dx = \frac{x_1 - x_0}{n-1}\f$, \f$k = 0 \ldots n-1.\f$

Периодическое заполнение (параметр `type=DSPL_PERIODIC`) согласно выражению: \n

\f$x(k) = x_0 + k \cdot dx\f$,
\f$dx = \frac{x_1 - x_0}{n}\f$, \f$k = 0 \ldots n-1.\f$

\param[in] x0
Начальное показателя \f$x_0\f$. \n \n

\param[in] x1
Конечное значение \f$x_1\f$. \n \n

\param[in] n
Количество точек массива `x`. \n \n

\param[in] type
Тип заполнения: \n

`DSPL_SYMMETRIC` --- симметричное заполнение, \n
`DSPL_PERIODIC` --- периодическое заполнение. \n \n

\param[in,out] x
Указатель на вектор равноотстоящих значений . \n
Размер вектора `[n x 1]`. \n
Память должна быть выделена. \n \n

\return
`RES_OK` --- функция выполнена успешно.    \n
В противном случае \ref ERROR_CODE_GROUP "код ошибки". \n \n

\note
Отличие периодического и симметричного заполнения можно
понять из следующих примеров.    \n
Пример 1. Периодическое заполнение.
\code{.cpp}
    double x[5];
    linspace(0, 5, 5, DSPL_PERIODIC, x);
\endcode
В массиве `x` будут лежать значения:
\verbatim
0,    1,    2,    3,    4
\endverbatim
 \n \n
Пример 2. Симметричное заполнение.
\code{.cpp}
        double x[5];
        linspace(0, 5, 5, DSPL_SYMMETRIC, x);
\endcode
В массиве `x` будут лежать значения:
\verbatim
0,    1.25,    2.5,    3.75,    5
\endverbatim

\author
Бахурин Сергей
www.dsplib.org
***************************************************************************** */
#endif
int DSPL_API linspace(double x0, double x1, int n, int type, double* x)
{
    double dx;
    int k;

    if(n < 2)
        return ERROR_SIZE;
    if(!x)
    return ERROR_PTR;

    switch (type)
    {
        case DSPL_SYMMETRIC:
            dx = (x1 - x0)/(double)(n-1);
            x[0] = x0;
            for(k = 1; k < n; k++)
                x[k] = x[k-1] + dx;
            break;
        case DSPL_PERIODIC:
            dx = (x1 - x0)/(double)n;
            x[0] = x0;
            for(k = 1; k < n; k++)
                x[k] = x[k-1] + dx;
            break;
        default:
            return ERROR_SYM_TYPE;
    }
    return RES_OK;
}





#ifdef DOXYGEN_ENGLISH
/*! ****************************************************************************
\ingroup ARRAY_GROUP
\fn int logspace(double x0, double x1, int n, int type, double* x)
\brief Function fills a vector with `n` logarithmically spaced elements
between \f$10^{x_0}\f$ and \f$10^{x_1}\f$.


Function supports two kinds of filling according to `type` parameter: \n

Symmetric fill (parameter `type=DSPL_SYMMETRIC`): \n

\f$x(k) = 10^{x_0} \cdot dx^k\f$, here \f$dx = \sqrt[n-1]{10^{x_1 - x_0}}\f$,
\f$k = 0 \ldots n-1.\f$

Periodic fill (parameter `type=DSPL_PERIODIC`): \n

\f$x(k) = 10^{x_0} \cdot dx^k\f$, here \f$dx = \sqrt[n]{10^{x_1 - x_0}}\f$,
\f$k = 0 \ldots n-1.\f$ \n


\param[in] x0
Start exponent value \f$x_0\f$. \n \n

\param[in] x1
End exponent value \f$x_1\f$. \n \n

\param[in] n
Number of points `x` (size of vector `x`). \n \n

\param[in] type
Fill type: \n
`DSPL_SYMMETRIC` --- symmetric, \n
`DSPL_PERIODIC` --- periodic. \n \n

\param[in,out] x
Pointer to the output logarithmically spaced vector `x` . \n
Vector size is `[n x 1]`. \n
Memory must be allocated. \n \n

\return
`RES_OK` if function returns successfully. \n
 Else \ref ERROR_CODE_GROUP "error code".

\note
Difference between symmetric and periodic filling we can
understand from the follow examples.    \n
Example 1. Periodic fill.
\code{.cpp}
    double x[5];
    logspace(-2, 3, 5, DSPL_PERIODIC, x);
\endcode

Values in the vector `x` are:

\verbatim
0.01,    0.1,    1,    10,    100
\endverbatim

\n \n

Example 2. Symmetric fill.
\code{.cpp}
        double x[5];
        logspace(-2, 3, 5, DSPL_SYMMETRIC, x);
\endcode

Values in the vector `x` are:

\verbatim
0.01    0.178    3.162    56.234    1000
\endverbatim

\author Sergey Bakhurin www.dsplib.org
***************************************************************************** */
#endif
#ifdef DOXYGEN_RUSSIAN
/*! ****************************************************************************
\ingroup ARRAY_GROUP
\fn int logspace(double x0, double x1, int n, int type, double* x)
\brief    Функция заполняет массив значениями логарифмической шкале

Заполняет массив `x` длиной `n` значениями в диапазоне
от \f$10^{x_0}\f$ до \f$10^{x_1}\f$. \n
Функция поддерживает два типа заполнения в соответствии с параметром `type`: \n

Симметричное заполнение согласно выражению: \n

\f$x(k) = 10^{x_0} \cdot dx^k\f$, где \f$dx = \sqrt[n-1]{10^{x_1 - x_0}}\f$,
\f$k = 0 \ldots n-1.\f$

Периодическое заполнение согласно выражению:

\f$x(k) = 10^{x_0} \cdot dx^k\f$, где \f$dx = \sqrt[n]{10^{x_1 - x_0}}\f$,
\f$k = 0 \ldots n-1.\f$ \n

\param[in] x0
Начальное значение показателя \f$x_0\f$. \n \n

\param[in] x1
Конечное значение показателя \f$x_1\f$. \n \n

\param[in] n
Количество точек массива `x`. \n \n

\param[in] type
Тип заполнения: \n
`DSPL_SYMMETRIC` --- симметричное заполнение, \n
`DSPL_PERIODIC` --- периодическое заполнение. \n \n

\param[in,out] x
Указатель на вектор значений в логарифмической шкале. \n
Размер вектора `[n x 1]`. \n
Память должна быть выделена. \n \n

\return
`RES_OK` --- функция выполнена успешно.    \n
В противном случае \ref ERROR_CODE_GROUP "код ошибки".

\note
Отличие периодического и симметричного заполнения можно
понять из следующих примеров.    \n

Пример 1. Периодическое заполнение.
\code{.cpp}
    double x[5];
    logspace(-2, 3, 5, DSPL_PERIODIC, x);
\endcode
В массиве `x` будут лежать значения:
\verbatim
0.01,    0.1,    1,    10,    100
\endverbatim
\n \n

Пример 2. Симметричное заполнение.
\code{.cpp}
    double x[5];
    logspace(-2, 3, 5, DSPL_SYMMETRIC, x);
\endcode

В массиве `x` будут лежать значения:

\verbatim
0.01    0.178    3.162    56.234    1000
\endverbatim

\author Бахурин Сергей www.dsplib.org
***************************************************************************** */
#endif
int DSPL_API logspace(double x0, double x1, int n, int type, double* x)
{
    double mx, a, b;
    int k;

    if(n < 2)
        return ERROR_SIZE;
    if(!x)
        return ERROR_PTR;

    a = pow(10.0, x0);
    b = pow(10.0, x1);

    switch (type)
    {
        case DSPL_SYMMETRIC:
            mx = pow(b/a, 1.0/(double)(n-1));
            x[0] = a;
            for(k = 1; k < n; k++)
                x[k] = x[k-1] * mx;
            break;
        case DSPL_PERIODIC:
            mx = pow(b/a, 1.0/(double)n);
            x[0] = a;
            for(k = 1; k < n; k++)
                x[k] = x[k-1] * mx;
            break;
        default:
            return ERROR_SYM_TYPE;
    }
    return RES_OK;
}




#ifdef DOXYGEN_ENGLISH
/*! ****************************************************************************
\ingroup ARRAY_GROUP
\fn int ones(double* x, int n)
\brief Function fills all real vector `x` by ones values.

\param[in, out] x
Pointer to the vector `x`. \n
Vector size is `[n x 1]`. \n
All elements on this vector will be set to one. \n
\n

\param[in] n
Size of vector `x`. \n
\n

\return
`RES_OK` if function returns successfully. \n
 Else \ref ERROR_CODE_GROUP "error code".

Example:
\code{.cpp}
double y[5] = {0};
int i;
ones(y, 5);
for(i = 0; i < 5; i++)
    printf("%6.1f%    ", y[i]);
\endcode
 \n
Vector `y` values are:
\verbatim
    1.0    1.0    1.0    1.0    1.0
\endverbatim

\author Sergey Bakhurin www.dsplib.org
***************************************************************************** */
#endif
#ifdef DOXYGEN_RUSSIAN
/*! ****************************************************************************
\ingroup ARRAY_GROUP
\fn int ones(double* x, int n)
\brief Функция заполнения вещественного массива единицами

\param[in, out] x
Указатель на вещественный вектор `x`. \n
Размер вектора `[n x 1]`. \n
Значения данного вектора будут установлены в единицу. \n
\n

\param[in] n
Размер вектора `x`. \n
\n

\return
`RES_OK` если функция выполнена успешно. \n
В противном случае \ref ERROR_CODE_GROUP "код ошибки".

Пример:
\code{.cpp}
double y[5] = {0};
int i;
ones(y, 5);
for(i = 0; i < 5; i++)
    printf("%6.1f%    ", y[i]);
\endcode
 \n
Результат выполнения:
\verbatim
    1.0    1.0    1.0    1.0    1.0
\endverbatim

\author Бахурин Сергей www.dsplib.org
***************************************************************************** */
#endif
int DSPL_API ones(double* x, int n)
{
    int i;
    if(!x)
        return ERROR_PTR;
    if(n<1)
        return ERROR_SIZE;
    for(i = 0; i < n; i++)
        x[i] = 1.0;
 return RES_OK;
}


#ifdef DOXYGEN_ENGLISH
/*! ****************************************************************************
\ingroup ARRAY_GROUP
\fn int verif(double* x,    double* y, size_t n, double eps, double* err)
\brief Real arrays verification

Function calculates a maximum relative error between two real arrays `x`
and `y` (both length equals `n`):

\f[
e = \max \left( \frac{|x(k) - y(k)| }{ |x(k)|} \right),
\quad if \quad |x(k)| > 0,
\f]
or
\f[
e = \max(|x(k) - y(k)| ), ~\qquad if \quad~|x(k)| = 0,
\f]

This function can be used for algorithms verification if vector `x` is user
algorithm result and vector `y` -- reference vector.

\param[in] x
Pointer to the first vector `x`. \n
Vector size is `[n x 1]`. \n \n

\param[in] y
Pointer to the second vector `y`. \n
Vector size is `[n x 1]`. \n \n

\param[in] n
Size of vectors `x` and `y`. \n \n

\param[in] eps
Relative error threshold. \n
If error less than `eps`, then function returns
`DSPL_VERIF_SUCCESS`, else `DSPL_VERIF_FAILED`.    \n \n

\param[in, out] err
Pointer to the variable which keep
maximum relative error. \n
Pointer can be `NULL`, maximum error will not be returned
in this case. \n \n

\return
`DSPL_VERIF_SUCCESS` if maximum relative error less than `eps`. \n
Otherwise `DSPL_VERIF_FAILED`.

\author Sergey Bakhurin www.dsplib.org
***************************************************************************** */
#endif
#ifdef DOXYGEN_RUSSIAN
/*! ****************************************************************************
\ingroup ARRAY_GROUP
\fn int verif(double* x, double* y, size_t n, double eps, double* err)
\brief Верификация вещественных массивов

Функция производит расчет максимальной относительной ошибки между вещественными
векторами `x` и `y` равной длины `n`:

\f[
e = \max \left( \frac{|x(k) - y(k)| }{ |x(k)|} \right), \quad    \quad |x(k)| > 0,
\f]
или
\f[
e = \max(|x(k) - y(k)| ), \qquad \quad~|x(k)| = 0,
\f]
и возвращает `DSPL_VERIF_SUCCESS` если
разница \f$ e\f$ меньше `eps`.
В противном случае возвращает `DSPL_VERIF_FAILED`. \n
Данная функция используется для верификации работы алгоритмов если вектор `x`
результат работы алгоритма пользователя, а `y` -- результат работы этого же
алгоритма сторонней функцией.

\param[in] x
Указатель на первый вектор `x`. \n
Размер вектора `[n x 1]`. \n \n

\param[in] y
Указатель на второй вектор `y`. \n
Размер вектора `[n x 1]`. \n \n

\param[in] n
Размер векторов `x` и `y`. \n \n

\param[in] eps
Допустимая относительная ошибка. \n
Если максимальная относительная ошибка меньше `eps`, то функция возвращает
`DSPL_VERIF_SUCCESS`, в противном случае `DSPL_VERIF_FAILED`.    \n \n

\param[in, out] err
Указатель на переменную максимальной относительной ошибки. \n
По данному адресу будет записано значение максимальной относительной ошибки. \n
Указатель может быть `NULL`, значение ошибки в этом случае
не возвращается. \n \n

\return
`DSPL_VERIF_SUCCESS` если относительная ошибка меньше `eps`. \n
 В противном случае `DSPL_VERIF_FAILED`.

\author Бахурин Сергей www.dsplib.org
***************************************************************************** */
#endif
int DSPL_API verif(double* x,    double* y, size_t n, double eps, double* err)
{
    double d, maxd;
    size_t k;
    int res;
    if(!x || !y)
        return ERROR_PTR;
    if(n < 1)
        return ERROR_SIZE;
    if(eps <= 0.0 )
        return ERROR_NEGATIVE;

    maxd = -100.0;

    for(k = 0; k < n; k++)
    {
        d = fabs(x[k] - y[k]);
        if(fabs(x[k]) > 0.0)
        {
            d = d / fabs(x[k]);
            if(d > maxd)
                maxd = d;
        }
    }
    if(err)
        *err = maxd;

    if(maxd > eps)
        res = DSPL_VERIF_FAILED;
    else
        res = DSPL_VERIF_SUCCESS;

    return res;
}



#ifdef DOXYGEN_ENGLISH
/*! ****************************************************************************
\ingroup ARRAY_GROUP
\fn int verif_cmplx(complex_t* x, complex_t* y, size_t n,
                    double eps, double* err)
\brief
Complex arrays verification

Function calculates a maximum relative error between two complex arrays `x`
and `y` (both length equals `n`):

\f[
e = \max \left( \frac{|x(k) - y(k)| }{ |x(k)|} \right),
\quad if \quad |x(k)| > 0,
\f]
or
\f[
e = \max(|x(k) - y(k)| ), ~\qquad if \quad~|x(k)| = 0,
\f]
Return `DSPL_VERIF_SUCCESS` if maximum relative error \f$ e\f$ less than `eps`.
Else returns `DSPL_VERIF_FAILED`. \n

This function can be used for algorithms verification if vector `x` is user
algorithm result and vector `y` -- reference vector.

\param[in] x
Pointer to the first vector `x`. \n
Vector size is `[n x 1]`. \n \n

\param[in] y
Pointer to the second vector `y`. \n
Vector size is `[n x 1]`. \n \n

\param[in] n
Size of vectors `x` and `y`. \n \n

\param[in] eps
Relative error threshold. \n
If error less than `eps`, then function returns
`DSPL_VERIF_SUCCESS`, else `DSPL_VERIF_FAILED`.    \n \n

\param[in, out] err
Pointer to the variable which keep
maximum relative error. \n
Pointer can be `NULL`, maximum error will not be returned
in this case. \n \n

\return
`DSPL_VERIF_SUCCESS` if maximum relative error less than `eps`. \n
Otherwise `DSPL_VERIF_FAILED`.

\author
Sergey Bakhurin
www.dsplib.org
***************************************************************************** */
#endif
#ifdef DOXYGEN_RUSSIAN
/*! ****************************************************************************
\ingroup ARRAY_GROUP
\fn int verif_cmplx(complex_t* x,    complex_t* y, size_t n,
                    double eps, double* err)
\brief Верификация комплексных массивов

Функция производит расчет максимальной относительной ошибки между комплексными
векторами `x` и `y` равной длины `n`:

\f[
e = \max \left( \frac{|x(k) - y(k)|}{|x(k)|} \right), \quad    \quad |x(k)| > 0,
\f]
или
\f[
e = \max(|x(k) - y(k)| ), ~\qquad    \quad~|x(k)| = 0,
\f]
и возвращает `DSPL_VERIF_SUCCESS` если
разница \f$ e\f$ меньше `eps`.
В противном случае возвращает `DSPL_VERIF_FAILED`. \n
Данная функция используется для верификации работы алгоритмов если вектор `x`
результат работы алгоритма пользователя, а `y` -- результат работы этого же
алгоритма сторонней функцией.

\param[in] x
Указатель на первый вектор `x`. \n
Размер вектора `[n x 1]`. \n \n

\param[in] y
Указатель на второй вектор `y`. \n
Размер вектора `[n x 1]`. \n \n

\param[in] n
Размер векторов `x` и `y`. \n \n

\param[in] eps
Допустимая относительная ошибка. \n
Если максимальная относительная ошибка меньше `eps`, то    функция возвращает
`DSPL_VERIF_SUCCESS`, в противном случае `DSPL_VERIF_FAILED`.    \n \n

\param[in, out] err
Указатель на переменную максимальной относительной ошибки. \n
По данному адресу будет записано значение максимальной относительной ошибки. \n
Указатель может быть `NULL`, значение ошибки в этом
случае не возвращается. \n \n

\return
`DSPL_VERIF_SUCCESS` если функция выполнена успешно. \n
 В противном случае `DSPL_VERIF_FAILED`.

\author
Бахурин Сергей www.dsplib.org
***************************************************************************** */
#endif
int DSPL_API verif_cmplx(complex_t* x,    complex_t* y, size_t n,
                 double eps, double* err)
{

    complex_t d;
    double mx, md, maxd;
    size_t k;
    int res;
    if(!x || !y)
        return ERROR_PTR;
    if(n < 1)
        return ERROR_SIZE;
    if(eps <= 0.0 )
        return ERROR_NEGATIVE;

    maxd = -100.0;

    for(k = 0; k < n; k++)
    {
        RE(d) = RE(x[k]) - RE(y[k]);
        IM(d) = IM(x[k]) - IM(y[k]);
        md = ABS(d);
        mx = ABS(x[k]);
        if(mx > 0.0)
        {
            md = md / mx;
            if(md > maxd)
                maxd = md;
        }
    }
    if(err)
        *err = maxd;

    if(maxd > eps)
        res = DSPL_VERIF_FAILED;
    else
        res = DSPL_VERIF_SUCCESS;

    return res;
}
