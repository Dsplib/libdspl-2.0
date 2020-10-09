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
#include <math.h>
#include <time.h>

#include "dspl.h"
#include "dspl_internal.h"
#include "mt19937.h"



#ifdef DOXYGEN_ENGLISH
/*! ****************************************************************************
\ingroup SPEC_MATH_RAND_GEN_GROUP
\fn int random_init(random_t* prnd, int type, void* seed)
\brief Pseudorandom numbers generators initialization. 

\param[in,out] prnd
Pointer to the pseudorandom generators parameters and state vectors.\n\n

\param[in]  type
Pseudorandom generator algorithm:
\verbatim
RAND_TYPE_MRG32K3A - 32-bits MRG32K3A generator
RAND_TYPE_MT19937  - 64-bits MT19937-64 generator
\endverbatim

\param[in]  seed
Pointer to the generator start initialization. \n
Type of this pointer is `void*` because different generators are using
different initial values. For example, if we initialize the MRG32K3A generator,
then the `type` parameter is specified as` RAND_TYPE_MRG32K3A`, and `seed`
converts to `double` pointer:
\code
    random_t rnd = {0};
    double seed = 1234.0;
    random_init(&rnd, RAND_TYPE_MRG32K3A, (void*)&seed);
\endcode
For 64-bits Mersenne Twister pseudorandom number generator 
(`type` sets as `RAND_TYPE_MT19937`), `seed` converts to the 
`unsigned long long` pointer:
\code
    random_t rnd = {0};
    unsigned long long seed = 1234353456;
    random_init(&rnd, RAND_TYPE_MT19937, (void*)&seed);
\endcode
Pseudorandom numbers will be repeated each program restart 
if `seed` value is the same.\n
The `seed` pointer can be `NULL`. Pseudorandom generators will be initialized
by pseudorandom numbers in this case and program will generate different 
pseudorandom numbers each restart.

\author Sergey Bakhurin. www.dsplib.org
****************************************************************************  */
#endif
#ifdef DOXYGEN_RUSSIAN
/*! ****************************************************************************
\ingroup SPEC_MATH_RAND_GEN_GROUP
\fn int random_init(random_t* prnd, int type, void* seed)
\brief Инициализация датчиков псевдослучайных чисел. 

\param[in,out] prnd
Указатель на структуру параметров и векторов состояния 
датчиков псевдослучайных чисел, которая будет инициализирована. \n\n

\param[in]  type
Тип датчика псевдослучайных чисел:
\verbatim
RAND_TYPE_MRG32K3A - 32-битный датчик MRG32K3A
RAND_TYPE_MT19937  - 64-битный датчик MT19937-64
\endverbatim

\param[in]  seed
Указатель на начальную инициализацию датчика. \n
Данный указатель имеет тип `void*`, поскольку параметр инициализации
зависит от типа датчика. Например если инициализируем датчик MRG32K3A,
т.е. параметр `type` задан как `RAND_TYPE_MRG32K3A`, то данный указатель 
приводится к типу `double`:
\code
    random_t rnd = {0};
    double seed = 1234.0;
    random_init(&rnd, RAND_TYPE_MRG32K3A, (void*)&seed);
\endcode
Если же используется 64-битный датчик Вихрь Мерсенна 
(`type` задан как `RAND_TYPE_MT19937`), то `seed` приводится к типу 
`unsigned long long`:
\code
    random_t rnd = {0};
    unsigned long long seed = 1234353456;
    random_init(&rnd, RAND_TYPE_MT19937, (void*)&seed);
\endcode
При фиксированном начальном значении датчика, псевдослучайные числа будут
повторяться при каждом запуске программы. \n
Указатель `seed` может быть `NULL`. В этом случае начальная инициализация 
датчиков будет задаваться случайными значениями и генерируемые псевдослучайные 
числа будут различными при каждом запуске программы.

\author Бахурин Сергей. www.dsplib.org
****************************************************************************  */
#endif
int DSPL_API random_init(random_t* prnd, int type, void* seed)
{
    srand(time(NULL));

    if(!prnd)
        return RES_OK;

    switch(type)
    {
        case RAND_TYPE_MRG32K3A:
            /* MRG32k3a init */
            prnd->mrg32k3a_x[0] = prnd->mrg32k3a_x[1] = 1.0;
            prnd->mrg32k3a_y[0] = prnd->mrg32k3a_y[1] = 
                                  prnd->mrg32k3a_y[2] = 1.0;
            if(seed)
                prnd->mrg32k3a_x[2] = *((double*)seed);
            else
                prnd->mrg32k3a_x[2] = (double) rand() * rand();
            break;
        case RAND_TYPE_MT19937:
            if(seed)
                mt19937_init_genrand64(*((unsigned long long*)seed), prnd);
            else
                mt19937_init_genrand64((unsigned long long)rand()*rand(), prnd);
            break;
        default:
            return ERROR_RAND_TYPE;
    }
    prnd->type = type;

    return RES_OK;
}




#ifdef DOXYGEN_ENGLISH
/*! ****************************************************************************
\ingroup SPEC_MATH_RAND_GEN_GROUP
\fn int randb(double* x, int n, random_t* prnd)
\brief Binary unipolar [0, 1] pseudorandom vector. 

The function generates a unipolar pseudo-random vector,
each element of which takes an equally probable value of 0 or 1

\param[in,out] x  
Pointer to the unipolar pseudo-random vector.  \n
Vector size is `[n x 1]`. \n
Memory must be allocated. \n\n

\param[in] n
Size of vector `x`. \n\n

\param[in] prnd
Pointer to the `random_t` structure. \n
The structure must be pre-filled with the \ref random_init function. \n
This pointer can be `NULL`, then it will be used
built-in pseudorandom generator defined by the C language standard.
However, this mode is not recommended, 
for example in cryptography and other tasks.
There is no guarantee of the quality of the pseudorandom numbers produced if
the `prnd` parameter is set to` NULL`. \n\n

\return
`RES_OK` --- if pseudorandom vector is calculated successfully. \n
Else \ref ERROR_CODE_GROUP "code error".

Example:

\include randb_test.c

Program genrates unipolar [0, 1] and bipolar[-1, 1] pseudorandom binary vectors.

As a result of the program run, you can see the graph:
\image html randb_test.png

\author Sergey Bakhurin. www.dsplib.org
***************************************************************************** */
#endif
#ifdef DOXYGEN_RUSSIAN
/*! ****************************************************************************
\ingroup SPEC_MATH_RAND_GEN_GROUP
\fn int randb(double* x, int n, random_t* prnd)
\brief Генерация бинарного униполярного [0, 1] псевдослучайного вектора 

Функция генерирует униполярный псевдослучайный вектор, 
каждый элемент которого принимает равновероятное значение 0 или 1. 

\param[in,out] x  
Указатель на вектор случайных бинарных чисел.  \n
Размер вектора `[n x 1]`. \n
Память должна быть выделена. \n\n

\param[in] n
Размер вектора `x`. \n\n

\param[in] prnd
Указатель на структуру `random_t` параметров датчиков
псевдослучайных чисел. \n
Структура должна быть предварительно заполнена функцией \ref random_init. \n
Данный указатель может быть `NULL`, тогда будет использоваться 
встроенный датчик, определенный стандартом языка Си. Однако для серьезных нужд, 
например в криптографии, данный режим использовать не рекомендуется. 
Нет гарантии в качестве произведенной случайной последовательности если 
параметр `prnd` задан как `NULL`. \n\n

\return
`RES_OK` --- вектор целых псевдослучайных чисел рассчитан успешно. \n
В противном случае \ref ERROR_CODE_GROUP "код ошибки".

Пример использования функции:

\include randb_test.c

Программа рассчитывает униполярный [0, 1] и биполярный [-1, 1] бинарные
псевдослучайные векторы.

В результате выполнения программы можно увидеть график:
\image html randb_test.png

\author Бахурин Сергей. www.dsplib.org
***************************************************************************** */
#endif
int DSPL_API randb(double* x, int n, random_t* prnd)
{
    double z[RAND_BUFSIZE];
    int i, cnt, err;
    if(!x)
        return ERROR_PTR;
    if(n < 1)
        return ERROR_SIZE;
    cnt = 0;
    while(cnt < n)
    {
        i = cnt % RAND_BUFSIZE;
        if(!i)
        {
            err = randu(z, RAND_BUFSIZE, prnd);
            if(err != RES_OK)
                return err;
        }
        x[cnt] = z[i] > 0.5 ? 1.0 : 0.0;
        cnt++;
    }
    return RES_OK;
}




#ifdef DOXYGEN_ENGLISH
/*! ****************************************************************************
\ingroup SPEC_MATH_RAND_GEN_GROUP
\fn int randb2(double* x, int n, random_t* prnd)
\brief Binary bipolar [-1, 1] pseudorandom vector. 

The function generates a unipolar pseudo-random vector,
each element of which takes an equally probable value of -1 or 1

\param[in,out] x  
Pointer to the bipolar pseudorandom vector.  \n
Vector size is `[n x 1]`. \n
Memory must be allocated. \n\n

\param[in] n
Size of vector `x`. \n\n

\param[in] prnd
Pointer to the `random_t` structure. \n
The structure must be pre-filled with the \ref random_init function. \n
This pointer can be `NULL`, then it will be used
built-in pseudorandom generator defined by the C language standard.
However, this mode is not recommended, 
for example in cryptography and other tasks.
There is no guarantee of the quality of the pseudorandom numbers produced if
the `prnd` parameter is set to` NULL`. \n\n

\return
`RES_OK` --- if pseudorandom vector is calculated successfully. \n
Else \ref ERROR_CODE_GROUP "code error".

Example:

\include randb_test.c

Program genrates unipolar [0, 1] and bipolar[-1, 1] pseudorandom binary vectors.

As a result of the program run, you can see the graph:
\image html randb_test.png

\author Sergey Bakhurin. www.dsplib.org
***************************************************************************** */
#endif
#ifdef DOXYGEN_RUSSIAN
/*! ****************************************************************************
\ingroup SPEC_MATH_RAND_GEN_GROUP
\fn int randb2(double* x, int n, random_t* prnd)
\brief Генерация бинарного биполярного [-1, 1] псевдослучайного вектора 

Функция генерирует биполярный псевдослучайный вектор, 
каждый элемент которого принимает равновероятное значение -1 или 1. 

\param[in,out] x  
Указатель на вектор случайных бинарных чисел.  \n
Размер вектора `[n x 1]`. \n
Память должна быть выделена. \n\n

\param[in] n
Размер вектора `x`. \n\n

\param[in] prnd
Указатель на структуру `random_t` параметров датчиков
псевдослучайных чисел. \n
Структура должна быть предварительно заполнена функцией \ref random_init. \n
Данный указатель может быть `NULL`, тогда будет использоваться 
встроенный датчик, определенный стандартом языка Си. Однако для серьезных нужд, 
например в криптографии, данный режим использовать не рекомендуется. 
Нет гарантии в качестве произведенной случайной последовательности если 
параметр `prnd` задан как `NULL`. \n\n

\return
`RES_OK` --- вектор целых псевдослучайных чисел рассчитан успешно. \n
В противном случае \ref ERROR_CODE_GROUP "код ошибки".

Пример использования функции:

\include randb_test.c

Программа рассчитывает униполярный [0, 1] и биполярный [-1, 1] бинарные
псевдослучайные векторы.

В результате выполнения программы можно увидеть график:
\image html randb_test.png

\author Бахурин Сергей. www.dsplib.org
***************************************************************************** */
#endif
int DSPL_API randb2(double* x, int n, random_t* prnd)
{
    double z[RAND_BUFSIZE];
    int i, cnt, err;
    if(!x)
        return ERROR_PTR;
    if(n < 1)
        return ERROR_SIZE;
    cnt = 0;
    while(cnt < n)
    {
        i = cnt % RAND_BUFSIZE;
        if(!i)
        {
            err = randu(z, RAND_BUFSIZE, prnd);
            if(err != RES_OK)
                return err;
        }
        x[cnt] = z[i] > 0.5 ? 1.0 : -1.0;
        cnt++;
    }
    return RES_OK;
}




#ifdef DOXYGEN_ENGLISH

#endif
#ifdef DOXYGEN_RUSSIAN

#endif
int randu_mrg32k3a (double* u, int n, random_t* prnd)
{

    if(!u || !prnd)
        return ERROR_PTR;
    if(n < 1)
        return ERROR_SIZE;

    long z;
    double xn, yn, *x, *y;
    int k;

    x = prnd->mrg32k3a_x;
    y = prnd->mrg32k3a_y;
    for(k = 0; k < n; k++)
    {
        /* Component x[n] */
        xn = MRG32K3A_A12 * x[1] - MRG32K3A_A13 * x[2];

        z = (long)(xn / MRG32K3A_M1);
        xn -= (double)z * MRG32K3A_M1;
        if (xn < 0.0)
            xn += MRG32K3A_M1;

        x[2] = x[1];
        x[1] = x[0];
        x[0] = xn;

        /* Component y[n] */
        yn = MRG32K3A_A21 * y[0] - MRG32K3A_A23 * y[2];
        z = (long)(yn / MRG32K3A_M2);
        yn -= (double)z * MRG32K3A_M2;
        if (yn < 0.0)
             yn += MRG32K3A_M2;

        y[2] = y[1];
        y[1] = y[0];
        y[0] = yn;

        /* Combination */
        u[k] = (xn <= yn) ? ((xn - yn + MRG32K3A_M1) * MRG32K3A_NORM):
                            (xn - yn) * MRG32K3A_NORM;
    }
    return RES_OK;
}








#ifdef DOXYGEN_ENGLISH

#endif
#ifdef DOXYGEN_RUSSIAN
/*! ****************************************************************************
\ingroup SPEC_MATH_RAND_GEN_GROUP
\fn int randi(int* x, int n, int start, int stop, random_t* prnd)
\brief  Генерация целочисленного вектора равномерно 
распределенных псевдослучайных чисел.

Функция генерирует псевдослучайный вектор целых чисел в диапазоне от `start`
до `stop` включительно. 

\param[in,out] x  
Указатель на вектор случайных чисел.  \n
Размер вектора `[n x 1]`. \n
Память должна быть выделена. \n\n

\param[in] n
Размер вектора `x`. \n\n

\param[in] start
Начало диапазона целых чисел. \n\n

\param[in] stop
Конец диапазона целых чисел. \n\n

\param[in] prnd
Указатель на структуру `random_t` параметров датчиков
псевдослучайных чисел. \n
Структура должна быть предварительно заполнена функцией \ref random_init. \n
Данный указатель может быть `NULL`, тогда будет использоваться 
встроенный датчик, определенный стандартом языка Си. Однако для серьезных нужд, 
например в криптографии, данный режим использовать не рекомендуется. 
Нет гарантии в качестве произведенной случайной последовательности если 
параметр `prnd` задан как `NULL`. \n\n

\return
`RES_OK` --- вектор целых псевдослучайных чисел рассчитан успешно. \n
В противном случае \ref ERROR_CODE_GROUP "код ошибки".

Пример использования функции:

\include randi_test.c

Программа рассчитывает целочисленный вектор 
псевдослучайных чисел в диапазоне [-4, 3].

В результате выполнения программы можно увидеть график:
\image html randi_test.png

\author Бахурин Сергей. www.dsplib.org
***************************************************************************** */
#endif
int DSPL_API randi(int* x, int n, int start, int stop, random_t* prnd)
{
    double z[RAND_BUFSIZE];
    double dx;
    int i, cnt, err;
    if(!x)
        return ERROR_PTR;
    if(n < 1)
        return ERROR_SIZE;
    
    dx = (double)stop - (double)start;
    cnt = 0;
    while(cnt < n)
    {
        i = cnt % RAND_BUFSIZE;
        if(!i)
        {
            err = randu(z, RAND_BUFSIZE, prnd);
            if(err != RES_OK)
                return err;
        }
        x[cnt] = start + (int)round(z[i] * dx);
        cnt++;
    }
    return RES_OK;
}







#ifdef DOXYGEN_ENGLISH

#endif
#ifdef DOXYGEN_RUSSIAN
/*! ****************************************************************************
\ingroup SPEC_MATH_RAND_GEN_GROUP
\fn int randn(double* x, int n, double mu, double sigma, random_t* prnd)
\brief Генерация вектора нормально распределенных псевдослучайных чисел.

Функция использует преобразование Бокса-Мюллера для приведения  
равномерно-распределенных псевдослучайных чисел к нормальному распределению 
с математическим ожиданием \f$\mu\f$ и среднеквадратическим 
отклонением \f$\sigma\f$.

\param[in,out] x  
Указатель на вектор случайных чисел.  \n
Размер вектора `[n x 1]`. \n
Память должна быть выделена. \n\n
      
\param[in] n
Размер вектора случайных чисел. \n\n

\param[in] mu
Математическое ожидание \f$\mu\f$. \n\n


\param[in] sigma
Среднеквадратическое отклонение \f$\sigma\f$. \n
Дисперсия сгенерированных чисел равна \f$\sigma^2\f$. \n\n

\param[in] prnd
Указатель на структуру `random_t` параметров датчиков
псевдослучайных чисел. \n
Структура должна быть предварительно заполнена функцией \ref random_init. \n
Данный указатель может быть `NULL`, тогда будет использоваться 
встроенный датчик, определенный стандартом языка Си. Однако для серьезных нужд, 
например в криптографии, данный режим использовать не рекомендуется. 
Нет гарантии в качестве произведенной случайной последовательности если 
параметр `prnd` задан как `NULL`. \n\n

\return
`RES_OK` --- вектор нормально распределенных 
псевдослучайных чисел рассчитан успешно. \n
В противном случае \ref ERROR_CODE_GROUP "код ошибки".

Пример использования функции:

\include randn_test.c

Программа рассчитывает независимые векторы нормально распределенных 
псевдослучайных чисел, \f$\mu = 0\f$ и \f$\sigma=1\f$.

В результате выполнения программы можно увидеть график:
\image html randn_test.png

\author Бахурин Сергей. www.dsplib.org
***************************************************************************** */
#endif
int DSPL_API randn(double* x, int n, double mu, double sigma, random_t* prnd)
{
    int k, m;
    double x1[RAND_BUFSIZE], x2[RAND_BUFSIZE];
    int res;
    if(!x)
        return ERROR_PTR;

     if(n<1)
        return ERROR_SIZE;

    if(sigma < 0.0)
        return ERROR_RAND_SIGMA;

    k=0;
    while(k < n)
    {
        if((res = randu(x1, RAND_BUFSIZE, prnd)) != RES_OK)
            goto exit_label;
        if((res = randu(x2, RAND_BUFSIZE, prnd)) != RES_OK)
            goto exit_label;
        m = 0;
        while(k < n && m < RAND_BUFSIZE)
        {
            if(x1[m] != 0.0)
            {
                x[k] = sqrt(-2.0*log(x1[m]))*cos(M_2PI*x2[m])*sigma + mu;
                k++;
                m++;
            }
        }
    }

    res = RES_OK;
exit_label:
    return res;
}




#ifdef DOXYGEN_ENGLISH

#endif
#ifdef DOXYGEN_RUSSIAN

#endif
int DSPL_API randn_cmplx(complex_t* x, int n, complex_t* mu, 
                         double sigma, random_t* prnd)
{
    int err, i;
    
    err = randn((double*)x, 2*n, 0.0, sigma / M_SQRT2, prnd);
    if(err!= RES_OK)
        return err;
    if(mu)
    {
        for(i = 0; i < n; i++)
        {
            RE(x[i]) += RE(mu[0]);
            IM(x[i]) += IM(mu[0]);
        }
    }
    return err;
}






#ifdef DOXYGEN_ENGLISH

#endif
#ifdef DOXYGEN_RUSSIAN
/*! ****************************************************************************
\ingroup SPEC_MATH_RAND_GEN_GROUP
\fn int randu(double* x, int n, random_t* prnd)
\brief Генерация вектора равномерно-распределенных в интервале 
от 0 до 1 псевдослучайных чисел.

\param[in,out] x  
Указатель на вектор случайных чисел.  \n
Размер вектора `[n x 1]`. \n
Память должна быть выделена. \n\n
      
\param[in] n
Размер вектора случайных чисел. \n\n

\param[in] prnd
Указатель на структуру `random_t` параметров датчиков
псевдослучайных чисел. \n
Структура должна быть предварительно заполнена функцией \ref random_init. \n
 Данный указатель может быть `NULL`, тогда будет использоваться 
встроенный датчик, определенный стандартом языка Си. Однако для серьезных нужд, 
например в криптографии, данный режим использовать не рекомендуется. 
Нет гарантии в качестве произведенной случайной последовательности если 
параметр `prnd` задан как `NULL`. \n\n

\return
`RES_OK` --- вектор равномерно-распределенных 
псевдослучайных чисел рассчитан успешно. \n
В противном случае \ref ERROR_CODE_GROUP "код ошибки".

Пример использования функции с различными датчиками псевдослучайных чисел
приведен в следующем листинге:

\include randu_test.c

Программа рассчитывает независимые векторы равномерно-распределенных 
от 0 до 1 псевдослучайных чисел и выводит их на график для трех различных 
датчиков: MRG32K3A, MT19937-64 и встроенный датчик, определенный 
стандартом языка Си.

В результате выполнения программы можно увидеть график:

\image html randu_test.png

Однако при детальном исследовании датчиков, можно обнаружить, что встроенный 
датчик, определенный стандартом языка Си, 
выдает значения на фиксированной сетке. 

Чтобы проверить это можно выполнить следующую программу:

\include randu_accuracy_test.c

Данная программа аккумулирует только значения датчиков в интервале 
от 0 до 0.001 и выводит их на график:

\image html randu_acc_test.png

Из графика хорошо видно, что данные встроенного датчика выдаются на 
равноотстоящей сетке значений, в отличии от датчиков MRG32K3A и MT19937-64,
которые сохранили псевдослучайный характер.

\author Бахурин Сергей. www.dsplib.org
***************************************************************************** */
#endif
int DSPL_API randu(double* x, int n, random_t* prnd)
{
    int i;

    if(!x)
        return ERROR_PTR;
    if(n < 0)
        return ERROR_SIZE;

    if(prnd)
    {
        switch(prnd->type)
        {
            case RAND_TYPE_MRG32K3A:
                return randu_mrg32k3a(x, n, prnd);
            case RAND_TYPE_MT19937:
                for(i = 0; i < n; i++)
                    x[i] = mt19937_genrand64_real1(prnd);
                return RES_OK;
            default:
                return ERROR_RAND_TYPE;
        }
    }
    else
    {
        if(!x)
            return ERROR_PTR;
        if(n<1)
            return ERROR_SIZE;
        for(i = 0; i < n; i++)
            x[i] = (double)rand()/RAND_MAX;
    }

    return RES_OK;
}

