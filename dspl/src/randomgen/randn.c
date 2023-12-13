/*
* Copyright (c) 2015-2024 Sergey Bakhurin
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
#include "randomgen.h"



#ifdef DOXYGEN_ENGLISH

#endif
#ifdef DOXYGEN_RUSSIAN
/*! ****************************************************************************
\ingroup SPEC_MATH_RAND_GEN_GROUP

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

