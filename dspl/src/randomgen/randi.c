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


