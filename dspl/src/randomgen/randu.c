/*
* Copyright (c) 2015-2022 Sergey Bakhurin
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

