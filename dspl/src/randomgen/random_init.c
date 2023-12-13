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
/*! ****************************************************************************
\ingroup SPEC_MATH_RAND_GEN_GROUP

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

