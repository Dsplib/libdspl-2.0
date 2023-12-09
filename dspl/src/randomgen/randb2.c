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


