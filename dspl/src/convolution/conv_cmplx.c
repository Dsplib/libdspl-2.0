/*
* Copyright (c) 2015-2020 Sergey Bakhurin
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
* MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
* GNU General Public License for more details.
*
* You should have received a copy of the GNU General Public License
* along with Foobar.  If not, see <http://www.gnu.org/licenses/>.
*/


#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "dspl.h"




#ifdef DOXYGEN_ENGLISH
/*! ****************************************************************************
\ingroup FILTER_CONV_GROUP
\fn int conv_cmplx(complex_t* a, int na, complex_t* b, int nb, complex_t* c) 
\brief Complex vectors linear convolution.

Function convolves two complex vectors \f$ c = a * b\f$ length `na` and `nb`.
The output convolution is a vector `c` with length equal to  `na + nb - 1`. 

\param[in]  a
Pointer to the first vector `a`. \n
Vector size is `[na x 1]`. \n \n

\param[in]  na
Size of the first vector `a`. \n \n

\param[in]  b
Pointer to the second vector `b`. \n
Vector size is `[nb x 1]`. \n \n

\param[in]  nb
Size of the second vector `b`. \n \n

\param[out] c
Pointer to the convolution output vector  \f$ c = a * b\f$. \n
Vector size is `[na + nb - 1  x  1]`. \n
Memory must be allocated. \n \n

\return `RES_OK` if convolution is calculated successfully. \n
Else \ref ERROR_CODE_GROUP "code error".

\note If vectors `a` and `b` are coefficients of two polynomials,
then convolution of the vectors `a` and `b` returns polynomial product 
coefficients.

Example:
\code{.cpp}
    complex_t ac[3] = {{0.0, 1.0}, {1.0, 1.0}, {2.0, 2.0}};
    complex_t bc[4] = {{3.0, 3.0}, {4.0, 4.0}, {5.0, 5.0}, {6.0, 6.0}};
    complex_t cc[6];
    
    int n;
    
    conv_cmplx(ac, 3, bc, 4, cc);
    
    for(n = 0; n < 6; n++)
        printf("cc[%d] = %5.1f%+5.1fj\n", n, RE(cc[n]),IM(cc[n]));
    
\endcode
 \n

Output:
\verbatim
cc[0] =  -3.0 +3.0j
cc[1] =  -4.0+10.0j
cc[2] =  -5.0+25.0j
cc[3] =  -6.0+32.0j
cc[4] =   0.0+32.0j
cc[5] =   0.0+24.0j
\endverbatim

\author Sergey Bakhurin www.dsplib.org
***************************************************************************** */
#endif
#ifdef DOXYGEN_RUSSIAN
/*! ****************************************************************************
\ingroup FILTER_CONV_GROUP
\fn int conv_cmplx(complex_t* a, int na, complex_t* b, int nb, complex_t* c) 
\brief Линейная свертка двух комплексных векторов

Функция рассчитывает линейную свертку двух векторов \f$ c = a * b\f$.

\param[in]  a
Указатель на первый вектор  \f$a\f$. \n 
Размер вектора `[na x 1]`. \n \n 

\param[in]  na
Размер первого вектора. \n  \n 

\param[in]  b
Указатель на второй вектор \f$b\f$. \n 
Размер вектора `[nb x 1]`. \n \n 

\param[in]  nb
Размер второго вектора. \n  \n 

\param[out] c
Указатель на вектор свертки \f$ c = a * b\f$. \n 
Размер вектора `[na + nb - 1  x  1]`. \n 
Память должна быть выделена. \n \n 

\return
`RES_OK` если свертка рассчитана успешно. \n 
В противном случае \ref ERROR_CODE_GROUP "код ошибки".

\note 
Если векторы `a` и `b` представляют собой коэффициенты двух полиномов,
то результат линейной свертки представляет собой коэффициенты произведения
исходных полиномов.

Пример использования функции:

\code{.cpp}
    complex_t ac[3] = {{0.0, 1.0}, {1.0, 1.0}, {2.0, 2.0}};
    complex_t bc[4] = {{3.0, 3.0}, {4.0, 4.0}, {5.0, 5.0}, {6.0, 6.0}};
    complex_t cc[6];
    
    int n;
    
    conv_cmplx(ac, 3, bc, 4, cc);
    
    for(n = 0; n < 6; n++)
        printf("cc[%d] = %5.1f%+5.1fj \n ", n, RE(cc[n]),IM(cc[n]));
    
\endcode
\n 

Результат работы:
\verbatim
cc[0] =  -3.0 +3.0j
cc[1] =  -4.0+10.0j
cc[2] =  -5.0+25.0j
cc[3] =  -6.0+32.0j
cc[4] =   0.0+32.0j
cc[5] =   0.0+24.0j
\endverbatim

\author Бахурин Сергей. www.dsplib.org
***************************************************************************** */
#endif
int DSPL_API conv_cmplx(complex_t* a, int na, complex_t* b,
                        int nb, complex_t* c)
{
    int k;
    int n;

    complex_t *t;
    size_t bufsize;

    if(!a || !b || !c)
        return ERROR_PTR;
    if(na < 1 || nb < 1)
        return ERROR_SIZE;

    bufsize = (na + nb - 1) * sizeof(complex_t);

    if((a != c) && (b != c))
        t = c;
    else
        t = (complex_t*)malloc(bufsize);

    memset(t, 0, bufsize);

    for(k = 0; k < na; k++)
    {
        for(n = 0; n < nb; n++)
                    {
            RE(t[k+n]) += CMRE(a[k], b[n]);
            IM(t[k+n]) += CMIM(a[k], b[n]);
        }
    }

    if(t!=c)
    {
        memcpy(c, t, bufsize);
        free(t);
    }

    return RES_OK;
}

