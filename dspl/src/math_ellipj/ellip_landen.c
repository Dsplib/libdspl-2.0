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
#include <math.h>
#include "dspl.h"



#ifdef DOXYGEN_ENGLISH
/*! ****************************************************************************
\ingroup SPEC_MATH_ELLIP_GROUP
\fn int  ellip_landen(double k, int n, double* y)
\brief  Function calculates complete elliptical integral 
coefficients  \f$ k_i \f$ 

Complete elliptical integral \f$ K(k) \f$ can be described as:

\f[
    K(k) = \frac{\pi}{2} \prod_{i = 1}^{\infty}(1+k_i),
\f]

here \f$ k_i \f$ -- coefficients which calculated 
iterative from \f$ k_0 = k\f$: 

\f[
    k_i = \left( \frac{k_{i-1}}{1+\sqrt{1-k_{i-1}^2}}\right)^2
\f]

This function calculates `n` fist coefficients \f$ k_i \f$, which can
be used for Complete elliptical integral.

\param[in]  k
Elliptical modulus \f$ k \f$. \n
Elliptical modulus is real parameter, which values can be from  0 to 1. \n \n

\param[in]  n
Number of \f$ k_i \f$ which need to calculate. \n 
Parameter `n` is size of output vector `y`. \n 

\param[out]  y
pointer to the real vector which keep \f$ k_i \f$. \n
Vector size is `[n x 1]`. \n
Memory must be allocated. \n \n

\return
  `RES_OK` -- successful exit, else \ref ERROR_CODE_GROUP "error code". \n
  
Example:

\include ellip_landen_test.c

Result:

\verbatim
 i        k[i]

 1    4.625e-01
 2    6.009e-02
 3    9.042e-04
 4    2.044e-07
 5    1.044e-14
 6    2.727e-29
 7    1.859e-58
 8   8.640e-117
 9   1.866e-233
10    0.000e+00
11    0.000e+00
12    0.000e+00
13    0.000e+00
\endverbatim

\note  Complete elliptical integral converges enough fast
 if modulus \f$ k<1 \f$. There are 10 to 20 coefficients \f$ k_i \f$ ​​
 are sufficient for practical applications
 to ensure  complete elliptic integral precision within EPS.

\author Sergey Bakhurin www.dsplib.org  
***************************************************************************** */
#endif
#ifdef DOXYGEN_RUSSIAN
/*! ****************************************************************************
\ingroup SPEC_MATH_ELLIP_GROUP
\fn int  ellip_landen(double k, int n, double* y)
\brief  Расчет коэффициентов \f$ k_i \f$ ряда полного эллиптического интеграла.  

Полный эллиптический интеграл \f$ K(k) \f$ может быть представлен рядом:

\f[
K(k) = \frac{\pi}{2} \prod_{i = 1}^{\infty}(1+k_i),
\f]

где \f$ k_i \f$ вычисляется итерационно при начальных условиях \f$ k_0 = k\f$: 

\f[
    k_i = \left( \frac{k_{i-1}}{1+\sqrt{1-k_{i-1}^2}}\right)^2
\f]

Данная функция рассчитывает ряд первых `n` значений \f$ k_i \f$, которые в
дальнейшем могут быть использованы для расчета эллиптического интеграла и 
эллиптических функций.

\param[in]  k
Эллиптический модуль \f$ k \f$. \n

\param[in]  n
Размер вектора `y` соответствующих коэффициентам \f$ k_i \f$.  \n \n 

\param[out]  y
Указатель на вектор значений коэффициентов \f$ k_i \f$. \n
Размер вектора `[n x 1]`. \n
Память должна быть выделена. \n \n


\return
`RES_OK` Расчет произведен успешно. \n
В противном случае \ref ERROR_CODE_GROUP "код ошибки". \n
      
Пример использования функции `ellip_landen`:

\include ellip_landen_test.c

Результат работы программы:

\verbatim
 i        k[i]

 1    4.625e-01
 2    6.009e-02
 3    9.042e-04
 4    2.044e-07
 5    1.044e-14
 6    2.727e-29
 7    1.859e-58
 8   8.640e-117
 9   1.866e-233
10    0.000e+00
11    0.000e+00
12    0.000e+00
13    0.000e+00
\endverbatim

\note
Ряд полного эллиптического интеграла сходится при значениях
эллиптического модуля \f$ k<1 \f$. При этом сходимость ряда достаточно
быстрая и для практический приложений достаточно от 10 до 20 значений 
\f$ k_i \f$ для обеспечения погрешности при расчете полного 
эллиптического интеграла в пределах машинной точности.

\author
Бахурин Сергей
www.dsplib.org  
***************************************************************************** */
#endif
int DSPL_API ellip_landen(double k, int n, double* y)
{
    int i;
    y[0] = k;

    if(!y)
        return ERROR_PTR;
    if(n < 1)
        return ERROR_SIZE;
    if(k < 0.0 || k>= 1.0)
        return ERROR_ELLIP_MODULE;

    for(i = 1; i < n; i++)
    {
        y[i] = y[i-1] / (1.0 + sqrt(1.0 - y[i-1] * y[i-1]));
        y[i] *= y[i];
    }

    return RES_OK;
}



