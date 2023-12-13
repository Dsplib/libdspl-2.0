/*
* Copyright (c) 2015-2024 Sergey Bakhurin
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
\ingroup IN_OUT_GROUP

\brief Save real data \f$y(x)\f$ to the text file `fn`. \n

File format is:

\verbatim
x[0]        y[0]
x[1]        y[1]
...         ...
x[n-1]    y[n-1]
\endverbatim

Text file can be used to plotting data with a third-party program
for example, the GNUPLOT package (see \ref PLOT_GROUP). \n

\param[in] x
Pointer to the vector `x`.    \n
Vector size is `[n x 1]`. \n \n

\param[in] y
Pointer to the vector `y`. \n
Vector size is    `[n x 1]`. \n
This pointer can be `NULL`.
File will have only one column corresponds to `x` vector in this case. \n \n

\param[in] n
Size of vectors `x` and `y`. \n \n

\param[in] fn
File name. \n \n


\return
`RES_OK` if file is saved successfully. \n
Else \ref ERROR_CODE_GROUP "code error".

\note
This function rounds data when writing to a file.
So, it is not recommended to use it to verify algorithms.

\author Sergey Bakhurin www.dsplib.org
***************************************************************************** */
#endif
#ifdef DOXYGEN_RUSSIAN
/*! ****************************************************************************
\ingroup IN_OUT_GROUP

\brief Сохранить вещественные данные в текстовый файл

Функция сохраняет вещественные данные в текстовый файл `fn`. \n
Файл имеет следующий формат: \n
\verbatim
x[0]        y[0]
x[1]        y[1]
...         ...
x[n-1]    y[n-1]
\endverbatim

Файл может быть использован для построения графика сторонней программой,
например пакетом GNUPLOT (см. раздел \ref PLOT_GROUP). \n

\param[in] x
Указатель на первый вектор. \n
Размер вектора `[n x 1]`. \n \n

\param[in] y
Указатель на второй вектор. \n
Размер вектора `[n x 1]`. \n
Может быть `NULL`. \n
Файл будет содержать только один столбец соответствующий
вектору `x` если `y == NULL`. \n \n

\param[in] n
Размер входных векторов. \n \n

\param[in] fn
Имя файла. \n \n

\return
`RES_OK` --- файл сохранен успешно. \n
В противном случае \ref ERROR_CODE_GROUP "код ошибки": \n

\note
Данная функция производит округление данных при записи в файл.
Поэтому не рекомендуется использовать ее для верификации данных DSPL.

\author Бахурин Сергей www.dsplib.org
***************************************************************************** */
#endif
int DSPL_API writetxt(double* x, double* y, int n, char* fn)
{
    int k, err;
    FILE* pFile = NULL;
    if(!x)
        return ERROR_PTR;
    if(n < 1)
        return ERROR_SIZE;
    if(!fn)
        return ERROR_FNAME;

    pFile = fopen(fn, "w+");
    if(pFile == NULL)
        return ERROR_FOPEN;

    if(y)
    {
        for(k = 0; k < n; k++)
        {
            if(isnan(x[k]) || isnan(y[k]))
            {
                err = ERROR_NAN;
                goto exit_label;
            }
            if(isinf(x[k]) || isinf(y[k]))
            {
                err = ERROR_INF;
                goto exit_label;
            }
            fprintf(pFile, "%+.12E\t%+.12E\n", x[k], y[k]);
        }
    }
    else
    {
        for(k = 0; k < n; k++)
        {
            if(isnan(x[k]))
            {
                err = ERROR_NAN;
                goto exit_label;
            }
            if(isinf(x[k]))
            {
                err = ERROR_INF;
                goto exit_label;
            }
            fprintf(pFile, "%+.12E\n", x[k]);
        }
    }
    err = RES_OK;
exit_label:
    fclose(pFile);
    return err;
}
