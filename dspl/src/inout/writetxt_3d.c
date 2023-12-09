/*
* Copyright (c) 2015-2022 Sergey Bakhurin
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
\fn    int writetxt_3d(double* x, int nx, double* y, int ny,
                       double* z, char* fn)

\brief Save data \f$ z(x,y)\f$    to the text file `fn` for 3D surface plotting.

Function \f$ z(x,y)\f$ describes as matrix `z[x[n], y[n]]` as it showed on the
follow figure:

\image html writetxt_3d_matrix.png

Matrix `z` writen in the memory by columns as it it showed on the figure by
red arrow. \n

Text file fas follow format: \n

\verbatim
x[0]        y[0]        z[0, 0]
x[1]        y[0]        z[1, 0]
x[2]        y[0]        z[2, 0]
...         ...         ...
x[nx-1] y[0]        z[nx-1, 0]

x[0]        y[1]        z[0, 1]
x[1]        y[1]        z[1, 1]
x[2]        y[1]        z[2, 1]
...         ...         ...
x[nx-1] y[1]        z[nx-1, 1]

...         ...         ...
...         ...         ...
...         ...         ...

x[0]        y[ny-1] z[0, ny-1]
x[1]        y[ny-1] z[1, ny-1]
x[2]        y[ny-1] z[2, ny-1]
...         ...         ...
x[nx-1] y[ny-1] z[nx-1, ny-1]
\endverbatim
Each `z` matrix value is writen on individual line corresponds to `x`
and `y` values. Matrix columns are separated from each other by an empty line.

The file can be used to build a 3D surface with a third-party program
for example, the GNUPLOT package (see \ref PLOT_GROUP). Also this format
supported by pgfplot3d packages of the Latex system. \n

\param[in] x
Pointer to the vector `x`. \n
Vector size is `[nx x 1]`. \n \n

\param[in] nx
Size of vector `x`. \n \n


\param[in] y
Pointer to the vector `y`. \n
Vector size is    `[ny x 1]`. \n\n

\param[in] ny
Size of vector `y`. \n\n

\param[in] z
Pointer to the matrix `z(x, y)`. \n
Size of matrix is    `[nx x ny]`. \n\n

\param[in] fn
3D data file name. \n\n


\return
`RES_OK` if file is saved successfully. \n
Else \ref ERROR_CODE_GROUP "code error".

Example of 3D surface plotting:

\include writetxt_3d_test.c

Program calcultes function

\f[
z(x,y) = x \exp(-x^2 -y^2)
\f]

and save data to the `dat/data3d.txt` file.\n
In addition, GNUPLOT built a 3D surface by data
 saved to the `dat/data3d.txt` file:

\image html writetxt_3d.png

\author Sergey Bakhurin www.dsplib.org
***************************************************************************** */
#endif
#ifdef DOXYGEN_RUSSIAN
/*! ****************************************************************************
\ingroup IN_OUT_GROUP
\fn int writetxt_3d(double* x, int nx, double* y, int ny, double* z, char* fn)

\brief Сохранить данные для построения 3D графика

Функция сохраняет вещественные данные в текстовый файл `fn`
для построения 3D поверхности. \n
Функция \f$ z(x,y)\f$ описывается матрицей значений `z[x[n], y[n]]`
как это показано на рисунке

\image html writetxt_3d_matrix.png

Матрица `z` хранится в памяти по столбцам, как это показано красной стрелкой. \n
Файл имеет следующий формат: \n

\verbatim
x[0]        y[0]        z[0, 0]
x[1]        y[0]        z[1, 0]
x[2]        y[0]        z[2, 0]
...         ...         ...
x[nx-1] y[0]        z[nx-1, 0]

x[0]        y[1]        z[0, 1]
x[1]        y[1]        z[1, 1]
x[2]        y[1]        z[2, 1]
...         ...         ...
x[nx-1] y[1]        z[nx-1, 1]

...         ...         ...
...         ...         ...
...         ...         ...

x[0]        y[ny-1] z[0, ny-1]
x[1]        y[ny-1] z[1, ny-1]
x[2]        y[ny-1] z[2, ny-1]
...         ...         ...
x[nx-1] y[ny-1] z[nx-1, ny-1]
\endverbatim
Таким образом, каждое значение матрицы `z` записано отдельной строкой со
соответствующими значениями `x` и `y`. Столбцы матрицы отделены пустой строкой.

Файл может быть использован для построения графика сторонней программой,
например пакетом GNUPLOT (см. раздел \ref PLOT_GROUP). Также данный формат
поддерживается пакетами pgfplot3d издательской системы Latex. \n

\param[in] x
Указатель на вектор значений оси `x`. \n
Размер вектора `[nx x 1]`. \n \n

\param[in] nx
Размер вектора оси `x`. \n \n


\param[in] y
Указатель на второй вектор значений оси `y`. \n
Размер вектора    `[ny x 1]`. \n

\param[in] ny
Размер вектора оси `y`. \n \n

\param[in] z
Указатель на матрицу значений функции `z(x, y)`. \n \n

\param[in] fn
Имя файла. \n \n

\return
`RES_OK` --- файл сохранен успешно. \n
В противном случае \ref ERROR_CODE_GROUP "код ошибки": \n

Пример использования функции и посторения 3D поверхности приведен
в следующем листинге:

\include writetxt_3d_test.c

Данная программа рассчитывает и строит поверхность функции

\f[
z(x,y) = x \exp(-x^2 -y^2)
\f]

В каталоге `dat` будет создан файл `data3d.txt`.\n
Кроме того программа GNUPLOT произведет построение 3D поверхности
по сохраненным в файл данным:

\image html writetxt_3d.png

\author Бахурин Сергей www.dsplib.org
***************************************************************************** */
#endif
int DSPL_API writetxt_3d(double* x, int nx, double* y, int ny,
                         double* z, char* fn)
{
    int k, n;
    FILE* pFile = NULL;

    if(!x || !y || !z)
        return ERROR_PTR;
    if(nx < 1 || ny < 1)
        return ERROR_SIZE;
    if(!fn)
        return ERROR_FNAME;

    pFile = fopen(fn, "w+");
    if(pFile == NULL)
        return ERROR_FOPEN;

    for(k = 0; k < ny; k++)
    {
        for(n = 0; n < nx; n++)
        {
            if(!isnan(z[n+k*nx]))
            {
                fprintf(pFile, "%+.12E\t%+.12E\t%+.12E\n",
                               x[n], y[k], z[n+k*nx]);

            }
        }
        fprintf(pFile, "\n");

    }
    fclose(pFile);
    return RES_OK;
}

