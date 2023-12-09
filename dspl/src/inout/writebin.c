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
 \fn    int writebin(void* x, int n, int m, int dtype, char* fn)
 \brief Save 1D vector or 2D matrix to the binary file
 
 Function saves real or complex 1D/2D data size `n x m` to the 
 binary file `fn`. \n
 \n
 
 File format supports 1D and 2D arrays and has follow format: \n
 \n
 
 \code {.unparsed}
 type     4 bytes type int.
          Can take on value:
          DAT_DOUBLE,  if x pointer to the real data;
          DAT_COMPLEX, if x pointer to the complex data;
 
 n        4 bytes type int.
          Number of array rows.
 
 m        4 bytes type int.
          Number of array columns.
          This parameter equals 1 means function saves 1D vector.
 
 data     Data in binary raw.
          Data size is:
          n * sizeof(double),    if dtype==DAT_DOUBLE;
          n * sizeof(complex_t), if dtype==DAT_COMPLEX.
 \endcode 
 
 Binary file can be used for algorithms verification by external packages like
 GNU Octave, Matlab, Python because the function writes to a file
 without loss of accuracy. \n 
 \n
 
 \param [in] x
 Pointer to the input vector. \n
 Vector size is `[n x 1]`. \n 
 \n
 
 \param [in] n
 Rows number. \n 
 \n
 
 \param [in] m
 Columns number. \n 
 \n
 
 \param [in] dtype
 Type of data. \n
 Can be one of follow: \n
 `DAT_DOUBLE` -- real data; \n
 `DAT_COMPLEX` -- complex data. \n 
 \n
 
 \param [in] fn
 File name. \n 
 \n
 
 \return
 `RES_OK` if file is saved successfully. \n
 Else \ref ERROR_CODE_GROUP "code error".
 
 Reading saved binary file from GNU Octave or Matlab:
 \code{.m}
 function [dat, n, m] = readbin(fn)
     fid = fopen(fn);
     if(~fid)
         error('cannot to open file');
     end
     type = fread(fid, 1, 'int32');
     n    = fread(fid, 1, 'int32');
     m    = fread(fid, 1, 'int32');
 
     if(type==0)
         dat = fread(fid, [n*m, 1], 'double');
     end
 
     if(type==1)
         y = fread(fid, [n*m*2, 1], 'double');
         dat = y(1:2:end) + 1i * y(2:2:end);
     end
 
     dat = reshape(dat, n, m);
 
     fclose(fid);
 end
 \endcode
 
 \author Sergey Bakhurin www.dsplib.org
***************************************************************************** */
#endif
#ifdef DOXYGEN_RUSSIAN
/*! ****************************************************************************
 \ingroup IN_OUT_GROUP
 \fn    int writebin(void* x, int n, int m, int dtype, char* fn)
 \brief Сохранить данные в бинарный файл
 
 Функция сохраняет реальный или комплексный массив данных
 размера `[n x m]` в бинарный файл `fn`. \n 
 \n
 
 Файл является универсальным для хранения как одномерных,
 так и двумерных массивов и имеет следующий формат: \n 
 \n
 
 \code {.unparsed}
 type     4 байта типа int.
          Может принимать значение:
          DAT_DOUBLE,  если x указатель на вектор вещественных чисел;
          DAT_COMPLEX, если x указатель на вектор комплексных чисел.
 
 n        4 байта типа int.
          Количество строк данных.
 
 m        4 байта типа int.
          Количество столбцов данных.

 data     Данные в бинарном виде.
          Размер данных:
          n * sizeof(double),    если dtype==DAT_DOUBLE;
          n * sizeof(complex_t), если dtype==DAT_COMPLEX.
 \endcode
 
 Файл может быть использован для верификации алгоритмов сторонними пакетами,
 такими как GNU Octave, Matlab, Python и т.д. \n 
 \n
 
 \param [in] x
 Указатель на массив данных. \n
 Размер вектора `[n x 1]`. \n
 \n
 
 \param [in] n
 Количество строк данных. \n
 \n
 
 \param [in] m
 Количество столбцов данных. \n
 \n
 
 
 \param [in] dtype
 Тип данных. \n
 Может принимать значения: \n
 `DAT_DOUBLE`  -- вещественные данные; \n
 `DAT_COMPLEX` -- комплексные данные. \n
 \n
 
 \param [in] fn
 Имя файла. \n
 \n
 
 \return
 `RES_OK` --- файл сохранен успешно. \n
 В противном случае \ref ERROR_CODE_GROUP "код ошибки". \n
 
 \note
 Данная функция производит запись в файл без потери точности,
 поэтому рекомендуется использовать ее для верификации данных DSPL. \n
 \n
 
 
 Функция для чтения бинарного файла в GNU Octave и Matlab:
 \code {.m}
 function [dat, n, m] = readbin(fn)
     fid = fopen(fn);
     if(~fid)
         error('cannot to open file');
     end
     type = fread(fid, 1, 'int32');
     n    = fread(fid, 1, 'int32');
     m    = fread(fid, 1, 'int32');
 
     if(type==0)
         dat = fread(fid, [n*m, 1], 'double');
     end
 
     if(type==1)
         y = fread(fid, [n*m*2, 1], 'double');
         dat = y(1:2:end) + 1i * y(2:2:end);
     end
 
     dat = reshape(dat, n, m);
 
     fclose(fid);
 end
 \endcode
 
 \author Бахурин Сергей www.dsplib.org
***************************************************************************** */
#endif
int DSPL_API writebin(void* x, int n, int m, int dtype, char* fn)
{
    int res;
    FILE* pFile = NULL;

    if(!x)
        return ERROR_PTR;
    if(n < 1)
        return ERROR_SIZE;
    if(!fn)
        return    ERROR_FNAME;

    pFile = fopen(fn, "wb");
    if(pFile == NULL)
            return ERROR_FOPEN;


    if(fwrite(&dtype, sizeof(int), 1, pFile) != 1)
    {
        res = ERROR_FWRITE_SIZE;
        goto exit_label;
    }



    if(fwrite(&n, sizeof(int), 1, pFile) != 1)
    {
        res = ERROR_FWRITE_SIZE;
        goto exit_label;
    }

    if(fwrite(&m, sizeof(int), 1, pFile) != 1)
    {
        res = ERROR_FWRITE_SIZE;
        goto exit_label;
    };

    switch(dtype)
    {
        case DAT_DOUBLE:
            if(fwrite((double*)x, sizeof(double), n*m, pFile) != n)
            {
                res = ERROR_FWRITE_SIZE;
                goto exit_label;
            }
            break;
        case DAT_COMPLEX:
            if(fwrite((complex_t*)x,
                    sizeof(complex_t), n*m, pFile) != n)
            {
                res = ERROR_FWRITE_SIZE;
                goto exit_label;
            }
            break;
        default:
            res = ERROR_DAT_TYPE;
            goto exit_label;
    }
    res = RES_OK;
exit_label:
    if(pFile)
        fclose(pFile);
    return res;
}

