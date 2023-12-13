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

 \brief Read array from the bin-file
 
 Function reads real or complex 1D/2D array size `n x m` from 
 the binary file `fn`. \n
 \n
 
 File supports 1D and 2D arrays and has follow format: \n\n
 
 \code {.unparsed}
 type     4 bytes type int.
          Can take on value:
          DAT_DOUBLE,  if x pointer to the real vector;
          DAT_COMPLEX, if x pointer to the complex vector;
 
 n        4 bytes type int.
          Number of array rows.
 
 m        4 bytes type int.
          Number of array columns.
          This parameter equals 1 because this function saves 1D vector.
 
 data     Data in binary raw.
          Data size is:
          n * sizeof(double),    if dtype==DAT_DOUBLE;
          n * sizeof(complex_t), if dtype==DAT_COMPLEX.
 \endcode
 
 Binary file can be used for algorithms verification by external packages like
 GNU Octave, Matlab, Python because the function reads data
 without loss of accuracy. \n 
 \n
 
 \param [in] fn
 File name. \n 
 \n
 
 \param [out] x
 Pointer to the data memory. \n
 Memory will be allocated for the output array. \n
 \n
 
 \param [out] pn
 Pointer to the rows number. \n 
 \n
 
 \param [out] pm
 Pointer to the columns number. \n 
 \n
 
 \param [out] dtype
 Pointer to the type of data. \n
 Can be one of follow: \n
 `DAT_DOUBLE` -- real data; \n
 `DAT_COMPLEX` -- complex data. \n 
 \n
 

 \return
 `RES_OK` if file read successfully. \n
 Else \ref ERROR_CODE_GROUP "code error".
 
 Write bin-file in GNU Octave and Matlab:
 \code {.m}
 function res = writebin(x, type, fn)
 
     if(type~=0 && type~=1)
         res = 2;
         return;
     end
 
     fid = fopen(fn, 'w');
     if(~fid)
         error('cannot to open file');
     end
 
     n = size(x, 1);
     m = size(x, 2);
         
     count = fwrite(fid,  type,  'int32');
     if(count ~= 1)
         res = 1;
         return;
     end
     count = fwrite(fid,  n,     'int32');
     if(count ~= 1)
         fclose(fid);
         res = 1;
         return;
     end
     count = fwrite(fid,  m,     'int32');
     if(count ~= 1)
         res = 1;
         fclose(fid);
         return;
     end
     
     flag = 0;
     if(type==0)
         
         count = fwrite(fid,  x,     'double');
         if(count ~= n*m)
             res = 1;
             fclose(fid);
             return;
         end
         flag = 1;
         
     else
         y = reshape(x, n*m, 1);
         z = zeros(2*n*m, 1);
         z(1:2:end) = real(y);
         z(2:2:end) = imag(y);
         count = fwrite(fid,  z,     'double');
         if(count ~= 2*n*m)
             res = 1;
             fclose(fid);
             return;
         end
         flag = 1;
     end
     if(flag == 0)
         res = 3;
     else 
         res = 0;
     end
     fclose(fid);
 end
 \endcode
 
 Example:
 
 \include readbin_test.c
 
 \author Sergey Bakhurin www.dsplib.org
***************************************************************************** */
#endif
#ifdef DOXYGEN_RUSSIAN
/*! ****************************************************************************
 \ingroup IN_OUT_GROUP

 \brief Считать данные из бинарного файла
 
 Функция считывает реальные или комплексные данные из бинарного файла 
 с именем `fn`. \n \n
 
 Файл является универсальным для хранения как одномерных,
 так и двумерных массивов и имеет следующий формат: \n 
 \n
 
 \code {.unparsed}
 type     4 байта типа int.
          Может принимать значение:
          DAT_DOUBLE,    если x указатель на вектор вещественных чисел;
          DAT_COMPLEX, если x указатель на вектор комплексных чисел.
 
 n        4 байта типа int.
          Количество строк данных.
 
 m        4 байта типа int.
          Количество столбцов данных.
          При сохранении вектора всегда равно 1.
 
 data     Данные в бинарном виде.
          Размер данных:
          n * m * sizeof(double),    если dtype==DAT_DOUBLE;
          n * m * sizeof(complex_t), если dtype==DAT_COMPLEX.
 \endcode
 
 Файл может быть использован для верификации алгоритмов сторонними пакетами,
 такими как GNU Octave, Matlab, Python и т.д. \n \n
 
 \param [in] fn
 Имя файла. \n
 \n
 
 \param [out] x
 Указатель на адрес массива данных.\n
 Данный указатель используется для возможности выделения памяти 
 под считываемые из файла данные. 
 \n
 
 \param [out] pn
 Указатель на количество строк возвращаемого массива. \n
 \n
 
 \param [out] pm
 Указатель на количество столбцов возвращаемого массива. \n
 \n
 
 
 \param[in] dtype
 Указатель на тип данных. \n
 Значение по данному адресу может быть: \n
 `DAT_DOUBLE`  -- вещественные данные; \n
 `DAT_COMPLEX` -- комплексные данные. \n
 \n
 
 
 
 \return
 `RES_OK` --- файл считан успешно. \n
 В противном случае \ref ERROR_CODE_GROUP "код ошибки". \n
 
 \note
 Данная функция производит считывание данных без потери точности,
 поэтому рекомендуется использовать ее для верификации DSPL. \n
 \n
 
 
 Функция для записи бинарного файла в GNU Octave и Matlab:
 \code {.m}
 function res = writebin(x, type, fn)
 
     if(type~=0 && type~=1)
         res = 2;
         return;
     end
 
     fid = fopen(fn, 'w');
     if(~fid)
         error('cannot to open file');
     end
 
     n = size(x, 1);
     m = size(x, 2);
         
     count = fwrite(fid,  type,  'int32');
     if(count ~= 1)
         res = 1;
         return;
     end
     count = fwrite(fid,  n,     'int32');
     if(count ~= 1)
         fclose(fid);
         res = 1;
         return;
     end
     count = fwrite(fid,  m,     'int32');
     if(count ~= 1)
         res = 1;
         fclose(fid);
         return;
     end
     
     flag = 0;
     if(type==0)
         
         count = fwrite(fid,  x,     'double');
         if(count ~= n*m)
             res = 1;
             fclose(fid);
             return;
         end
         flag = 1;
         
     else
         y = reshape(x, n*m, 1);
         z = zeros(2*n*m, 1);
         z(1:2:end) = real(y);
         z(2:2:end) = imag(y);
         count = fwrite(fid,  z,     'double');
         if(count ~= 2*n*m)
             res = 1;
             fclose(fid);
             return;
         end
         flag = 1;
     end
     if(flag == 0)
         res = 3;
     else 
         res = 0;
     end
     fclose(fid);
 end
 \endcode
 
 Пример записи данных в бинарный файл с последующим считыванием функцией 
 \ref writebin :
 
 \include readbin_test.c
  
 \author Бахурин Сергей www.dsplib.org
***************************************************************************** */
#endif
int DSPL_API readbin(char* fn, void** x, int* pn, int* pm, int* dtype)
{
    FILE* pFile = NULL;
    int n, m, t, err;

    if(!x || !pn || !pm)
        return ERROR_PTR;
    if(!fn)
        return ERROR_FNAME;

    pFile = fopen(fn, "rb");
    if(pFile == NULL)
            return ERROR_FOPEN;

    if(fread(&t, sizeof(int), 1, pFile) != 1)
    {
        err = ERROR_FREAD_SIZE;
        goto exit_label;
    }

    if(dtype)
        *dtype = t;

    if(fread(&n, sizeof(int), 1, pFile) != 1)
    {
        err = ERROR_FREAD_SIZE;
        goto exit_label;
    }
    *pn = n;

    if(fread(&m, sizeof(int), 1, pFile) != 1)
    {
        err = ERROR_FREAD_SIZE;
        goto exit_label;
    }
    *pm = m;

    switch(t)
    {
        case DAT_DOUBLE:
            (*x) = (*x) ? realloc(*x, n*m*sizeof(double)) :
                          malloc(n*m*sizeof(double));
            if(fread(*x, sizeof(double), n*m, pFile) != n*m)
            {
                err = ERROR_FREAD_SIZE;
                goto exit_label;
            }

            break;
        case DAT_COMPLEX:
            (*x) = (*x) ? realloc(*x, n*m*sizeof(complex_t)) :
                          malloc(n*m*sizeof(complex_t));
            if(fread(*x, sizeof(complex_t), n*m, pFile) != n*m)
            {
                err = ERROR_FREAD_SIZE;
                goto exit_label;
            }
            break;
        default:
            err = ERROR_DAT_TYPE;
            goto exit_label;
    }
    err = RES_OK;
exit_label:
    if(pFile)
        fclose(pFile);
    return err;
}

