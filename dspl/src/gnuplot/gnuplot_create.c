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
#include <unistd.h>
#include "dspl.h"




#define GNUPLOT_NO     1
#define GNUPLOT_WIN    3
#define GNUPLOT_PNG    4


#ifdef DOXYGEN_ENGLISH
/*! ****************************************************************************
\ingroup PLOT_GROUP

\brief Create GNUPLOT chart. 

This function opens the GNUPLOT package depending on the `argv` parameters 
passed to the program. After calling this function, the handle of the graph 
will be written to the address `hplot` and it becomes possible to send GNUPLOT 
commands to display the graphs.

\note From a system point of view, `hplot` is a pointer to an open file 
in which you can write commands for execution by the GNUPLOT package.

\param[in]  argc
The number of arguments to invoke the program. \n
(number of `argv` variables). \n
This value cannot be less than one, because `argv [0]`
keeps the program executable file name. \n
\n

\param[in]  argv
Pointer to a strings array of the executable file parameters . \n
The size of the array is `argc`. \n
`argv [0]` keeps the executable file name. \n
`argv [1]` sets the graph display format: \n
\verbatim

--plotwin   display the graph in a individual window
            (this parameter is used by default);

--plotpng   save the plot to a `fn_png` png file;

--noplot    plot not to create, ignore all GNUPLOT commands.
\endverbatim
\n

\param[in]  w
The width of the graph window or png file in pixels. \n 
\n

\param[in]  h
The height of the graph window or png file in pixels. \n 
\n

\param[in]  fn_png
Plot png-file name. \n 

\param[in, out] hplot
Pointer to the handle address of the GNUPLOT package. \n
A pointer to the current graph will be recorded at this address. This pointer 
is required to send GNUPLOT commands. \n 
\n

\return
`RES_OK` if function is calculated successfully. \n 
Else \ref ERROR_CODE_GROUP "code error".
The `hplot` pointer sets in `NULL` if function returns error. \n
\note
If `argv[1]` parameter equals `--noplot`, then function returns `RES_OK`, but 
`hplot` will be set as `NULL`. \n

GNUPLOT handle must be closed by \ref gnuplot_close after plotting.\n

An example of plotting sine and cosine is given in the following listing:

\include gnuplot_script_test.c

This program calculates the values of the sine and cosine functions,
and also saves the text files `dat / cosine.txt` and` dat / sine.txt`. \n

The saved data is displayed on three graphs 
(with sine, cosine data by both data):

\image html gnuplot_script_sin.png

\image html gnuplot_script_cos.png

\image html gnuplot_script_sincos.png

\author Sergey Bakhurin www.dsplib.org
***************************************************************************** */
#endif
#ifdef DOXYGEN_RUSSIAN
/*! ****************************************************************************
\ingroup PLOT_GROUP

\brief Создать график GNUPLOT.

Данная функция открывает пакет GNUPLOT в зависимости от передаваемых в программу
параметров `argv`. После вызова данной функции по адресу `hplot` будет записан
handle графика и появляется возможность посылать GNUPLOT команды для отображения
графиков. 

\note С точки зрения системы, `hplot` является указателем на открытый файл, 
 в который можно записывать команды для исполнения пакетом GNUPLOT.

\param[in]  argc
Количество аргументов вызова программы. \n
(количество переменных `argv`). \n 
Данное значение не должно быть меньше единицы, потому что `argv[0]` 
хранит имя исполняемого файла программы. \n 
\n

\param[in]  argv
Указатель на массив строк параметров исполняемого файла. \n
Размер массива `argc` строк. \n
`argv[0]` хранит имя исполняемого файла программы. \n
`argv[1]` задает формат отображения графика: \n
\verbatim

--plotwin   отображать график в отдельном окне
            (данный параметр используется по умолчанию);

--plotpng   сохранить график в png-файл, заданный именем fn_png; 

--noplot    график не создавать, игнорировать все команды GNUPLOT.
\endverbatim
\n

\param[in]  w
Ширина окна графика или png-файла в пикселях. \n 
\n

\param[in]  h
Высота окна графика или png-файла в пикселях. \n 
\n

\param[in]  fn_png
Имя png-файла, в который будет сохранен график. \n 
Путь сохранения графика совпадает с путем исполняемого файла программы. \n 

\param[in, out] hplot
Указатель на адрес хэндла пакета GNUPLOT. \n 
По данному адресу будет записан указатель на текщий график. Данный указатель 
необходим для посылки команд GNUPLOT для построения графика. \n 
\n

\return
`RES_OK` --- функция выполнена успешно. \n
В противном случае \ref ERROR_CODE_GROUP "код ошибки". \n
При возникновении ошибки по адресу  `hplot` будет записан `NULL`. \n
\note
Если в `argv[1]` задан параметр `--noplot`, то функция вернет `RES_OK`,
но по адресу  `hplot` также  будет записан `NULL`. \n

После построения графиков необходимо закрыть хэндл GNUPLOT функцией 
\ref gnuplot_close. \n

Пример построения графиков синуса и косинуса привден в следующем листинге:

\include gnuplot_script_test.c

Данная программа рассчитвает значения функций синуса и косинуса, 
а также сохрянет текстовые файлы `dat/cosine.txt` и `dat/sine.txt`. \n

Сохраненные данные отображаются на трех графиках (с данными синуса, косинуса
и обоими данными):

\image html gnuplot_script_sin.png

\image html gnuplot_script_cos.png

\image html gnuplot_script_sincos.png

\author Бахурин Сергей www.dsplib.org
***************************************************************************** */
#endif
int DSPL_API gnuplot_create(int argc, char* argv[], 
                            int w, int h, char* fn_png, 
                            void** hplot)
{
    FILE* hp;
    char    str[1024] = {0};
    int state = 0;
    int err;
    if(argc>1 && !argv)
        return ERROR_ARG_PARAM;
    
    if(argc < 2)
        state = GNUPLOT_WIN;
    else
    {
        if(!strcmp(argv[1], "--noplot"))
            state = GNUPLOT_NO;
        if(!strcmp(argv[1], "--plotwin"))
            state = GNUPLOT_WIN;
        if(!strcmp(argv[1], "--plotpng"))
            state = GNUPLOT_PNG;
    }
    if(!hplot)
        return ERROR_PTR;

    switch(state)
    {
        case GNUPLOT_NO:
            hp = NULL;
            
            err = RES_OK;
            break;
            
        case GNUPLOT_WIN:
            hp = popen("gnuplot -p",    "w");
            if(!hp)
                return ERROR_GNUPLOT_CREATE;
            memset(str, 0, 1024*sizeof(char));
            sprintf(str, "set terminal wxt size %d,%d", w,h);
            gnuplot_cmd(hp, str);
            
            err = RES_OK;
            break;
        
        case GNUPLOT_PNG:
            if(!fn_png)
                return ERROR_GNUPLOT_FNPNG;
            hp = popen("gnuplot -p",    "w");
            if(!hp)
                return ERROR_GNUPLOT_CREATE;
            memset(str, 0, 1024*sizeof(char));
            sprintf(str, "set terminal pngcairo size %d,%d\
                          enhanced font 'Verdana,8'", w, h);

            gnuplot_cmd(hp, str);
            
            memset(str, 0, 1024*sizeof(char));
            sprintf(str, "set output '%s'", fn_png); 
            gnuplot_cmd(hp, str);
            err = RES_OK;
            break;
            
        default:
            err = ERROR_GNUPLOT_TERM;
            hp = NULL;
    }
    *hplot = hp; 
    return err;
}

