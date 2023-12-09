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



#ifdef DOXYGEN_ENGLISH
/*! ****************************************************************************
\ingroup PLOT_GROUP
\fn void gnuplot_cmd(void* h, char* cmd)
\brief Function sends `cmd` command to GNUPLOT corresponds to `h` handle.

\note From a system point of view, `h` is a pointer to an open file,
into which you can write commands for execution by the GNUPLOT package. This
 function writes the string `cmd` to the given file. \n \n
 
From the user's point of view, calling the `gnuplot_cmd` function is 
GNUPLOT script one line executing. \n \n

A comprehensive description of the graphing capabilities of the GNUPLOT package,
with examples of commands can be found on
  <a href = "http://gnuplot.sourceforge.net/demo_5.2/"> GNUPLOT project </a>.

\param[in]  h
GNUPLOT handle. \n \n

\param[in]  cmd
GNUPLOT command string. \n \n

\author Sergey Bakhurin www.dsplib.org
***************************************************************************** */
#endif
#ifdef DOXYGEN_RUSSIAN
/*! ****************************************************************************
\ingroup PLOT_GROUP
\fn void gnuplot_cmd(void* h, char* cmd)
\brief Функция посылает команду `cmd` пакету GNUPLOT, для построения 
или оформления графика, соответсвующего хэндлу `h`.

Данная функция закрывает открытый ранее хэндл GNUPLOT после передачи всех 
команд построения графиков.

\note С точки зрения системы, `h` является указателем на открытый файл, 
в который можно записывать команды для исполнения пакетом GNUPLOT. Данная 
функция записывает в данный файл строку  `cmd`. \n \n

С точки зрения пользователя, вызов функции `gnuplot_cmd` равносильно выполнению
одной строки скрипта GNUPLOT. \n \n

Исчерпывающее описание возможностей построения графиков пакетом GNUPLOT, 
с примерами команд можно найти на 
 <a href = "http://gnuplot.sourceforge.net/demo_5.2/">сайте проекта GNUPLOT</a>. 

\param[in]  h
Хэндл графика GNUPLOT. \n \n

\param[in]  cmd
Тектовая строка команды построения или оформления графика. \n \n

\author Бахурин Сергей www.dsplib.org
***************************************************************************** */
#endif
void DSPL_API gnuplot_cmd(void* h, char* cmd)
{
    if(h)
    {
        FILE *pf = (FILE*)h;
        fprintf(pf, cmd);
        fflush(pf);
        fprintf(pf, "\n");
        fflush(pf);
    }
}

