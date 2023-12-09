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
\fn void gnuplot_close(void* h) 
\brief Close GNUPLOT handle. 

This function closes the previously opened GNUPLOT handle after transferring all
commands.

\note 
From a system point of view, `h` is a pointer to an open file in which
You can write commands for execution by the GNUPLOT package.
This function closes this file. \n
It is important to note that closing the handle does not mean closing the 
plot window.
After closing the package handle, the user cannot send build commands, 
but the plot window remains open, as it is processed by the package
GNUPLOT regardless.

\param[in]  h
GNUPLOT handle which will be closed. \n 
\n

\author Sergey Bakhurin www.dsplib.org
***************************************************************************** */
#endif
#ifdef DOXYGEN_RUSSIAN
/*! ****************************************************************************
\ingroup PLOT_GROUP
\fn void gnuplot_close(void* h) 
\brief Закрыть хэндл GNUPLOT. 

Данная функция закрывает открытый ранее хэндл GNUPLOT после передачи всех 
команд построения графиков.

\note 
С точки зрения системы, `h` является указателем на открытый файл, в который 
можно записывать команды для исполнения пакетом GNUPLOT. 
Данная функция закрывает этот файл. \n
Важно отметить, что закрытие хэндла не означает закрытия окна графика. 
После закрытия хэндла пакета, пользователь не может посылать команды построения
графика, но окно графика остается открытым, так как обрабатывается пакетом 
GNUPLOT независимо.

\param[in]  h
Хэндл графика GNUPLOT, который будет закрыт. \n 
\n

\author Бахурин Сергей 
www.dsplib.org
***************************************************************************** */
#endif
void DSPL_API gnuplot_close(void* h)
{
    if(h) 
        fclose((FILE*)h);
}

