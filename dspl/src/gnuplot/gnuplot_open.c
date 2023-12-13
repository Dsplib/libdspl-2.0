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
#include <unistd.h>
#include "dspl.h"




#ifdef DOXYGEN_ENGLISH
/*! ****************************************************************************
\ingroup PLOT_GROUP

\brief Open GNUPLOT program. 

This function opens the GNUPLOT package. 
After calling this function, the handle of the GNUPLOT 
will be written to the address `hplot` and it becomes possible to send GNUPLOT 
commands.

\note From a system point of view, `hplot` is a pointer to an open file 
in which you can write commands for execution by the GNUPLOT package.

\param[in, out] hplot
Pointer to the handle address of the GNUPLOT package. \n
This pointer is required to send GNUPLOT commands. \n 
\n

\return
`RES_OK` if function is calculated successfully. \n 
Else \ref ERROR_CODE_GROUP "code error".
The `hplot` pointer sets in `NULL` if function returns error. \n

GNUPLOT handle must be closed by \ref gnuplot_close after plotting.\n

An example of plotting sine and cosine is given in the following listing:

\note
The difference between `gnuplot_open` and` gnuplot_create` is that 
`gnuplot_create` processes the program execution parameters and 
creates a GNUPLOT terminal. \n
The gnuplot_open function opens the GNUPLOT handle to be able to 
send commands, regardless of program execution parameters, 
but does not create terminals.

\author Sergey Bakhurin www.dsplib.org
***************************************************************************** */
#endif
#ifdef DOXYGEN_RUSSIAN
/*! ****************************************************************************
\ingroup PLOT_GROUP

\brief Открыть пакет GNUPLOT.

Данная функция открывает пакет GNUPLOT. 
После вызова данной функции по адресу `hplot` будет записан
handle GNUPLOT и появляется возможность посылать GNUPLOT команды. 

\note С точки зрения системы, `hplot` является указателем на открытый файл, 
 в который можно записывать команды для исполнения пакетом GNUPLOT.

\param[in, out] hplot
Указатель на адрес хэндла пакета GNUPLOT. \n 
По данному адресу будет записан указатель на текщий график. Данный указатель 
необходим для посылки команд GNUPLOT для построения графика. \n 
\n

\return
`RES_OK` --- функция выполнена успешно. \n
В противном случае \ref ERROR_CODE_GROUP "код ошибки". \n
При возникновении ошибки по адресу  `hplot` будет записан `NULL`. \n

После построения графиков необходимо закрыть хэндл GNUPLOT функцией 
\ref gnuplot_close. \n

\note
Отличие функции `gnuplot_open` от `gnuplot_create` заключается в том, 
что `gnuplot_create` обрабатывает параметры выполнения программы 
и создает терминал GNUPLOT.\n
Функция `gnuplot_open` открывает GNUPLOT хэндл для возможности посылки команд, 
вне зависимости от параметров выполнения программы, но не создает терминалов.

\author Бахурин Сергей www.dsplib.org
***************************************************************************** */
#endif
int DSPL_API gnuplot_open(void** hplot)
{
    if(!hplot)
        return ERROR_PTR;
    *hplot = popen("gnuplot -p", "w");

    if(!(*hplot))
        return ERROR_GNUPLOT_CREATE;
    return RES_OK;
}




