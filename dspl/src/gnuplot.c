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
* MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
* GNU General Public License for more details.
*
* You should have received a copy of the GNU General Public License
* along with Foobar.  If not, see <http://www.gnu.org/licenses/>.
*/
#include <stdio.h>
#include <unistd.h>
#include "dspl.h"

#define GNUPLOT_NO   1
#define GNUPLOT_WIN  2
#define GNUPLOT_PNG  3



/*******************************************************************************
Create Gnuplot header
*******************************************************************************/
int DSPL_API gnuplot_create(int argc, char* argv[], 
                            int w, int h, char* fn_png, 
                            void** hplot)
{
  FILE* hp;
  char  str[1024] = {0};
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

  switch(state)
  {
    case GNUPLOT_NO:
      hp = NULL;
      
      err = RES_OK;
      break;
      
    case GNUPLOT_WIN:
      hp = popen("gnuplot -p",  "w");
      if(!hp)
        return ERROR_GNUPLOT_CREATE;
      memset(str, 0, 1024*sizeof(char));
      sprintf(str, "set terminal wxt size %d,%d", w,h);
      gnuplot_cmd(hp, str);
      
      err =  RES_OK;
      break;
    
    case GNUPLOT_PNG:
      if(!fn_png)
        return ERROR_GNUPLOT_FNPNG;
      hp = popen("gnuplot -p",  "w");
      if(!hp)
        return ERROR_GNUPLOT_CREATE;
      memset(str, 0, 1024*sizeof(char));
      sprintf(str, 
             "set terminal pngcairo size %d,%d enhanced font 'Verdana,8'", w,h);
      gnuplot_cmd(hp, str);
      
      memset(str, 0, 1024*sizeof(char));      
      sprintf(str, "set output '%s'", fn_png); 
      gnuplot_cmd(hp, str);
      err =  RES_OK;
      break;
      
    default:
      err = ERROR_GNUPLOT_TERM;
      hp = NULL;   
  }
  *hplot = hp; 
  return err;
}





/*******************************************************************************
Write Guplot Command
*******************************************************************************/
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






/*******************************************************************************
Close Gnuplot
*******************************************************************************/
void DSPL_API gnuplot_close(void* h)
{
  if(h) 
    fclose((FILE*)h);
}

