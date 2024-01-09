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
#include "contour2d.h"

int DSPL_API contour2d(double* z, double* x, double* y, 
                       int n, int m, double lev, 
                       contour2d_t* c)
{
     
    linseg_t *ls = NULL;
    int nls, i;
    
    // Создаем полный массив отезков
    linseg_create(z, x, y, n, m, lev, &ls, &nls);
    
	//for(i =0; i< nls; i++)
    //    printf("%d, [%.1f %.1f] -- [%.1f %.1f]\n", i, ls[i].p[0][0], ls[i].p[0][1], ls[i].p[1][0], ls[i].p[1][1]);
    
    // очищаем контур при необходимости
	contour2d_free(c);

    // Создаем массив линий
    if (c->lines == NULL)
    {
        c->nlines = BSIZE;
        c->lines = (line2d_t*)malloc(c->nlines*sizeof(line2d_t));
        memset(c->lines, 0, c->nlines*sizeof(line2d_t));
    }
    
    // Главный цикл. На каждой итерации создаем 
    // одну линию и добавляем ее в контур
    i = -1;
    do
    {
        i++;
        // при необходимости увеличиваем размер массива
        if(i >= c->nlines)
        {
            
            c->lines = (line2d_t*)realloc(c->lines, (c->nlines + BSIZE)*sizeof(line2d_t));
			memset(c->lines + c->nlines, 0, BSIZE*sizeof(line2d_t));
			c->nlines += BSIZE;
			
        }
        // создаем очередную линию путем стыковки отезков
        line_create(ls, nls, c->lines+i);
        //for(j =0; j< c->lines[i].npoints; j++)
        //    printf("[%.1f %.1f] -- ", c->lines[i].points[j][0], c->lines[i].points[j][1]);
        //printf("\n");
    }while(c->lines[i].npoints);
    
    c->lines = (line2d_t*)realloc(c->lines, i*sizeof(line2d_t));
    c->nlines = i;
    c->level = lev;
    
    if(ls)
        free(ls);
    return RES_OK;
}




/*
int line_create(linseg_t* ls, int nls, point2d_t** line, int* np)
{
    int i, j, c, n, k;
    if(!line || !ls ||  !np)
        return ERROR_PTR;
    // printf("*line = %x, *np = %d \n", *line, *np);
    i = 0;
    while(!(ls[i].flag) && i < nls)
        i++;
    if(i==nls)
    {
        *np = 0;
        return RES_OK;
    }
    //printf("i = %d, ls[%d] = [%.1f %.1f] -- [%.1f, %.1f]\n", 
    //        i, i,   ls[i].p[0][0], ls[i].p[0][1], 
    //                ls[i].p[1][0], ls[i].p[1][1]);
    n = BSIZE;
    if((*line) == NULL)
        (*line) = (point2d_t*)malloc(n*sizeof(point2d_t));
    else
        (*line) = (point2d_t*)realloc((*line), n*sizeof(point2d_t));

    c = 0;
    (*line)[c][0] = ls[i].p[0][0];
    (*line)[c][1] = ls[i].p[0][1];
    c++;
    (*line)[c][0] = ls[i].p[1][0];
    (*line)[c][1] = ls[i].p[1][1];
    c++;
    ls[i].flag = 0;
    
    for(i = 0; i < nls; i++)
    {
        for(j = 0; j < nls; j++)
        {
            if(ls[j].flag)
            {
                //сравниваем последнюю точку линии с первой точкой отрезка ls[j]
                if(dist((*line)+c-1, ls[j].p) < DIST_TH)
                {
        
                    // проверяем выделение памяти.
                    // при необходимости увеличиваем
                    if(c>=n)
                    {
                        n += BSIZE;
                        (*line) = (point2d_t*)realloc(*line, n*sizeof(point2d_t));
                    }
                    //если первая точка отрезка совпадает, то добавляем
                    //в линию вторую точку СЗАДИ
                    (*line)[c][0] = ls[j].p[1][0];
                    (*line)[c][1] = ls[j].p[1][1];
                    ls[j].flag = 0;
                    c++;
                }
            }
            if(ls[j].flag)
            {
                // сравниваем последнюю точку линии 
                // со второй точкой отрезка ls[j]
                if(dist((*line)+c-1, ls[j].p+1) < DIST_TH)
                {
                    // проверяем выделение памяти.
                    // при необходимости увеличиваем
                    if(c>=n)
                    {
                        n += BSIZE;
                        (*line) = (point2d_t*)realloc(*line, n*sizeof(point2d_t));
                    }
                    //если вторая точка совпадает, то добавляем
                    //в линию первую точку СЗАДИ
                    (*line)[c][0] = ls[j].p[0][0];
                    (*line)[c][1] = ls[j].p[0][1];
                    ls[j].flag = 0;
                    c++;
                }
            }
            
            if(ls[j].flag)
            {
                //сравниваем последнюю точку линии с первой точкой отрезка ls[j]
                if(dist((*line), ls[j].p) < DIST_TH)
                {
        
                    // проверяем выделение памяти.
                    // при необходимости увеличиваем
                    if(c>=n)
                    {
                        n += BSIZE;
                        (*line) = (point2d_t*)realloc(*line, n*sizeof(point2d_t));
                    }
                    //если первая точка отрезка совпадает, то добавляем
                    //в линию вторую точку СПЕРЕДИ
                    for(k = c; k>0; k--)
                    {
                        (*line)[k][0] = (*line)[k-1][0];
                        (*line)[k][1] = (*line)[k-1][1]; 
                    }
                    (*line)[0][0] = ls[j].p[1][0];
                    (*line)[0][1] = ls[j].p[1][1];
                    ls[j].flag = 0;
                    c++;
                }
            }
            if(ls[j].flag)
            {
                // сравниваем последнюю точку линии 
                // со второй точкой отрезка ls[j]
                if(dist((*line), ls[j].p+1) < DIST_TH)
                {
                    // проверяем выделение памяти.
                    // при необходимости увеличиваем
                    if(c>=n)
                    {
                        n += BSIZE;
                        (*line) = (point2d_t*)realloc(*line, n*sizeof(point2d_t));
                    }
                    //если вторая точка совпадает, то добавляем
                    //в линию первую точку СПЕРЕДИ
                    
                    for(k = c; k>0; k--)
                    {
                        (*line)[k][0] = (*line)[k-1][0];
                        (*line)[k][1] = (*line)[k-1][1]; 
                    }
                    
                    (*line)[0][0] = ls[j].p[0][0];
                    (*line)[0][1] = ls[j].p[0][1];
                    ls[j].flag = 0;
                    c++;
                }
            }
            
        }
    }
    (*line) = (point2d_t*)realloc((*line), c*sizeof(point2d_t));
    *np = c;
    return RES_OK;
}

*/



