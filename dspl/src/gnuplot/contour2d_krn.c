 
 
 
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "dspl.h"
#include "contour2d.h"

/*******************************************************************************
Добавляем отрезок в массив отрезков
*******************************************************************************/
int add_linseg(linseg_t** ls, int* lsnum, int* lscnt, 
                point2d_t* p0, point2d_t* p1)

{
    int n, c;
    n = *lsnum;
    c = *lscnt;
    // проверяем выделение памяти.
    if((n == 0) && ((*ls)==NULL))
    {
        n = BSIZE;
        (*ls) = (linseg_t*)malloc(n * sizeof(linseg_t));
    }
    else
    {
        // при необходимости увеличиваем
        if(c >= n)
        {
            n += BSIZE;
            (*ls) = (linseg_t*)realloc((*ls), n * sizeof(linseg_t));
        }
    }
    (*ls)[c].p[0][0] = p0[0][0];
    (*ls)[c].p[0][1] = p0[0][1];
    (*ls)[c].p[1][0] = p1[0][0];
    (*ls)[c].p[1][1] = p1[0][1];
    (*ls)[c].flag = 1;
    c++;
    (*lsnum) = n;
    (*lscnt) = c;
    
    return RES_OK;
}

/*******************************************************************************
Создаем массив отрезков
*******************************************************************************/
int linseg_create(double* z, double* x, double* y, 
                  int n, int m, double lev, 
                  linseg_t** ls, int* sz)
{
    int lsnum, lscnt, t, in, im, i;
    point2d_t p0 = {0};
    point2d_t p1 = {0};

    double dx;
    double dy;
    
    if((ls== NULL)||(z==NULL))
        return ERROR_PTR;
    lsnum = 0;
    lscnt = 0;
    
    for(in = 0; in < n-1; in++)
    {
        for(im = 0; im < m-1; im++)
        {
            i = in + im * n;
            t = 0;
            t += z[i]     > lev ? 8 : 0;
            t += z[i+n]   > lev ? 4 : 0;
            t += z[i+n+1] > lev ? 2 : 0;
            t += z[i+1]   > lev ? 1 : 0;
            
            //printf("%d, %d, %d\n", in, im, t);
            switch(t)
            {
                case 0:
                case 15:
                    break;
                case 1:
                case 14:
                    dx = (lev - z[i]) / (z[i+1] - z[i]);
                    p0[0] = x[in] + dx * (x[in+1] - x[in]);
                    p0[1] = y[im];
                    
                    dy = (lev - z[i+1]) / (z[i+n+1] - z[i+1]);
                    p1[0] = x[in+1];
                    p1[1] = y[im] + dy * (y[im+1] - y[im]);
                    add_linseg(ls, &lsnum, &lscnt, &p0, &p1);
                    break;
                case 2:
                case 13:
                    dx = (lev - z[i+n]) / (z[i+n+1] - z[i+n]);
                    p0[0] = x[in] + dx * (x[in+1] - x[in]);
                    p0[1] = y[im + 1];
                    
                    dy = (lev - z[i+1]) / (z[i+n+1] - z[i+1]);
                    p1[0] = x[in + 1];
                    p1[1] = y[im] + dy * (y[im+1] - y[im]);
                    add_linseg(ls, &lsnum, &lscnt, &p0, &p1);
                    break;
                case 3:
                case 12:
                    dx = (lev - z[i]) / (z[i+1] - z[i]);
                    p0[0] = x[in] + dx * (x[in+1] - x[in]);
                    p0[1] = y[im];
                    
                    dx = (lev - z[i+n]) / (z[i+n+1] - z[i+n]);
                    p1[0] = x[in] + dx * (x[in+1] - x[in]);
                    p1[1] = y[im+1];
                    add_linseg(ls, &lsnum, &lscnt, &p0, &p1);
                    break;
                case 4:
                case 11:
                    dy = (lev - z[i]) / (z[i+n] - z[i]);
                    p0[0] = x[in];
                    p0[1] = y[im] + dy * (y[im+1] - y[im]);
                    
                    dx = (lev - z[i+n]) / (z[i+n+1] - z[i+n]);
                    p1[0] = x[in] + dx * (x[in+1] - x[in]);
                    p1[1] = y[im + 1];
                    add_linseg(ls, &lsnum, &lscnt, &p0, &p1);
                    break;
                case 5:
                    dy = (lev - z[i]) / (z[i+n] - z[i]);
                    p0[0] = x[in];
                    p0[1] = y[im] + dy * (y[im+1] - y[im]);
                    dx = (lev - z[i]) / (z[i+1] - z[i]);
                    p1[0] = x[in] + dx * (x[in+1] - x[in]);
                    p1[1] = y[im];
                    add_linseg(ls, &lsnum, &lscnt, &p0, &p1);
                    
                    dx = (lev - z[i+n]) / (z[i+n+1] - z[i+n]);
                    p0[0] = x[in]+dx * (x[in+1] - x[in]);
                    p0[1] = y[im+1];
                    dy = (lev - z[i+1]) / (z[i+n+1] - z[i+1]);
                    p1[0] = x[in+1];
                    p1[1] = y[im] + dy * (y[im+1] - y[im]);
                    add_linseg(ls, &lsnum, &lscnt, &p0, &p1);
                    break;
                case 6:
                case 9:
                    dy = (lev - z[i]) / (z[i+n] - z[i]);
                    p0[0] = x[in];
                    p0[1] = y[im] + dy * (y[im+1] - y[im]);
                    dy = (lev - z[i+1]) / (z[i+n+1] - z[i+1]);
                    p1[0] = x[in+1];
                    p1[1] = y[im]+dy * (y[im+1] - y[im]);
                    add_linseg(ls, &lsnum, &lscnt, &p0, &p1);
                    break;
                
                case 7:
                case 8:
                    dy = (lev - z[i]) / (z[i+n] - z[i]);
                    p0[0] = x[in];
                    p0[1] = y[im] + dy * (y[im+1] - y[im]);
                    dx = (lev - z[i]) / (z[i+1] - z[i]);
                    p1[0] = x[in] + dx * (x[in+1] - x[in]);
                    p1[1] = y[im];
                    add_linseg(ls, &lsnum, &lscnt, &p0, &p1);
                    break;
                case 10:
                    dy = (lev - z[i]) / (z[i+n] - z[i]);
                    p0[0] = x[in];
                    p0[1] = y[im] + dy * (y[im+1] - y[im]);
                    dx = (lev - z[i+n]) / (z[i+n+1] - z[i+1]);
                    p1[0] = x[in] + dx * (x[in+1] - x[in]);
                    p1[1] = y[im+1];
                    add_linseg(ls, &lsnum, &lscnt, &p0, &p1);
                    
                    dx = (lev - z[i]) / (z[i+1] - z[i]);
                    p0[0] = x[in]+dx * (x[in+1] - x[in]);
                    p0[1] = y[im];
                    dy = (lev - z[i+1]) / (z[i+n+1] - z[i+1]);
                    p1[0] = x[in+1];
                    p1[1] = y[im] + dy * (y[im+1] - y[im]);
                    add_linseg(ls, &lsnum, &lscnt, &p0, &p1);
                    break;
                default:
                    break;
            }
        }
    }
    *ls = (linseg_t*)realloc(*ls, lscnt * sizeof(linseg_t));
    *sz = lscnt;
    return RES_OK;
}

/*******************************************************************************
Расчет дистанции между точками
*******************************************************************************/
double dist(point2d_t* p0, point2d_t* p1)
{
    double dx, dy;
    dx = p0[0][0] - p1[0][0];
    dy = p0[0][1] - p1[0][1];
    return sqrt(dx*dx + dy*dy);
}

/*******************************************************************************
Стыковка концов отрезков в линию
*******************************************************************************/
int line_create(linseg_t* ls, int nls, line2d_t* line)
{
    int i, j, c, n, k;
    if(!line || !ls)
        return ERROR_PTR;
    // printf("*line = %x, *np = %d \n", *line, *np);
    i = 0;
    while(!(ls[i].flag) && i < nls)
        i++;
    if(i==nls)
    {
        line->npoints = 0;
        return RES_OK;
    }
    //printf("i = %d, ls[%d] = [%.1f %.1f] -- [%.1f, %.1f]\n", 
    //        i, i,   ls[i].p[0][0], ls[i].p[0][1], 
    //                ls[i].p[1][0], ls[i].p[1][1]);
    n = BSIZE;
    if(line->points == NULL)
        line->points = (point2d_t*)malloc(n*sizeof(point2d_t));
    else
        line->points = (point2d_t*)realloc(line->points, n*sizeof(point2d_t));

    c = 0;
    line->points[c][0] = ls[i].p[0][0];
    line->points[c][1] = ls[i].p[0][1];
    c++;
    line->points[c][0] = ls[i].p[1][0];
    line->points[c][1] = ls[i].p[1][1];
    c++;
    ls[i].flag = 0;
    
    for(i = 0; i < nls; i++)
    {
        for(j = 0; j < nls; j++)
        {
            if(ls[j].flag)
            {
                //сравниваем последнюю точку линии с первой точкой отрезка ls[j]
                if(dist(line->points+c-1, ls[j].p) < DIST_TH)
                {
        
                    // проверяем выделение памяти.
                    // при необходимости увеличиваем
                    if(c>=n)
                    {
                        n += BSIZE;
                        line->points = (point2d_t*)realloc(line->points, n*sizeof(point2d_t));
                    }
                    //если первая точка отрезка совпадает, то добавляем
                    //в линию вторую точку СЗАДИ
                    line->points[c][0] = ls[j].p[1][0];
                    line->points[c][1] = ls[j].p[1][1];
                    ls[j].flag = 0;
                    c++;
                }
            }
            if(ls[j].flag)
            {
                // сравниваем последнюю точку линии 
                // со второй точкой отрезка ls[j]
                if(dist(line->points+c-1, ls[j].p+1) < DIST_TH)
                {
                    // проверяем выделение памяти.
                    // при необходимости увеличиваем
                    if(c>=n)
                    {
                        n += BSIZE;
                        line->points = (point2d_t*)realloc(line->points, n*sizeof(point2d_t));
                    }
                    //если вторая точка совпадает, то добавляем
                    //в линию первую точку СЗАДИ
                    line->points[c][0] = ls[j].p[0][0];
                    line->points[c][1] = ls[j].p[0][1];
                    ls[j].flag = 0;
                    c++;
                }
            }
            
            if(ls[j].flag)
            {
                //сравниваем последнюю точку линии с первой точкой отрезка ls[j]
                if(dist(line->points, ls[j].p) < DIST_TH)
                {
        
                    // проверяем выделение памяти.
                    // при необходимости увеличиваем
                    if(c>=n)
                    {
                        n += BSIZE;
                        line->points = (point2d_t*)realloc(line->points, n*sizeof(point2d_t));
                    }
                    //если первая точка отрезка совпадает, то добавляем
                    //в линию вторую точку СПЕРЕДИ
                    for(k = c; k>0; k--)
                    {
                        line->points[k][0] = line->points[k-1][0];
                        line->points[k][1] = line->points[k-1][1]; 
                    }
                    line->points[0][0] = ls[j].p[1][0];
                    line->points[0][1] = ls[j].p[1][1];
                    ls[j].flag = 0;
                    c++;
                }
            }
            if(ls[j].flag)
            {
                // сравниваем последнюю точку линии 
                // со второй точкой отрезка ls[j]
                if(dist(line->points, ls[j].p+1) < DIST_TH)
                {
                    // проверяем выделение памяти.
                    // при необходимости увеличиваем
                    if(c>=n)
                    {
                        n += BSIZE;
                        line->points = (point2d_t*)realloc(line->points, n*sizeof(point2d_t));
                    }
                    //если вторая точка совпадает, то добавляем
                    //в линию первую точку СПЕРЕДИ
                    
                    for(k = c; k>0; k--)
                    {
                        line->points[k][0] = line->points[k-1][0];
                        line->points[k][1] = line->points[k-1][1]; 
                    }
                    
                    line->points[0][0] = ls[j].p[0][0];
                    line->points[0][1] = ls[j].p[0][1];
                    ls[j].flag = 0;
                    c++;
                }
            }
            
        }
    }
    line->points = (point2d_t*)realloc(line->points, c*sizeof(point2d_t));
    line->npoints = c;
    return RES_OK;
}

