#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "dspl.h"

#define N  3
#define M  4
#define BSIZE 256

typedef double    point2d_t[2];

typedef struct
{
    point2d_t p[2];
    int flag;
} linseg_t;
 
 
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
            
            printf("%d, %d, %d\n", in, im, t);
            switch(t)
            {
                case 0:
                case 15:
                    break;
                case 1:
                case 14:
                    // z[i,j] * (1-dx) + z[i+1, j] * dx   = lev
                    // z[i,j]  + dx *(z[i+1, j] - z[i,j]) = lev
                    dx = (lev - z[i]) / (z[i+1] - z[i]);
                    p0[0] = x[in] + dx;
                    p0[1] = y[im];
                    
                    dy = (lev - z[i+1]) / (z[i+n+1] - z[i+1]);
                    p1[0] = x[in+1];
                    p1[1] = y[im] + dy;
                    add_linseg(ls, &lsnum, &lscnt, &p0, &p1);
                    break;
                case 2:
                case 13:
                    dx = (lev - z[i+n]) / (z[i+n+1] - z[i+n]);
                    p0[0] = x[in] + dx;
                    p0[1] = y[im + 1];
                    
                    dy = (lev - z[i+1]) / (z[i+n+1] - z[i+1]);
                    p1[0] = x[in + 1];
                    p1[1] = y[im] + dy;
                    add_linseg(ls, &lsnum, &lscnt, &p0, &p1);
                    break;
                case 3:
                case 12:
                    dx = (lev - z[i]) / (z[i+1] - z[i]);
                    p0[0] = x[in] + dx;
                    p0[1] = y[im];
                    
                    dx = (lev - z[i+n]) / (z[i+n+1] - z[i+n]);
                    p1[0] = x[in] + dx;
                    p1[1] = y[im+1];
                    add_linseg(ls, &lsnum, &lscnt, &p0, &p1);
                    break;
                case 4:
                case 11:
                    dy = (lev - z[i]) / (z[i+n] - z[i]);
                    p0[0] = x[in];
                    p0[1] = y[im] + dy;
                    
                    dx = (lev - z[i+n]) / (z[i+n+1] - z[i+n]);
                    p1[0] = x[in] + dx;
                    p1[1] = y[im + 1];
                    add_linseg(ls, &lsnum, &lscnt, &p0, &p1);
                    break;
                case 5:
                    dy = (lev - z[i]) / (z[i+n] - z[i]);
                    p0[0] = x[in];
                    p0[1] = y[im] + dy;
                    dx = (lev - z[i]) / (z[i+1] - z[i]);
                    p1[0] = x[in] + dx;
                    p1[1] = y[im];
                    add_linseg(ls, &lsnum, &lscnt, &p0, &p1);
                    
                    dx = (lev - z[i+n]) / (z[i+n+1] - z[i+n]);
                    p0[0] = x[in]+dx;
                    p0[1] = y[im+1];
                    dy = (lev - z[i+1]) / (z[i+n+1] - z[i+1]);
                    p1[0] = x[in+1];
                    p1[1] = y[im] + dy;
                    add_linseg(ls, &lsnum, &lscnt, &p0, &p1);
                    break;
                case 6:
                case 9:
                    dy = (lev - z[i]) / (z[i+n] - z[i]);
                    p0[0] = x[in];
                    p0[1] = y[im] + dy;
                    dy = (lev - z[i+1]) / (z[i+n+1] - z[i+1]);
                    p1[0] = x[in+1];
                    p1[1] = y[im];
                    add_linseg(ls, &lsnum, &lscnt, &p0, &p1);
                    break;
                
                case 7:
                case 8:
                    dy = (lev - z[i]) / (z[i+n] - z[i]);
                    p0[0] = x[in];
                    p0[1] = y[im] + dy;
                    dx = (lev - z[i]) / (z[i+1] - z[i]);
                    p1[0] = x[in] + dx;
                    p1[1] = y[im];
                    add_linseg(ls, &lsnum, &lscnt, &p0, &p1);
                    break;
                case 10:
                    dy = (lev - z[i]) / (z[i+n] - z[i]);
                    p0[0] = x[in];
                    p0[1] = y[im] + dy;
                    dx = (lev - z[i+n]) / (z[i+n+1] - z[i+1]);
                    p1[0] = x[in] + dx;
                    p1[1] = y[im+1];
                    add_linseg(ls, &lsnum, &lscnt, &p0, &p1);
                    
                    dx = (lev - z[i]) / (z[i+1] - z[i]);
                    p0[0] = x[in]+dx;
                    p0[1] = y[im];
                    dy = (lev - z[i+1]) / (z[i+n+1] - z[i+1]);
                    p1[0] = x[in+1];
                    p1[1] = y[im] + dy;
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


double dist(point2d_t* p0, point2d_t* p1)
{
    double dx, dy;
    dx = p0[0][0] - p1[0][0];
    dy = p0[0][1] - p1[0][1];
    return sqrt(dx*dx + dy*dy);
}



int line_create(linseg_t* ls, int nls, point2d_t** line, int* np)
{
    int i, j, c, n;
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
                //сравниваем с первой точкой отрезка ls[j]
                if(dist((*line)+c-1, ls[j].p) < 1E-8)
                {
                    // printf("c0 = %d, j0 = %d\n", c, j);
                    
                    // проверяем выделение памяти.
                    // при необходимости увеличиваем
                    if(c>=n)
                    {
                        n += BSIZE;
                        (*line) = (point2d_t*)realloc(*line, n*sizeof(point2d_t));
                    }
                    //если первая точка совпадает, то добавляем
                    //в линию вторую точку
                    (*line)[c][0] = ls[j].p[1][0];
                    (*line)[c][1] = ls[j].p[1][1];
                    ls[j].flag = 0;
                    c++;
                }
            }
            if(ls[j].flag)
            {
                //сравниваем со второй точкой отрезка ls[j]
                if(dist((*line)+c-1, ls[j].p+1) < 1E-8)
                {
                    // printf("c1 = %d, j1 = %d\n", c, j);
                    
                    // проверяем выделение памяти.
                    // при необходимости увеличиваем
                    if(c>=n)
                    {
                        n += BSIZE;
                        (*line) = (point2d_t*)realloc(*line, n*sizeof(point2d_t));
                    }
                    //если вторая точка совпадает, то добавляем
                    //в линию первую точку
                    (*line)[c][0] = ls[j].p[0][0];
                    (*line)[c][1] = ls[j].p[0][1];
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





int main(int argc, char* argv[])
{
    /* Matrix 
    z = [ 0   0   1   0;
          0   1   1   0;
          0   0   0   0];
    in array a by columns
    */
    double z[N*M] = { 0, 0, 0,   0, 1, 0,   1, 1, 0,   0, 0, 0};
    double x[N] = {0.0, 1.0, 2.0};
    double y[M] = {0.0, 1.0, 2.0, 3.0};
    linseg_t *ls = NULL;
    int nls, i, j;
    point2d_t** line = NULL;
    int *np = NULL;
    
    line = (point2d_t**)malloc(10*sizeof(point2d_t*));
    memset(line,0,10*sizeof(point2d_t*));
    
    np = (int*)malloc(10*sizeof(int));
    memset(np, 0, 10*sizeof(int));

    linseg_create(z, x, y, N, M, 0.5, &ls, &nls);
    for(i =0; i< nls; i++)
        printf("%d, [%.1f %.1f] -- [%.1f %.1f]\n", i, ls[i].p[0][0], ls[i].p[0][1], ls[i].p[1][0], ls[i].p[1][1]);


    i = -1;
    do
    {
        i++;
        line_create(ls, nls, line+i, np);
        for(j =0; j< np[i]; j++)
            printf("%[%.1f %.1f] -- ", line[i][j][0], line[i][j][1]);
        printf("\n");
    }while(np[i]);


    return 0;
}
