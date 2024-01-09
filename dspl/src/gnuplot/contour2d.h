#ifndef CONTOUR2D_H
#define CONTOUR2D_H


#define BSIZE 256
#define DIST_TH 1E-8

typedef struct
{
    point2d_t p[2];
    int flag;
} linseg_t;


int add_linseg(linseg_t** ls, int* lsnum, int* lscnt, 
                point2d_t* p0, point2d_t* p1);

int linseg_create(double* z, double* x, double* y, 
                  int n, int m, double lev, 
                  linseg_t** ls, int* sz);

double dist(point2d_t* p0, point2d_t* p1);

int line_create(linseg_t* ls, int nls, line2d_t* line);  

#endif