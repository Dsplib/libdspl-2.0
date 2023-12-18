#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "dspl.h"

#define N  3
#define M  4



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
	if((n == 0) && ((*ls)==NULL))
	{
		n = 128;
		(*ls) = (linseg_t*)malloc(n * sizeof(linseg_t));
	}
	else
	{
		if(c >= n)
		{
			n += 128;
			(*ls) = (linseg_t*)realloc((*ls), n * sizeof(linseg_t));
		}
	}
	(*ls)[c].p[0][0] = p0[0][0];
	(*ls)[c].p[0][1] = p0[0][1];
	(*ls)[c].p[1][0] = p1[0][0];
	(*ls)[c].p[1][1] = p1[0][1];
	(*ls)[c].flag = 0;
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
					break;
				case 1:
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
					dx = (lev - z[i+n]) / (z[i+n+1] - z[i+n]);
					p0[0] = x[in] + dx;
					p0[1] = y[im + 1];
					
					dy = (lev - z[i+1]) / (z[i+n+1] - z[i+1]);
					p1[0] = x[in + 1];
					p1[1] = y[im] + dy;
					add_linseg(ls, &lsnum, &lscnt, &p0, &p1);
					break;
				case 3:
					dx = (lev - z[i]) / (z[i+1] - z[i]);
					p0[0] = (double) in + dx;
					p0[1] = y[im];
					p1[0] = (double) in + dx;
					p1[1] = y[im + 1];
					add_linseg(ls, &lsnum, &lscnt, &p0, &p1);
					break;
				case 4:
					dy = (lev - z[i]) / (z[i+n] - z[i]);
					p0[0] = x[in];
					p0[1] = y[im] + dy;
					dx = (lev - z[i+n]) / (z[i+n+1] - z[i+n]);
					p1[0] = x[in] + dx;
					p1[1] = y[im + 1];
					add_linseg(ls, &lsnum, &lscnt, &p0, &p1);
					break;
				case 5:
					break;
				case 6:
					break;
				case 7:
					break;
				case 8:
					dy = (lev - z[i]) / (z[i+n] - z[i]);
					p0[0] = x[in];
					p0[1] = y[im] + dy;
					dx = (lev - z[i]) / (z[i+1] - z[i]);
					p1[0] = x[in] + dx;
					p1[1] = y[im];
					add_linseg(ls, &lsnum, &lscnt, &p0, &p1);
					break;
				case 9:
					break;
				case 10:
					break;
				case 12:
					dx = (lev - z[i]) / (z[i+1] - z[i]);
					p0[0] = x[in] + dx;
					p0[1] = y[im];
					
					p1[0] = x[in] + dx;
					p1[1] = y[im+1];
					add_linseg(ls, &lsnum, &lscnt, &p0, &p1);
					break;
				case 13:
					break;
				case 14:
					break;
				case 15:
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

int main(int argc, char* argv[])
{
    /* Matrix 
    z = [ 0   0   0   0;
          0   1   1   0;
          0   0   0   0];
    in array a by columns
    */
    double z[N*M] = { 0, 0, 0,   0, 1, 0,   0, 1, 0,   0, 0, 0};
	double x[N] = {0.0, 1.0, 2.0};
    double y[M] = {0.0, 1.0, 2.0, 3.0};
	linseg_t *ls = NULL;
	int nls, i;
	
	linseg_create(z, x, y, N, M, 0.1, &ls, &nls);
	for(i =0; i< nls; i++)
		printf("%d, [%.1f %.1f] -- [%.1f %.1f]\n", i, ls[i].p[0][0], ls[i].p[0][1], ls[i].p[1][0], ls[i].p[1][1]);
    return 0;
}
