#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "dspl.h"

#define N  3
#define M  4


typedef double    point2d_t[2];
typedef point2d_t linseg_t[2];
 
 
int add_linseg(linseg_t** ls, int* lsnum, int* lscnt, point2d_t* p0, point2d_t* p1)
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
	(*ls)[c][0][0] = p0[0][0];
	(*ls)[c][0][1] = p0[0][1];
	(*ls)[c][1][0] = p1[0][0];
	(*ls)[c][1][1] = p1[0][1];
	c++;
	(*lsnum) = n;
	(*lscnt) = c;
	
	return RES_OK;
}

int linseg_create(int* z, int n, int m, linseg_t** ls, int* sz)
{
    int lsnum, lscnt, t, in, im, i;
	point2d_t p0 = {0};
	point2d_t p1 = {0};
	
	if((ls== NULL)||(z==NULL))
		return ERROR_PTR;
	lsnum = 0;
	lscnt = 0;
	
	for(in = 0; in < n-1; in++)
	{
		for(im = 0; im < m-1; im++)
		{
			i = in + im * n;
			t = z[i]*8 + z[i+n]*4 + z[i+n+1]*2 + z[i+1];
			printf("%d, %d, %d\n", in, im, t);
			switch(t)
			{
				case 0:
					break;
				case 1:
					p0[0] = (double) in + 0.5;
					p0[1] = (double) im;
					p1[0] = (double) in + 1.0;
					p1[1] = (double) im + 0.5;
					add_linseg(ls, &lsnum, &lscnt, &p0, &p1);
					break;
				case 2:
					p0[0] = (double) in + 0.5;
					p0[1] = (double) im + 1.0;
					p1[0] = (double) in + 1.0;
					p1[1] = (double) im + 0.5;
					add_linseg(ls, &lsnum, &lscnt, &p0, &p1);
					break;
				case 3:
					p0[0] = (double) in + 0.5;
					p0[1] = (double) im;
					p1[0] = (double) in + 0.5;
					p1[1] = (double) im + 1.0;
					add_linseg(ls, &lsnum, &lscnt, &p0, &p1);
					break;
				case 4:
					p0[0] = (double) in;
					p0[1] = (double) im + 0.5;
					p1[0] = (double) in + 0.5;
					p1[1] = (double) im + 1.0;
					add_linseg(ls, &lsnum, &lscnt, &p0, &p1);
					break;
				case 5:
					break;
				case 6:
					break;
				case 7:
					break;
				case 8:
					p0[0] = (double) in;
					p0[1] = (double) im + 0.5;
					p1[0] = (double) in + 0.5;
					p1[1] = (double) im;
					add_linseg(ls, &lsnum, &lscnt, &p0, &p1);
					break;
				case 9:
					break;
				case 10:
					break;
				case 12:
					p0[0] = (double) in + 0.5;
					p0[1] = (double) im;
					p1[0] = (double) in + 0.5;
					p1[1] = (double) im + 1.0;
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
    int z[N*M] = { 0, 0, 0,   0, 1, 0,   0, 1, 0,   0, 0, 0};
    linseg_t *ls = NULL;
	int nls, i;
	
	linseg_create(z, N, M, &ls, &nls);
	for(i =0; i< nls; i++)
		printf("%d, [%.1f %.1f] -- [%.1f %.1f]\n", i, ls[i][0][0], ls[i][0][1], ls[i][1][0], ls[i][1][1]);
    return 0;
}
