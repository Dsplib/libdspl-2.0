/*
* Copyright (c) 2015-2019 Sergey Bakhurin
* Digital Signal Processing Library [http://dsplib.org]
*
* This file is part of libdspl-2.0.
*
* is free software: you can redistribute it and/or modify
* it under the terms of the GNU Lesser  General Public License as published by
* the Free Software Foundation, either version 3 of the License, or
* (at your option) any later version.
*
* DSPL is distributed in the hope that it will be useful,
* but WITHOUT ANY WARRANTY; without even the implied warranty of
* MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
* GNU General Public License for more details.
*
* You should have received a copy of the GNU Lesser General Public License
* along with Foobar.  If not, see <http://www.gnu.org/licenses/>.
*/
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "dspl.h"
#include "dspl_internal.h"






/*******************************************************************************
matrix_create
*******************************************************************************/
int DSPL_API matrix_create(matrix_t* a, int n, int m, int type)
{
  if(!a)
    return ERROR_PTR;
  if(n < 1 || m < 1)
    return ERROR_MATRIX_SIZE;

  if(a->dat)
  {
    a->dat = (type & DAT_MASK) ?
      (void*) realloc(a->dat, n*m*sizeof(complex_t)):
      (void*) realloc(a->dat, n*m*sizeof(double));
  }
  else
  {
    a->dat = (type & DAT_MASK) ?
      (void*) malloc(n*m*sizeof(complex_t)):
      (void*) malloc(n*m*sizeof(double));
  }

  a->n = n;
  a->m = m;
  a->type = type;

  return RES_OK;
}


/*******************************************************************************
matrix_create eye
*******************************************************************************/
int DSPL_API matrix_create_eye(matrix_t* a, int n, int type)
{
  double    *pr;
  complex_t *pc;
  int err, m, k;
  
  err = matrix_create(a, n, n, type);
  if(err != RES_OK)
    return RES_OK;
  
  k = 0;
  if((a->type & DAT_MASK) == DAT_DOUBLE)
  {
    pr = (double*) a->dat;
    memset(pr, 0, n*n*sizeof(double));
    for(m = 0; m < n; m++)
    {
      pr[k] = 1.0;
      k += n+1;
    }
  }
  
  if((a->type & DAT_MASK) == DAT_COMPLEX)
  {
    pc = (complex_t*) a->dat;
    memset(pc, 0, n*n*sizeof(complex_t));
    for(m = 0; m < n; m++)
    {
      RE(pc[k]) = 1.0;
      k += n+1;      
    }
  }
  return RES_OK;
}



/*******************************************************************************
matrix_free
*******************************************************************************/
void DSPL_API matrix_free(matrix_t* a)
{
  if(!a)
    return;
  if(a->dat)
    free(a->dat);
  a->n = a->m = a->type = 0;
}


/*******************************************************************************
matrix LU decomposition
*******************************************************************************/
int DSPL_API matrix_lu(matrix_t* a, matrix_t* L, matrix_t* U, matrix_t* P)
{
  int err, k, n, m, N, ind;
  double *rl, *ru, mu, ukk, gmax;
   
  if(!a || !L || !U || !P)
    return ERROR_PTR;
  
  if(a->n != a->m || a->n < 1)
    return ERROR_MATRIX_SIZE;
  
  N = a->n;
  err = matrix_create(L, N, N, a->type);
  if(err != RES_OK)
    return err;
  
  err = matrix_create(U, N, N, a->type);
  if(err != RES_OK)
    return err;
  err = matrix_create_eye(P, N, a->type);
  if(err != RES_OK)
    return err;

  if((a->type & DAT_MASK) == DAT_DOUBLE)
  {
    rl = (double*)L->dat;
    ru = (double*)U->dat;
    
    memcpy(ru, (double*)a->dat, N*N*sizeof(double));
    memset(rl, 0,  N*N*sizeof(double));
    
    find_max_abs(ru, N*N, &gmax, NULL);
    for(k = 0; k < N; k++)
    {
      find_max_abs(ru+k*N+k, N-k, NULL, &ind);
      ind += k;
      matrix_swap_rows(U, k, ind);
      matrix_swap_rows(L, k, ind);
      matrix_swap_rows(P, k, ind);
      ukk = ru[N*k+k];
      if(fabs(ukk / gmax) < MATRIX_SINGULAR_THRESHOLD)
        return ERROR_MATRIX_SINGULAR;
      
      for(m = k+1; m < N; m++)
      {
        mu = ru[m+k*N] / ukk;
        rl[m+k*N] = mu;
        for(n = k; n < N; n++)
        {
          ru[m + n*N] -= ru[k + n*N] * mu;
        }
      }
    }
    for(n =0; n < N; n++)
      rl[n+n*N] = 1.0;
  }
  
  return RES_OK;
}



/*******************************************************************************
matrix transposition
*******************************************************************************/
int DSPL_API matrix_print(matrix_t* a, const char* name, const char* format)
{
  int n,m;
  if(!a)
    return ERROR_PTR;
  if(!a->dat)
    return ERROR_PTR;

  if((a->type & DAT_MASK) == DAT_DOUBLE)
  {
    printf("\nMatrix %s size [%d x %d] type: real\n",
              name, a->n, a->m);
    double* p = (double*)(a->dat);
    for(n = 0; n < a->n; n++)
    {
      for(m = 0; m < a->m; m++)
      {
        printf(format, p[m*a->n + n]);
      }
      printf("\n");
    }
  }


  if((a->type & DAT_MASK) == DAT_COMPLEX)
  {
    printf("\nMatrix %s size [%d x %d] type: complex\n",
              name, a->n, a->m);
    complex_t* p = (complex_t*)(a->dat);
    for(n = 0; n < a->n; n++)
    {
      for(m = 0; m < a->m; m++)
      {
        printf(format, RE(p[m*a->n + n]), IM(p[m*a->n + n]));
      }
      printf("\n");
    }
  }
  return RES_OK;
}


/*******************************************************************************
matrix swap 2 elements
*******************************************************************************/
int  DSPL_API matrix_swap(matrix_t* a, int r0, int c0, int r1, int c1)
{
  double    tr;
  complex_t tc;
  double *pr;
  complex_t *pc;

  if(!a)
    return ERROR_PTR;
  if(r0 >= a->n || r1 >= a->n || c0 >= a->m || c1 >= a->m)
    return ERROR_MATRIX_INDEX;
  
  if((a->type & DAT_MASK) == DAT_DOUBLE)
  { 
    pr = (double*)(a->dat);
    tr = pr[r0 + c0 * a->n];
    pr[r0 + c0 * a->n] = pr[r1 + c1 * a->n];
    pr[r1 + c1 * a->n] = tr;
  }
  
  if((a->type & DAT_MASK) == DAT_COMPLEX)
  {
    pc = (complex_t*)(a->dat);
    RE(tc) = RE(pc[r0 + c0 * a->n]);
    IM(tc) = IM(pc[r0 + c0 * a->n]);
    
    RE(pc[r0 + c0 * a->n]) = RE(pc[r1 + c1 * a->n]);
    IM(pc[r0 + c0 * a->n]) = IM(pc[r1 + c1 * a->n]);
    
    RE(pc[r1 + c1 * a->n]) = RE(tc);
    IM(pc[r1 + c1 * a->n]) = IM(tc);
  }
  return RES_OK;
}



/*******************************************************************************
matrix swap 2 rows
*******************************************************************************/
int  DSPL_API matrix_swap_rows(matrix_t* a, int r0, int r1)
{
  int c, err;
  if(!a)
    return ERROR_PTR;
  if(r0 >= a->n || r1 >= a->n)
    return ERROR_MATRIX_INDEX;
  
  for(c = 0; c < a->m; c++)
  {
    err = matrix_swap(a, r0, c, r1, c);
    if(err != RES_OK)
      break;
  }
  return err;
}




/*******************************************************************************
matrix transposition
*******************************************************************************/
int DSPL_API matrix_transpose(matrix_t* a, matrix_t* b)
{
  int err;
  if(!a || !b)
    return ERROR_PTR;

  err = matrix_create(b, a->m, a->n, a->type);
  if(err != RES_OK)
    return err;

  if((a->type & DAT_MASK) == DAT_DOUBLE)
    transpose((double*)(a->dat), a->n, a->m, (double*)(b->dat));
  if((a->type & DAT_MASK) == DAT_COMPLEX)
    transpose_cmplx((complex_t*)(a->dat), a->n, a->m, (complex_t*)(b->dat));

  return RES_OK;
}





/*******************************************************************************
matrix Hermite transposition
*******************************************************************************/
int DSPL_API matrix_transpose_hermite(matrix_t* a, matrix_t* b)
{
  int err;
  if(!a || !b)
    return ERROR_PTR;

  err = matrix_create(b, a->m, a->n, a->type);
  if(err != RES_OK)
    return err;

  if((a->type & DAT_MASK) == DAT_DOUBLE)
    transpose((double*)(a->dat), a->n, a->m, (double*)(b->dat));
  if((a->type & DAT_MASK) == DAT_COMPLEX)
    transpose_hermite((complex_t*)(a->dat), a->n, a->m, (complex_t*)(b->dat));

  return RES_OK;
}







/*******************************************************************************
Real matrx transpose
*******************************************************************************/
void transpose(double* a, int n, int m, double* b)
{
  int p, q, i, j, aind, bind;

  for(p = 0; p < n - DSPL_MATRIX_BLOCK; p+=DSPL_MATRIX_BLOCK)
  {
    for(q = 0; q < m - DSPL_MATRIX_BLOCK; q+=DSPL_MATRIX_BLOCK)
    {
      for(i = 0; i < DSPL_MATRIX_BLOCK; i++)
      {
        for(j = 0; j < DSPL_MATRIX_BLOCK; j++)
        {
          aind = (q+j) * n + p + i;
          bind = (p+i) * m + q + j;
          b[bind] = a[aind];
        }
      }
    }
  }
  for(i = p; i < n; i++)
    for(j = 0; j < m; j++)
      b[i*m + j] = a[j*n+i];

  for(i = 0; i < p; i++)
    for(j = q; j < m; j++)
      b[i*m + j] = a[j*n+i];

}





/*******************************************************************************
Complex matrx transpose
*******************************************************************************/
void transpose_cmplx(complex_t* a, int n, int m, complex_t* b)
{
  int p, q, i, j, aind, bind;

  for(p = 0; p < n - DSPL_MATRIX_BLOCK; p+=DSPL_MATRIX_BLOCK)
  {
    for(q = 0; q < m - DSPL_MATRIX_BLOCK; q+=DSPL_MATRIX_BLOCK)
    {
      for(i = 0; i < DSPL_MATRIX_BLOCK; i++)
      {
        for(j = 0; j < DSPL_MATRIX_BLOCK; j++)
        {
          aind = (q+j) * n + p + i;
          bind = (p+i) * m + q + j;
          RE(b[bind]) = RE(a[aind]);
          IM(b[bind]) = IM(a[aind]);
        }
      }
    }
  }
  for(i = p; i < n; i++)
  {
    for(j = 0; j < m; j++)
    {
      RE(b[i*m + j]) = RE(a[j*n+i]);
      IM(b[i*m + j]) = IM(a[j*n+i]);
    }
  }

  for(i = 0; i < p; i++)
  {
    for(j = q; j < m; j++)
    {
      RE(b[i*m + j]) = RE(a[j*n+i]);
      IM(b[i*m + j]) = IM(a[j*n+i]);
    }
  }
}





/*******************************************************************************
Hermite matrx transpose
*******************************************************************************/
void transpose_hermite(complex_t* a, int n, int m, complex_t* b)
{
  int p, q, i, j, aind, bind;

  for(p = 0; p < n - DSPL_MATRIX_BLOCK; p+=DSPL_MATRIX_BLOCK)
  {
    for(q = 0; q < m - DSPL_MATRIX_BLOCK; q+=DSPL_MATRIX_BLOCK)
    {
      for(i = 0; i < DSPL_MATRIX_BLOCK; i++)
      {
        for(j = 0; j < DSPL_MATRIX_BLOCK; j++)
        {
          aind = (q+j) * n + p + i;
          bind = (p+i) * m + q + j;
          RE(b[bind]) =  RE(a[aind]);
          IM(b[bind]) = -IM(a[aind]);
        }
      }
    }
  }
  for(i = p; i < n; i++)
  {
    for(j = 0; j < m; j++)
    {
      RE(b[i*m + j]) = RE(a[j*n+i]);
      IM(b[i*m + j]) = -IM(a[j*n+i]);
    }
  }

  for(i = 0; i < p; i++)
  {
    for(j = q; j < m; j++)
    {
      RE(b[i*m + j]) =  RE(a[j*n+i]);
      IM(b[i*m + j]) = -IM(a[j*n+i]);
    }
  }
}

