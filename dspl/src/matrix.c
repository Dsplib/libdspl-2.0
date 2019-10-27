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
#include "blas.h"





int DSPL_API matrix_eig_cmplx(complex_t* a, int n, complex_t* v, int* info)
{
  int err;
  int sdim = 0;
  int ldvs = 1;
  int lwork = 2*n;
  if(!a || !v)
    return ERROR_PTR;
  
  if(n<1)
    return ERROR_MATRIX_SIZE;
  
  complex_t *work=(complex_t*)malloc(lwork*sizeof(complex_t)); 
  double *rwork = (double*)malloc(n*sizeof(double)); 

  zgees_("N", "N", NULL, &n, a, &n, &sdim, v, NULL, &ldvs, work, &lwork, 
         rwork, NULL, &err);
  
  if(err!=0)
  {
    if(info)
      *info = err;
    err = ERROR_LAPACK;
  }
  else
    err = RES_OK;
  
  free(work);
  free(rwork);
  return err;
}




/*******************************************************************************
Real matrix eye
*******************************************************************************/
int DSPL_API matrix_eye(double* a, int n, int m)
{
  int  p, k;
  if(!a)
    return ERROR_PTR;
  if (n < 1 || m < 1)
    return ERROR_MATRIX_SIZE;
    
  k = 0;
  memset(a, 0, n*m*sizeof(double));
  for(p = 0; p < m; p++)
  {
    a[k] = 1.0;
    k += n+1;
  }

  return RES_OK;
}


/*******************************************************************************
Complex matrix eye 
*******************************************************************************/
int DSPL_API matrix_eye_cmplx(complex_t* a, int n, int m)
{
  int p, k;
  if(!a)
    return ERROR_PTR;
  if (n < 1 || m < 1)
    return ERROR_MATRIX_SIZE;
    
  k = 0;
  memset(a, 0, n*m*sizeof(complex_t));
  for(p = 0; p < m; p++)
  {
    RE(a[k]) = 1.0;
    k += n+1;
  }

  return RES_OK;
}




/*******************************************************************************
matrix LU decomposition
******************************************************************************
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
*/


/*******************************************************************************
real matrix multiplication
*******************************************************************************/
int DSPL_API matrix_mul(double* a, int na, int ma, 
                        double* b, int nb, int mb,
                        double* c)
{
  
  double alpha = 1;
  double beta = 0.0;
  
  if(!a || !b || !c)
    return ERROR_PTR;
  if(na < 1 || ma < 1 || nb < 1 || mb < 1 || ma != nb)
    return ERROR_MATRIX_SIZE;

  /* BLAS DGEMM */
  dgemm_("N", "N", &na, &mb, &ma, &alpha, a, &na, b, &nb, &beta, c, &na);

  return RES_OK;
}



/*******************************************************************************
real matrix print
*******************************************************************************/
int DSPL_API matrix_print(double* a, int n, int m, 
                          const char* name, const char* format)
{
  int p,q;
  
  if(!a)
    return ERROR_PTR;
  if(n < 1 || m < 1)
    return ERROR_SIZE;
    
  printf("\n%s = [ %% size [%d x %d] type: real", name, n, m);
  
  for(p = 0; p < n; p++)
  {
    printf("\n");
    for(q = 0; q < m; q++)
    {
      printf(format, a[q*n + p]);
      if(q == m-1)
        printf(";");
      else
        printf(", ");
    }
  }
  printf("];\n");
  
  return RES_OK;
}


/*******************************************************************************
complex matrix print
*******************************************************************************/
int DSPL_API matrix_print_cmplx(complex_t* a, int n, int m, 
                                const char* name, const char* format)
{
  int p,q;

  if(!a)
    return ERROR_PTR;
  if(n < 1 || m < 1)
    return ERROR_MATRIX_SIZE;

  if(!a)
    return ERROR_PTR;
  if(n < 1 || m < 1)
    return ERROR_SIZE;
    
  printf("\n%s = [ %% size [%d x %d] type: complex", name, n, m);
  
  for(p = 0; p < n; p++)
  {
    printf("\n");
    for(q = 0; q < m; q++)
    {
      printf(format, RE(a[q*n + p]), IM(a[q*n + p]));
      if(q == m-1)
        printf(";");
      else
        printf(", ");
    }
  }
  printf("];\n");

  return RES_OK;
}






/*******************************************************************************
Real matrix transpose
*******************************************************************************/
int DSPL_API matrix_transpose(double* a, int n, int m, double* b)
{
  int p, q, i, j, aind, bind;
  if(!a || !b)
    return ERROR_PTR;
  if(n < 1 || m < 1)
    return ERROR_MATRIX_SIZE;
    
    
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
      
  return RES_OK;
}





/*******************************************************************************
Complex matrix transpose
*******************************************************************************/
int DSPL_API matrix_transpose_cmplx(complex_t* a, int n, int m, complex_t* b)
{
  int p, q, i, j, aind, bind;
  
  if(!a || !b)
    return ERROR_PTR;
  if(n < 1 || m < 1)
    return ERROR_MATRIX_SIZE;
    
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
  return RES_OK;
}





/*******************************************************************************
Hermite matrix transpose
*******************************************************************************/
int DSPL_API matrix_transpose_hermite(complex_t* a, int n, int m, complex_t* b)
{
  int p, q, i, j, aind, bind;
  
  if(!a || !b)
    return ERROR_PTR;
  if(n < 1 || m < 1)
    return ERROR_MATRIX_SIZE;
  
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
  
  return RES_OK;
}





/*******************************************************************************
 * Vector dot product
 ******************************************************************************/
int DSPL_API vector_dot(double* x, double* y, int n, double* p)
{
  int inc = 1;
  
  if(!x || !y || !p)
    return ERROR_PTR;
  if(n<1)
    return ERROR_SIZE;
    
  *p = ddot_(&n, x, &inc, y, &inc); 
  
  return RES_OK;
}

