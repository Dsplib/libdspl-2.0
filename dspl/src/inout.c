/*
* Copyright (c) 2015-2018 Sergey Bakhurin
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
#include <stdlib.h>
#include "dspl.h"



/*******************************************************************************
Print DSPL info
*******************************************************************************/
void DSPL_API dspl_info()
{
  printf("\n\n        D S P L - 2.0\n");
  printf("        version 2.18.10.05\n");
  printf("\n        Copyright (C) 2015-2018\n");
  printf("        Sergey Bakhurin        www.dsplib.org\n\n");
}






/******************************************************************************
Write a real array to the binary file "fn"
*******************************************************************************/
int DSPL_API writebin(void* x, int n, int dtype, char* fn)
{
  int k, res;
  FILE* pFile = NULL;

  if(!x)
    return ERROR_PTR;
  if(n < 1)
    return ERROR_SIZE;
  if(!fn)
    return  ERROR_FNAME;

  pFile = fopen(fn, "wb");
  if(pFile == NULL)
      return ERROR_FOPEN;


  if(fwrite(&dtype, sizeof(int), 1, pFile) != 1)
  {
    res = ERROR_FWRITE_SIZE;
    goto exit_label;
  }


  if(fwrite(&n, sizeof(int),   1, pFile) != 1)
  {
    res = ERROR_FWRITE_SIZE;
    goto exit_label;
  }

  k = 1;
  if(fwrite(&k, sizeof(int),   1, pFile) != 1)
  {
    res = ERROR_FWRITE_SIZE;
    goto exit_label;
  };

  switch(dtype)
  {
    case DAT_DOUBLE:
      if(fwrite((double*)x, sizeof(double), n, pFile) != n)
      {
        res = ERROR_FWRITE_SIZE;
        goto exit_label;
      }
      break;
    case DAT_COMPLEX:
      if(fwrite((complex_t*)x,
          sizeof(complex_t), n, pFile) != n)
      {
        res = ERROR_FWRITE_SIZE;
        goto exit_label;
      }
      break;
    default:
      res = ERROR_DAT_TYPE;
      goto exit_label;
  }
  res = RES_OK;
exit_label:
  if(pFile)
    fclose(pFile);
  return res;
}






/******************************************************************************
Write a real arrays to the text file "fn"
*******************************************************************************/
int DSPL_API writetxt(double* x, double *y, int n, char* fn)
{
  int k;
  FILE* pFile = NULL;

  if(!x)
    return ERROR_PTR;
  if(n < 1)
    return ERROR_SIZE;
  if(!fn)
    return ERROR_FNAME;

  pFile = fopen(fn, "w+");
  if(pFile == NULL)
    return ERROR_FOPEN;

  if(y)
    for(k = 0; k < n; k++)
      fprintf(pFile, "%+.12E\t%+.12E\n", x[k], y[k]);
  else
    for(k = 0; k < n; k++)
      fprintf(pFile, "%+.12E\n", x[k]);

  fclose(pFile);
  return RES_OK;
}






/******************************************************************************
 * Write a 3d plot data to file "fn" (pgfplots3d accepteble)
 ******************************************************************************/
int DSPL_API writetxt_3d(double* x, int nx, double *y, int ny,
               double* z, char* fn)
{
  int k, n;
  FILE* pFile = NULL;

  if(!x || !y || !z)
    return ERROR_PTR;
  if(nx < 1 || ny < 1)
    return ERROR_SIZE;
  if(!fn)
    return ERROR_FNAME;

  pFile = fopen(fn, "w+");
  if(pFile == NULL)
    return ERROR_FOPEN;

  for(k = 0; k < ny; k++)
  {
          for(n = 0; n < nx; n++)
      fprintf(pFile, "%+.12E\t%+.12E\t%+.12E\n",
              x[n], y[k], z[n+k*nx]);
    fprintf(pFile, "\n");

  }
  fclose(pFile);
  return RES_OK;
}







/******************************************************************************
 * Write a 3d line data to file "fn" (pgfplots3d accepteble)
 ******************************************************************************/
int DSPL_API writetxt_3dline(double* x, double *y, double* z, int n, char* fn)
{
  int k;
  FILE* pFile = NULL;

  if(!x || !y || !z)
    return ERROR_PTR;
  if(n < 1)
    return ERROR_SIZE;
  if(!fn)
    return ERROR_FNAME;

  pFile = fopen(fn, "w+");
  if(pFile == NULL)
    return ERROR_FOPEN;

  for(k = 0; k < n; k++)
    fprintf(pFile, "%+.12E\t%+.12E\t%+.12E\n", x[k], y[k], z[k]);
  fprintf(pFile, "\n");
  fclose(pFile);
  return RES_OK;
}




/******************************************************************************
Write a real part of coplex array to the text file "fn"
*******************************************************************************/
int DSPL_API writetxt_cmplx_re(double* x, complex_t *y, int n, char* fn)
{
  int k;
  FILE* pFile = NULL;

  if(!x)
    return ERROR_PTR;
  if(n < 1)
    return ERROR_SIZE;
  if(!fn)
    return ERROR_FNAME;

  pFile = fopen(fn, "w+");
  if(pFile == NULL)
    return ERROR_FOPEN;

  if(x)
    for(k = 0; k < n; k++)
      fprintf(pFile, "%+.12E\t%+.12E\n", x[k], RE(y[k]));
  else
    for(k = 0; k < n; k++)
      fprintf(pFile, "%+.12E\n", RE(y[k]));

  fclose(pFile);
  return RES_OK;
}



/******************************************************************************
Write a image part of coplex array to the text file "fn"
*******************************************************************************/
int DSPL_API writetxt_cmplx_im(double* x, complex_t *y, int n, char* fn)
{
  int k;
  FILE* pFile = NULL;

  if(!x)
    return ERROR_PTR;
  if(n < 1)
    return ERROR_SIZE;
  if(!fn)
    return ERROR_FNAME;

  pFile = fopen(fn, "w+");
  if(pFile == NULL)
    return ERROR_FOPEN;

  if(x)
    for(k = 0; k < n; k++)
      fprintf(pFile, "%+.12E\t%+.12E\n", x[k], IM(y[k]));
  else
    for(k = 0; k < n; k++)
      fprintf(pFile, "%+.12E\n", IM(y[k]));

  fclose(pFile);
  return RES_OK;
}

