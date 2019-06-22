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


#include <stdlib.h>
#include <string.h>
#include "dspl.h"




/*******************************************************************************
\ingroup FILTER_CONV_GROUP
\fn int conv(double* a, int na, double* b, int nb, double* c) 
\brief Real vectors linear convolution.

Function convolves two real vectors \f$ c = a * b\f$ length `na` and `nb`.
The output convolution is a vector `c` with length equal to  `na + nb - 1`. 

\param[in]  a   Pointer to the first vector `a`.<BR>
                Vector size is `[na x 1]`.<BR><BR>

\param[in]  na  Size of the first vector `a`.<BR><BR>

\param[in]  b   Pointer to the second vector `b`.<BR>
                Vector size is `[nb x 1]`.<BR><BR>

\param[in]  nb  Size of the second vector `b`.<BR><BR>

\param[out] c   Pointer to the convolution output vector  \f$ c = a * b\f$.<BR>
                Vector size is `[na + nb - 1  x  1]`.<BR>
                Memory must be allocated.<BR><BR>

\return `RES_OK` if convolution is calculated successfully.<BR>
Else \ref ERROR_CODE_GROUP "code error".

\note If vectors `a` and `b` are coefficients of two polynomials,
then convolution of the vectors `a` and `b` returns polynomial product 
coefficients.

Example:
\code{.cpp}
  double ar[3] = {1.0, 2.0, 3.0};
  double br[4] = {3.0, -1.0, 2.0, 4.0};
  double cr[6];
  
  int n;
  
  conv(ar, 3, br, 4, cr);
  
  for(n = 0; n < 6; n++)
    printf("cr[%d] = %5.1f\n", n, cr[n]);

\endcode
<BR>

Output:
\verbatim
cr[0] =   3.0
cr[1] =   5.0
cr[2] =   9.0
cr[3] =   5.0
cr[4] =  14.0
cr[5] =  12.0
\endverbatim

\author Sergey Bakhurin www.dsplib.org
*******************************************************************************/
int DSPL_API conv(double* a, int na, double* b, int nb, double* c)
{
  int k;
  int n;

  double *t;
  size_t bufsize;

  if(!a || !b || !c)
    return ERROR_PTR;
  if(na < 1 || nb < 1)
    return ERROR_SIZE;


  bufsize = (na + nb - 1) * sizeof(double);

  if((a != c) && (b != c))
    t = c;
  else
    t = (double*)malloc(bufsize);

  memset(t, 0, bufsize);

  for(k = 0; k < na; k++)
    for(n = 0; n < nb; n++)
      t[k+n] += a[k]*b[n];

  if(t!=c)
  {
    memcpy(c, t, bufsize);
    free(t);
  }
  return RES_OK;
}






/******************************************************************************
\ingroup FILTER_CONV_GROUP
\fn int conv_cmplx(complex_t* a, int na, complex_t* b, int nb, complex_t* c) 
\brief Complex vectors linear convolution.

Function convolves two complex vectors \f$ c = a * b\f$ length `na` and `nb`.
The output convolution is a vector `c` with length equal to  `na + nb - 1`. 

\param[in]  a   Pointer to the first vector `a`.<BR>
                Vector size is `[na x 1]`.<BR><BR>

\param[in]  na  Size of the first vector `a`.<BR><BR>

\param[in]  b   Pointer to the second vector `b`.<BR>
                Vector size is `[nb x 1]`.<BR><BR>

\param[in]  nb  Size of the second vector `b`.<BR><BR>

\param[out] c   Pointer to the convolution output vector  \f$ c = a * b\f$.<BR>
                Vector size is `[na + nb - 1  x  1]`.<BR>
                Memory must be allocated.<BR><BR>

\return `RES_OK` if convolution is calculated successfully.<BR>
Else \ref ERROR_CODE_GROUP "code error".

\note If vectors `a` and `b` are coefficients of two polynomials,
then convolution of the vectors `a` and `b` returns polynomial product 
coefficients.

Example:
\code{.cpp}
  complex_t ac[3] = {{0.0, 1.0}, {1.0, 1.0}, {2.0, 2.0}};
  complex_t bc[4] = {{3.0, 3.0}, {4.0, 4.0}, {5.0, 5.0}, {6.0, 6.0}};
  complex_t cc[6];
  
  int n;
  
  conv_cmplx(ac, 3, bc, 4, cc);
  
  for(n = 0; n < 6; n++)
    printf("cc[%d] = %5.1f%+5.1fj\n", n, RE(cc[n]),IM(cc[n]));
  
\endcode
<BR>

Output:
\verbatim
cc[0] =  -3.0 +3.0j
cc[1] =  -4.0+10.0j
cc[2] =  -5.0+25.0j
cc[3] =  -6.0+32.0j
cc[4] =   0.0+32.0j
cc[5] =   0.0+24.0j
\endverbatim

\author Sergey Bakhurin www.dsplib.org
*******************************************************************************/
int DSPL_API conv_cmplx(complex_t* a, int na, complex_t* b,
                        int nb, complex_t* c)
{
  int k;
  int n;

  complex_t *t;
  size_t bufsize;

  if(!a || !b || !c)
    return ERROR_PTR;
  if(na < 1 || nb < 1)
    return ERROR_SIZE;

  bufsize = (na + nb - 1) * sizeof(complex_t);

  if((a != c) && (b != c))
    t = c;
  else
    t = (complex_t*)malloc(bufsize);

  memset(t, 0, bufsize);

  for(k = 0; k < na; k++)
  {
    for(n = 0; n < nb; n++)
          {
      RE(t[k+n]) += CMRE(a[k], b[n]);
      IM(t[k+n]) += CMIM(a[k], b[n]);
    }
  }

  if(t!=c)
  {
    memcpy(c, t, bufsize);
    free(t);
  }

    return RES_OK;
}






/*******************************************************************************
 Complex vectors FFT linear convolution
 ******************************************************************************/
int DSPL_API conv_fft_cmplx(complex_t* a, int na, complex_t* b, int nb,
                            fft_t* pfft, complex_t* c)
{
  complex_t *pa = NULL;
  complex_t *pb = NULL;
  complex_t *pc = NULL;
  complex_t *pA = NULL;
  complex_t *pB = NULL;
  complex_t *pC = NULL;

  int nfft, nfft2, n, npos, err;
  int ma, mb;
  complex_t *ta, *tb;

  if(!a || !b || !c)
    return ERROR_PTR;
  if(na < 1 || nb < 1)
    return ERROR_SIZE;


  if(na > nb)
  {
    ma = na;
    mb = nb;
    ta = a;
    tb = b;
  }
  else
  {
    ma = nb;
    mb = na;
    ta = b;
    tb = a;
  }
  if(ma > 2*mb)
  {
    nfft = 4;
    n = mb-1;
    while(n>>=1)
      nfft <<= 1;
    nfft2 = nfft >> 1;

    pa = (complex_t*)malloc(nfft * sizeof(complex_t));
    pb = (complex_t*)malloc(nfft * sizeof(complex_t));
    pc = (complex_t*)malloc(nfft * sizeof(complex_t));
    pA = (complex_t*)malloc(nfft * sizeof(complex_t));
    pB = (complex_t*)malloc(nfft * sizeof(complex_t));
    pC = (complex_t*)malloc(nfft * sizeof(complex_t));

    npos = -nfft2;
    memset(pa, 0, nfft*sizeof(complex_t));
    memset(pb, 0, nfft*sizeof(complex_t));

    memcpy(pa + nfft2, ta, nfft2 * sizeof(complex_t));
    memcpy(pb,     tb,    mb * sizeof(complex_t));

    err = fft_cmplx(pa, nfft, pfft, pA);
    if(err != RES_OK)
      goto exit_label;

    err = fft_cmplx(pb, nfft, pfft, pB);
    if(err != RES_OK)
      goto exit_label;

    for(n = 0; n < nfft; n++)
    {
      RE(pC[n]) = CMRE(pA[n], pB[n]);
      IM(pC[n]) = CMIM(pA[n], pB[n]);
    }

    err = ifft_cmplx(pC, nfft, pfft, pc);
    if(err != RES_OK)
      goto exit_label;

    memcpy(c, pc+nfft2, nfft2*sizeof(complex_t));

    npos = 0;
    while(npos < ma)
    {
      if(npos+nfft > ma)
      {
        memset(pa, 0, nfft * sizeof(complex_t));
        memcpy(pa, ta+npos, (ma - npos) * sizeof(complex_t));
        err = fft_cmplx(pa, nfft, pfft, pA);


      }
      else
        err = fft_cmplx(ta+npos, nfft, pfft, pA);
      if(err != RES_OK)
        goto exit_label;
      for(n = 0; n < nfft; n++)
      {
        RE(pC[n]) = CMRE(pA[n], pB[n]);
        IM(pC[n]) = CMIM(pA[n], pB[n]);
      }

      err = ifft_cmplx(pC, nfft, pfft, pc);
      if(err != RES_OK)
        goto exit_label;
      if(npos+nfft <= ma+mb-1)
        memcpy(c+npos+nfft2, pc+nfft2,
            nfft2*sizeof(complex_t));
      else
      {
        if(ma+mb-1-npos-nfft2 > 0)
        {
          memcpy(c+npos+nfft2, pc+nfft2,(ma+mb-1-npos-nfft2)*sizeof(complex_t));
        }
      }
      npos+=nfft2;
    }
  }
  else
  {
    nfft = 4;
    n = ma - 1;
    while(n>>=1)
      nfft <<= 1;

    pa = (complex_t*)malloc(nfft * sizeof(complex_t));
    pb = (complex_t*)malloc(nfft * sizeof(complex_t));
    pc = (complex_t*)malloc(nfft * sizeof(complex_t));
    pA = (complex_t*)malloc(nfft * sizeof(complex_t));
    pB = (complex_t*)malloc(nfft * sizeof(complex_t));
    pC = (complex_t*)malloc(nfft * sizeof(complex_t));


    memset(pa, 0, nfft*sizeof(complex_t));
    memset(pb, 0, nfft*sizeof(complex_t));

    memcpy(pa, ta, ma * sizeof(complex_t));
    memcpy(pb, tb, mb * sizeof(complex_t));

    err = fft_cmplx(pa, nfft, pfft, pA);
    if(err != RES_OK)
      goto exit_label;
    err = fft_cmplx(pb, nfft, pfft, pB);
    if(err != RES_OK)
      goto exit_label;
    for(n = 0; n < nfft; n++)
    {
      RE(pC[n]) = CMRE(pA[n], pB[n]);
      IM(pC[n]) = CMIM(pA[n], pB[n]);
    }
    err = ifft_cmplx(pC, nfft, pfft, pc);
    if(err != RES_OK)
      goto exit_label;
    memcpy(c, pc, (ma+mb-1)*sizeof(complex_t));
  }
exit_label:
  if(pa)  free(pa);
  if(pb)  free(pb);
  if(pc)  free(pc);
  if(pA)  free(pA);
  if(pB)  free(pB);
  if(pB)  free(pC);

  return err;
}





/*******************************************************************************
IIR FILTER for real vector
*******************************************************************************/
int DSPL_API filter_iir(double* b, double* a, int ord,
                        double* x, int n, double* y)
{
  double* buf = NULL;
  double* an  = NULL;
  double  u;
  int   k;
  int   m;
  int   count;

  if(!b || !x || !y)
    return  ERROR_PTR;

  if(ord < 1 || n < 1)
      return ERROR_SIZE;

  if(a && a[0]==0.0)
    return ERROR_FILTER_A0;

  count = ord + 1;
  buf = (double*) malloc(count*sizeof(double));
  an =  (double*) malloc(count*sizeof(double));

  memset(buf, 0, count*sizeof(double));

  if(!a)
    memset(an, 0, count*sizeof(double));
  else
    for(k = 0; k < count; k++)
      an[k] = a[k] / a[0];

  for(k = 0; k < n; k++)
  {
    for(m = ord; m > 0; m--)
      buf[m] = buf[m-1];
    u = 0.0;
    for(m = ord; m > 0; m--)
      u += buf[m]*an[m];

    buf[0] = x[k] - u;
    y[k] = 0.0;
    for(m = 0; m < count; m++)
      y[k] += buf[m] * b[m];
  }

  if(buf)
    free(buf);
  if(an)
    free(an);
  return RES_OK;
}

