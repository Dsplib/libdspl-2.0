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

\param[in]  a   Pointer to the first vector `a`. /n
                Vector size is `[na x 1]`. /n /n

\param[in]  na  Size of the first vector `a`. /n /n

\param[in]  b   Pointer to the second vector `b`. /n
                Vector size is `[nb x 1]`. /n /n

\param[in]  nb  Size of the second vector `b`. /n /n

\param[out] c   Pointer to the convolution output vector  \f$ c = a * b\f$. /n
                Vector size is `[na + nb - 1  x  1]`. /n
                Memory must be allocated. /n /n

\return `RES_OK` if convolution is calculated successfully. /n
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
 /n

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

\param[in]  a   Pointer to the first vector `a`. /n
                Vector size is `[na x 1]`. /n /n

\param[in]  na  Size of the first vector `a`. /n /n

\param[in]  b   Pointer to the second vector `b`. /n
                Vector size is `[nb x 1]`. /n /n

\param[in]  nb  Size of the second vector `b`. /n /n

\param[out] c   Pointer to the convolution output vector  \f$ c = a * b\f$. /n
                Vector size is `[na + nb - 1  x  1]`. /n
                Memory must be allocated. /n /n

\return `RES_OK` if convolution is calculated successfully. /n
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
 /n

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





/******************************************************************************
\ingroup FILTER_CONV_GROUP
\fn int conv_fft(double* a, int na, double* b, int nb,
                 fft_t* pfft, int nfft, double* c) 
\brief Real vectors fast linear convolution by using fast Fourier
transform algorithms

Function convolves two real vectors \f$ c = a * b\f$ length `na` and `nb`
in the frequency domain by using FFT algorithms. This approach provide 
high-performance convolution which increases with `na` and `nb` increasing.
The output convolution is a vector `c` with length equal to  `na + nb - 1`. 

\param[in]  a     Pointer to the first vector `a`. /n
                  Vector size is `[na x 1]`. /n /n

\param[in]  na    Size of the first vector `a`. /n /n

\param[in]  b     Pointer to the second vector `b`. /n
                  Vector size is `[nb x 1]`. /n /n

\param[in]  nb    Size of the second vector `b`. /n /n

\param[in]  pfft  Pointer to the structure `fft_t`. /n
                  Function changes `fft_t` structure fields so `fft_t` must
                  be clear before program returns. /n /n
                  
\param[in] nfft   FFT size.  /n
                  This parameter set which FFT size will be used 
                  for overlapped frequency domain convolution. /n
                  FFT size must be more of minimal `na` and `nb` value.
                  For example if `na = 10`, `nb = 4` then `nfft` parameter must 
                  be more than 4.  /n

\param[out] c     Pointer to the convolution output vector  \f$ c = a * b\f$. /n
                  Vector size is `[na + nb - 1  x  1]`. /n
                  Memory must be allocated. /n /n

\return `RES_OK` if convolution is calculated successfully. /n
Else \ref ERROR_CODE_GROUP "code error".  /n /n

Example:
\include conv_fft_test.c

Program output:

\verbatim
conv_fft error: 0x00000000
conv error:     0x00000000
c[  0] =     -0.00    d[  0] =      0.00
c[  1] =     -0.00    d[  1] =      0.00
c[  2] =      1.00    d[  2] =      1.00
c[  3] =      4.00    d[  3] =      4.00
c[  4] =     10.00    d[  4] =     10.00
c[  5] =     20.00    d[  5] =     20.00
c[  6] =     35.00    d[  6] =     35.00
c[  7] =     56.00    d[  7] =     56.00
c[  8] =     77.00    d[  8] =     77.00
c[  9] =     98.00    d[  9] =     98.00
c[ 10] =    119.00    d[ 10] =    119.00
c[ 11] =    140.00    d[ 11] =    140.00
c[ 12] =    161.00    d[ 12] =    161.00
c[ 13] =    182.00    d[ 13] =    182.00
c[ 14] =    190.00    d[ 14] =    190.00
c[ 15] =    184.00    d[ 15] =    184.00
c[ 16] =    163.00    d[ 16] =    163.00
c[ 17] =    126.00    d[ 17] =    126.00
c[ 18] =     72.00    d[ 18] =     72.00
\endverbatim

\author Sergey Bakhurin www.dsplib.org
*******************************************************************************/
int DSPL_API conv_fft(double* a, int na,   double* b, int nb,
                      fft_t* pfft,  int nfft, double* c)
{
  complex_t *pa = NULL, *pb = NULL, *pc = NULL;
  int err;
  
  if(!a || !b || !c || !pfft)
    return ERROR_PTR;
  if(na<1 || nb < 1)
    return ERROR_SIZE;
  if(nfft<2)
    return ERROR_FFT_SIZE;
  
  pa = (complex_t*) malloc(na*sizeof(complex_t));
  pb = (complex_t*) malloc(nb*sizeof(complex_t));
  pc = (complex_t*) malloc((na+nb-1)*sizeof(complex_t));
  
  re2cmplx(a, na, pa);
  re2cmplx(b, nb, pb);
  
  err = conv_fft_cmplx(pa, na, pb, nb, pfft, nfft, pc);
  if(err != RES_OK)
    goto exit_label;
  
  err = cmplx2re(pc, na+nb-1, c, NULL);
  
exit_label:
  if(pa) free(pa);
  if(pb) free(pb);
  if(pc) free(pc);
  
  return err;
}




/******************************************************************************
\ingroup FILTER_CONV_GROUP
\fn int conv_fft_cmplx(complex_t* a, int na, complex_t* b, int nb,
                       fft_t* pfft, int nfft, complex_t* c) 
\brief Complex vectors fast linear convolution by using fast Fourier
transform algorithms

Function convolves two complex vectors \f$ c = a * b\f$ length `na` and `nb`
in the frequency domain by using FFT algorithms. This approach provide 
high-performance convolution which increases with `na` and `nb` increasing.
The output convolution is a vector `c` with length equal to  `na + nb - 1`. 

\param[in]  a   Pointer to the first vector `a`. /n
                Vector size is `[na x 1]`. /n /n

\param[in]  na  Size of the first vector `a`. /n /n

\param[in]  b   Pointer to the second vector `b`. /n
                Vector size is `[nb x 1]`. /n /n

\param[in]  nb  Size of the second vector `b`. /n /n

\param[in]  pfft  Pointer to the structure `fft_t`. /n
                  Function changes `fft_t` structure fields so `fft_t` must
                  be clear before program returns. /n /n
                  
\param[in] nfft   FFT size.  /n
                  This parameter set which FFT size will be used 
                  for overlapped frequency domain convolution. /n
                  FFT size must be more of minimal `na` and `nb` value.
                  For example if `na = 10`, `nb = 4` then `nfft` parameter must 
                  be more than 4.  /n

\param[out] c   Pointer to the convolution output vector  \f$ c = a * b\f$. /n
                Vector size is `[na + nb - 1  x  1]`. /n
                Memory must be allocated. /n /n

\return `RES_OK` if convolution is calculated successfully. /n
Else \ref ERROR_CODE_GROUP "code error".  /n /n

Example:
\include conv_fft_cmplx_test.c

Program output:

\verbatim
c[  0] =     -1.00    -0.00j    d[  0] =     -1.00    +0.00j
c[  1] =     -6.00    +4.00j    d[  1] =     -6.00    +4.00j
c[  2] =    -15.00   +20.00j    d[  2] =    -15.00   +20.00j
c[  3] =    -28.00   +56.00j    d[  3] =    -28.00   +56.00j
c[  4] =    -45.00  +120.00j    d[  4] =    -45.00  +120.00j
c[  5] =    -55.00  +210.00j    d[  5] =    -55.00  +210.00j
c[  6] =    -65.00  +300.00j    d[  6] =    -65.00  +300.00j
c[  7] =    -75.00  +390.00j    d[  7] =    -75.00  +390.00j
c[  8] =    -85.00  +480.00j    d[  8] =    -85.00  +480.00j
c[  9] =    -95.00  +570.00j    d[  9] =    -95.00  +570.00j
c[ 10] =   -105.00  +660.00j    d[ 10] =   -105.00  +660.00j
c[ 11] =   -115.00  +750.00j    d[ 11] =   -115.00  +750.00j
c[ 12] =   -125.00  +840.00j    d[ 12] =   -125.00  +840.00j
c[ 13] =   -135.00  +930.00j    d[ 13] =   -135.00  +930.00j
c[ 14] =   -145.00 +1020.00j    d[ 14] =   -145.00 +1020.00j
c[ 15] =   -124.00 +1080.00j    d[ 15] =   -124.00 +1080.00j
c[ 16] =    -99.00 +1016.00j    d[ 16] =    -99.00 +1016.00j
c[ 17] =    -70.00  +820.00j    d[ 17] =    -70.00  +820.00j
c[ 18] =    -37.00  +484.00j    d[ 18] =    -37.00  +484.00j
\endverbatim

\author Sergey Bakhurin www.dsplib.org
*******************************************************************************/
int DSPL_API conv_fft_cmplx(complex_t* a, int na,   complex_t* b, int nb,
                            fft_t* pfft,  int nfft, complex_t* c)
{
  
  int La, Lb, Lc, Nz, n, p0, p1, ind, err;
  complex_t *pa, *pb;
  complex_t *pt, *pA, *pB, *pC;
  
  if(!a || !b || !c)
    return ERROR_PTR;
  if(na < 1 || nb < 1)
    return ERROR_SIZE;
  
  if(na >= nb)
  {
    La = na;
    Lb = nb;
    pa = a; 
    pb = b;
  }
  else
  {
    La = nb;
    pa = b;
    Lb = na;
    pb = a;
  }
    
  Lc = La + Lb - 1;
  Nz = nfft - Lb;

  if(Nz <= 0)
    return ERROR_FFT_SIZE;

  pt = (complex_t*)malloc(nfft*sizeof(complex_t));
  pB = (complex_t*)malloc(nfft*sizeof(complex_t));  
  pA = (complex_t*)malloc(nfft*sizeof(complex_t));  
  pC = (complex_t*)malloc(nfft*sizeof(complex_t));

  memset(pt,    0,  nfft*sizeof(complex_t));
  memcpy(pt+Nz, pb, Lb*sizeof(complex_t));

  err = fft_cmplx(pt, nfft, pfft, pB);
  if(err != RES_OK)
    goto exit_label;

  p0 = -Lb;
  p1 = p0 + nfft;
  ind = 0;
  while(ind < Lc)
  {
    if(p0 >=0)
    {
      if(p1 < La)
        err = fft_cmplx(pa + p0, nfft, pfft, pA);
      else
      {
        memset(pt, 0, nfft*sizeof(complex_t));
        memcpy(pt, pa+p0, (nfft+La-p1)*sizeof(complex_t));
        err = fft_cmplx(pt, nfft, pfft, pA);
      }
    }
    else
    {
      memset(pt, 0, nfft*sizeof(complex_t));
      if(p1 < La)        
        memcpy(pt - p0, pa, (nfft+p0)*sizeof(complex_t));
      else
        memcpy(pt - p0, pa, La * sizeof(complex_t));
      err = fft_cmplx(pt, nfft, pfft, pA);
    }
    
    if(err != RES_OK)
      goto exit_label;

    for(n = 0; n < nfft; n++)
    {
      RE(pC[n]) = CMRE(pA[n], pB[n]);
      IM(pC[n]) = CMIM(pA[n], pB[n]);
    }


    if(ind+nfft < Lc)
      err = ifft_cmplx(pC, nfft, pfft, c+ind);
    else
    {
      err = ifft_cmplx(pC, nfft, pfft, pt);
      memcpy(c+ind, pt, (Lc-ind)*sizeof(complex_t));
    }
    if(err != RES_OK)
      goto exit_label;
    
    p0  += Nz;
    p1  += Nz;
    ind += Nz;
  }
 
exit_label: 
  if(pt) free(pt);
  if(pB) free(pB);
  if(pA) free(pA);
  if(pC) free(pC);
  
  return err;
}



/*******************************************************************************
\ingroup FILTER_CONV_GROUP
\fn int filter_iir(double* b, double* a, int ord, double* x, int n, double* y)
\brief Real IIR filtration

Function calculates real IIR filter output for real signal. The real filter
contains real coefficients of the transfer function \f$H(z)\f$
 numerator and denominator:
\f[
  H(z) = \frac{\sum_{n = 0}^{N} b_n  z^{-n}}
  {1+{\frac{1}{a_0}}\sum_{m = 1}^{M} a_m  z^{-n}},
\f]
here \f$a_0\f$ cannot be equals zeros, \f$N=M=\f$`ord`.


\param[in]  b     Pointer to the vector \f$b\f$ of IIR filter 
                  transfer function numerator coefficients. /n 
                  Vector size is `[ord + 1 x 1]`. /n /n 

\param[in]  a     Pointer to the vector \f$a\f$ of IIR filter 
                  transfer function denominator coefficients. /n 
                  Vector size is `[ord + 1 x 1]`. /n 
                  This pointer can be `NULL` if filter is FIR. /n /n 

\param[in]  ord   Filter order. Number of the transfer function 
                  numerator and denominator coefficients 
                  (length of vectors `b` and `a`) is `ord + 1`. /n /n 

\param[in]  x     Pointer to the input signal vector. /n 
                  Vector size is `[n x 1]`. /n /n 

\param[in]  n     Size of the input signal vector `x`. /n /n 

\param[out] y     Pointer to the IIR filter output vector. /n 
                  Vector size is `[n x  1]`. /n 
                  Memory must be allocated. /n /n 
\return
`RES_OK` if filter output is calculated successfully. /n
Else \ref ERROR_CODE_GROUP "code error": /n

Example:

\include filter_iir_test.c

Input signal is
\f$s(t) = \sin(2\pi \cdot 0.05 t) + n(t)\f$, here \f$n(t)\f$ white Gaussian
noise with zero mean value and unit standard deviation. \n

Input signal is filtered by elliptic LPF order 6 and output signal and data 
saves in the txt-files  

\verbatim
dat/s.txt  - input signal + noise
dat/sf.txt - filter output.
\endverbatim

Plots:

\image html filter_iir_test.png

GNUPLOT script for make plots is:
\include filter_iir.plt

\author Sergey Bakhurin www.dsplib.org
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

