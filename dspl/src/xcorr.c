/*
* Copyright (c) 2015-2019 Sergey Bakhurin
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
* MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.    See the
* GNU General Public License for more details.
*
* You should have received a copy of the GNU General Public License
* along with Foobar.    If not, see <http://www.gnu.org/licenses/>.
*/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "dspl.h"


int xcorr_fft_size(int nx, int ny, int* pnfft, int* pndata);

int xcorr_get_lag_cmplx(complex_t* x, int nd, int nr, complex_t* r, double* t);

int xcorr_krn(complex_t* x, int nx, complex_t* y, int ny, fft_t* pfft,
                       int flag, int nr, complex_t* r, double* t);

int xcorr_scale_cmplx(complex_t* x, int nd, int flag);



#ifdef DOXYGEN_ENGLISH
/*! ****************************************************************************
\ingroup SPEC_MATH_STAT_GROUP
\fn int xcorr(double* x, int nx, double* y, int ny, 
                   int flag, int nr, double* r, double* t)
\brief Estimates the cross-correlation vector for real 
discrete-time sequences `x` and `y`.

Estimate the cross correlation \f$\widehat{r}_{xy}(k)\f$ of vector arguments 
`x` and `y` or estimate autocorrelation vector if \f$x = y \f$.

Cross-correlation vector is defined as:

\f[
\widehat{r}_{xy}(k) = 
\frac{1}{N-k} \sum\limits_{n = 0}^{N-k-1} x(n+k)y^*(n), \qquad 0 \leq k <N;
\f]

\f[
\widehat{r}_{xy}(k) = 
\frac{1}{N+k} \sum\limits_{n = 0}^{N+k-1} y^*(n-k)x(n), \qquad  -N < k < 0.
\f]
here \f$ N = \max(n_x, n_y) \f$.

This function uses the FFT algorithm to estimate the scaled correlation vector:
\f[
\breve{r}_{xy}(k) = \operatorname{IFFT}\Big[ 
\operatorname{FFT}\big[ \mathbf{x} \big] 
\operatorname{FFT}^*\big[ \mathbf{y} \big] 
\Big]
\f]


\param[in] x
Pointer to the discrete-time vector `x`. \n
Vector size is `[nx x 1]`. \n
\n

\param[in] nx
Size of vector `x`. \n
\n

\param[in] y
Pointer to the discrete-time vector `y`. \n
Vector size is `[ny x 1]`. \n
\n

\param[in] ny
Size of vector `y`. \n
\n

\param[in] flag
Flag specifies the type of scaling applied to the correlation vector
\f$\breve{r}_{xy}(k)\f$.\n 
Is one of:\n
`DSPL_XCORR_NOSCALE` unscaled correlation vector from IFFT output 
\f$\breve{r}_{xy}(k)\f$;\n
`DSPL_XCORR_BIASED`  biased correlation vector \f$\breve{r}_{xy}(k)/N \f$;\n
`DSPL_XCORR_UNBIASED` unbiased correlation vector 
\f$\widehat{r}_{xy}(k) = \frac{\breve{r}_{xy}(k)}{N-|k|} \f$;\n
\n

\param[in] nr
Maximum correlation lag.\n
Correlation vector \f$\widehat{r}_{xy}(k)\f$ is calculated for 
\f$ k= -n_r,\,\, -n_r +1, \ldots n_r\f$.\n
\n

\param[out] r
Pointer to the cross-correlation or autocorrelation vector. \n
Vector size is `[(2*nr+1) x 1]`. \n
Memory must be allocated. \n
\n

\param[out] r
Pointer to the cross-correlation argument vector  
\f$ k= -n_r,\,\, -n_r +1, \ldors n_r\f$.\n
Vector size is `[(2*nr+1) x 1]`. \n
Pointer can be `NULL`. \n
\n

\return
`RES_OK` if function returns successfully. \n
Else \ref ERROR_CODE_GROUP "code error".

\author Sergey Bakhurin www.dsplib.org
***************************************************************************** */
#endif
#ifdef DOXYGEN_RUSSIAN
/*! ****************************************************************************
\ingroup SPEC_MATH_STAT_GROUP
\fn int xcorr(double* x, int nx, double* y, int ny, 
                   int flag, int nr, double* r, double* t)
\brief Оценка вектора взаимной корреляции для дискретных 
вещественных последовательностей `x` и `y`.

Функция производит оценку вектора взаимной корреляции \f$\widehat{r}_{xy}(k)\f$ 
для векторов `x` и `y` или вектора автокорреляции
\f$\widehat{r}_{xx}(k)\f$ если \f$x = y \f$.

Несмещенная оценка вектора взаимной корреляции:
\f[
\widehat{r}_{xy}(k) = 
\frac{1}{N-k} \sum\limits_{n = 0}^{N-k-1} x(n+k)y^*(n), \qquad 0 \leq k <N;
\f]

\f[
\widehat{r}_{xy}(k) = 
\frac{1}{N+k} \sum\limits_{n = 0}^{N+k-1} y^*(n-k)x(n), \qquad  -N < k < 0.
\f]
где \f$ N = \max(n_x, n_y) \f$.

Данная функция использует алгоритм FFT для вычислительной эффективности оценки:
\f[
\breve{r}_{xy}(k) = \operatorname{IFFT}\Big[ 
\operatorname{FFT}\big[ \mathbf{x} \big] 
\operatorname{FFT}^*\big[ \mathbf{y} \big] 
\Big]
\f]
 

\param[in] x
Указатель на первую дискретную последовательность `x`. \n
Размер вектора `[nx x 1]`. \n
\n

\param[in] nx
Размер вектора первой дискретной последовательности `x`. \n
\n

\param[in] y
Указатель на вторую дискретную последовательность `y`. \n
Размер вектора `[ny x 1]`. \n
\n

\param[in] ny
Размер вектора второй дискретной последовательности `y`. \n
\n

\param[in] flag
Флаг задает способ масштабирования выходного корреляционного вектора
\f$\breve{r}_{xy}(k)\f$.\n 
Может принимать одно из следующих значений:\n
`DSPL_XCORR_NOSCALE` немасштабированный выход алгоритма IFFT 
\f$\breve{r}_{xy}(k)\f$;\n
`DSPL_XCORR_BIASED`  Смещенная оценка \f$\breve{r}_{xy}(k)/N \f$;\n
`DSPL_XCORR_UNBIASED` Несмещенная оценка
\f$\widehat{r}_{xy}(k) = \frac{\breve{r}_{xy}(k)}{N-|k|} \f$;\n
\n

\param[in] nr
Диапазон оценки вектора корреляции относительно нуля.\n
Вектор \f$\widehat{r}_{xy}(k)\f$ рассчитывается для значений аргумента 
\f$ k= -n_r,\,\, -n_r +1, \ldots n_r\f$.\n
\n

\param[out] r
Указатель на масштабированный вектор взаимной корреляции. \n
Размер вектора `[(2*nr+1) x 1]`. \n
Память должна быть выделена. \n
\n

\param[out] t
Указатель на значения аргумента вектора взаимной корреляции  
\f$ k= -n_r,\,\, -n_r +1, \ldors n_r\f$.\n
Размер вектора `[(2*nr+1) x 1]`. \n
Указатель может быть `NULL`. В этом случае значения аргумента не возвращаются.\n
\n

\return
`RES_OK` Если функция рассчитана успешно. \n
В противном случае \ref ERROR_CODE_GROUP "код ошибки".

\author Бахурин Сергей www.dsplib.org
***************************************************************************** */
#endif
int DSPL_API xcorr(double* x, int nx, double* y, int ny, 
                   int flag, int nr, double* r, double* t)
{
    fft_t fft = {0};
    int err;
    complex_t *cx = NULL;
    complex_t *cy = NULL;
    complex_t *cr = NULL;
    
    cr = (complex_t*)malloc((2 * nr + 1) * sizeof(complex_t));
    if(!cr)
    {
        err = ERROR_MALLOC;
        goto exit_label;
    }
    cx = (complex_t*)malloc( nx * sizeof(complex_t));
    if(!cx)
    {
        err = ERROR_MALLOC;
        goto exit_label;
    }
    cy = (complex_t*)malloc( ny * sizeof(complex_t));
    if(!cy)
    {
        err = ERROR_MALLOC;
        goto exit_label;
    }
    
    err = re2cmplx(x, nx, cx);
    if(err != RES_OK)
        goto exit_label;
      
    err = re2cmplx(y, ny, cy);
    
    if(err != RES_OK)
        goto exit_label;
      
    err = xcorr_krn(cx, nx, cy, ny, &fft, flag, nr, cr, t);
    if(err != RES_OK)
        goto exit_label;

    err = cmplx2re(cr, 2*nr+1, r, NULL);

exit_label:
    if(cr)
        free(cr);
    if(cx)
        free(cx);
    if(cy)
        free(cy);
      
    fft_free(&fft);
    return err;
}





#ifdef DOXYGEN_ENGLISH
/*! ****************************************************************************
\ingroup SPEC_MATH_STAT_GROUP
int xcorr_cmplx(complex_t* x, int nx, complex_t* y, int ny, 
                         int flag, int nr, complex_t* r, double* t)
\brief Estimates the cross-correlation vector for complex 
discrete-time sequences `x` and `y`.

Estimate the cross correlation \f$\widehat{r}_{xy}(k)\f$ of vector arguments 
`x` and `y` or estimate autocorrelation vector if \f$x = y \f$.

Cross-correlation vector is defined as:

\f[
\widehat{r}_{xy}(k) = 
\frac{1}{N-k} \sum\limits_{n = 0}^{N-k-1} x(n+k)y^*(n), \qquad 0 \leq k <N;
\f]

\f[
\widehat{r}_{xy}(k) = 
\frac{1}{N+k} \sum\limits_{n = 0}^{N+k-1} y^*(n-k)x(n), \qquad  -N < k < 0.
\f]
here \f$ N = \max(n_x, n_y) \f$.

This function uses the FFT algorithm to estimate the scaled correlation vector:
\f[
\breve{r}_{xy}(k) = \operatorname{IFFT}\Big[ 
\operatorname{FFT}\big[ \mathbf{x} \big] 
\operatorname{FFT}^*\big[ \mathbf{y} \big] 
\Big]
\f]


\param[in] x
Pointer to the discrete-time vector `x`. \n
Vector size is `[nx x 1]`. \n
\n

\param[in] nx
Size of vector `x`. \n
\n

\param[in] y
Pointer to the discrete-time vector `y`. \n
Vector size is `[ny x 1]`. \n
\n

\param[in] ny
Size of vector `y`. \n
\n

\param[in] flag
Flag specifies the type of scaling applied to the correlation vector
\f$\breve{r}_{xy}(k)\f$.\n 
Is one of:\n
`DSPL_XCORR_NOSCALE` unscaled correlation vector from IFFT output 
\f$\breve{r}_{xy}(k)\f$;\n
`DSPL_XCORR_BIASED`  biased correlation vector \f$\breve{r}_{xy}(k)/N \f$;\n
`DSPL_XCORR_UNBIASED` unbiased correlation vector 
\f$\widehat{r}_{xy}(k) = \frac{\breve{r}_{xy}(k)}{N-|k|} \f$;\n
\n

\param[in] nr
Maximum correlation lag.\n
Correlation vector \f$\widehat{r}_{xy}(k)\f$ is calculated for 
\f$ k= -n_r,\,\, -n_r +1, \ldots n_r\f$.\n
\n

\param[out] r
Pointer to the cross-correlation or autocorrelation vector. \n
Vector size is `[(2*nr+1) x 1]`. \n
Memory must be allocated. \n
\n

\param[out] r
Pointer to the cross-correlation argument vector  
\f$ k= -n_r,\,\, -n_r +1, \ldors n_r\f$.\n
Vector size is `[(2*nr+1) x 1]`. \n
Pointer can be `NULL`. \n
\n

\return
`RES_OK` if function returns successfully. \n
Else \ref ERROR_CODE_GROUP "code error".

\author Sergey Bakhurin www.dsplib.org
***************************************************************************** */
#endif
#ifdef DOXYGEN_RUSSIAN
/*! ****************************************************************************
\ingroup SPEC_MATH_STAT_GROUP
int xcorr_cmplx(complex_t* x, int nx, complex_t* y, int ny, 
                         int flag, int nr, complex_t* r, double* t)
\brief Оценка вектора взаимной корреляции для дискретных 
комплексных последовательностей `x` и `y`.

Функция производит оценку вектора взаимной корреляции \f$\widehat{r}_{xy}(k)\f$ 
для векторов `x` и `y` или вектора автокорреляции
\f$\widehat{r}_{xx}(k)\f$ если \f$x = y \f$.

Несмещенная оценка вектора взаимной корреляции:
\f[
\widehat{r}_{xy}(k) = 
\frac{1}{N-k} \sum\limits_{n = 0}^{N-k-1} x(n+k)y^*(n), \qquad 0 \leq k <N;
\f]

\f[
\widehat{r}_{xy}(k) = 
\frac{1}{N+k} \sum\limits_{n = 0}^{N+k-1} y^*(n-k)x(n), \qquad  -N < k < 0.
\f]
где \f$ N = \max(n_x, n_y) \f$.

Данная функция использует алгоритм FFT для вычислительной эффективности оценки:
\f[
\breve{r}_{xy}(k) = \operatorname{IFFT}\Big[ 
\operatorname{FFT}\big[ \mathbf{x} \big] 
\operatorname{FFT}^*\big[ \mathbf{y} \big] 
\Big]
\f]
 

\param[in] x
Указатель на первую дискретную последовательность `x`. \n
Размер вектора `[nx x 1]`. \n
\n

\param[in] nx
Размер вектора первой дискретной последовательности `x`. \n
\n

\param[in] y
Указатель на вторую дискретную последовательность `y`. \n
Размер вектора `[ny x 1]`. \n
\n

\param[in] ny
Размер вектора второй дискретной последовательности `y`. \n
\n

\param[in] flag
Флаг задает способ масштабирования выходного корреляционного вектора
\f$\breve{r}_{xy}(k)\f$.\n 
Может принимать одно из следующих значений:\n
`DSPL_XCORR_NOSCALE` немасштабированный выход алгоритма IFFT 
\f$\breve{r}_{xy}(k)\f$;\n
`DSPL_XCORR_BIASED`  Смещенная оценка \f$\breve{r}_{xy}(k)/N \f$;\n
`DSPL_XCORR_UNBIASED` Несмещенная оценка
\f$\widehat{r}_{xy}(k) = \frac{\breve{r}_{xy}(k)}{N-|k|} \f$;\n
\n

\param[in] nr
Диапазон оценки вектора корреляции относительно нуля.\n
Вектор \f$\widehat{r}_{xy}(k)\f$ рассчитывается для значений аргумента 
\f$ k= -n_r,\,\, -n_r +1, \ldots n_r\f$.\n
\n

\param[out] r
Указатель на масштабированный вектор взаимной корреляции. \n
Размер вектора `[(2*nr+1) x 1]`. \n
Память должна быть выделена. \n
\n

\param[out] t
Указатель на значения аргумента вектора взаимной корреляции  
\f$ k= -n_r,\,\, -n_r +1, \ldors n_r\f$.\n
Размер вектора `[(2*nr+1) x 1]`. \n
Указатель может быть `NULL`. В этом случае значения аргумента не возвращаются.\n
\n

\return
`RES_OK` Если функция рассчитана успешно. \n
В противном случае \ref ERROR_CODE_GROUP "код ошибки".

\author Бахурин Сергей www.dsplib.org
***************************************************************************** */
#endif
int DSPL_API xcorr_cmplx(complex_t* x, int nx, complex_t* y, int ny, 
                         int flag, int nr, complex_t* r, double* t)
{
    fft_t fft = {0};
    int err;
    err = xcorr_krn(x, nx, y, ny, &fft, flag, nr, r, t);
    fft_free(&fft);
    return err;
}




#ifdef DOXYGEN_ENGLISH

#endif
#ifdef DOXYGEN_RUSSIAN

#endif
int xcorr_get_lag_cmplx(complex_t* x, int nd, int nr, complex_t* r, double* t)
{
    int i;
    if(!x || !r)
        return ERROR_PTR;
    if(nd < 1 || nr < 1)
        return ERROR_SIZE;
     
    if(nr < nd)
        memcpy(r, x+nd-1-nr, (2*nr+1)*sizeof(complex_t));
    else
    {
        memset(r, 0, (2*nr+1) * sizeof(complex_t));
        memcpy(r + nr - nd + 1, x, (2*nd-1)*sizeof(complex_t));
    }
    if(t)
        for(i = 0; i < 2*nr+1; i++)
            t[i] = (double)i - (double)nr;
    return RES_OK;
}





#ifdef DOXYGEN_ENGLISH

#endif
#ifdef DOXYGEN_RUSSIAN

#endif
int xcorr_krn(complex_t* x, int nx, complex_t* y, int ny, fft_t* pfft,
              int flag, int nr, complex_t* r, double* t)
{
    complex_t *px = NULL;
    complex_t *py = NULL;
    complex_t *pc = NULL;
    complex_t *pX = NULL;
    complex_t *pY = NULL;
    complex_t *pC = NULL;

    int nfft, ndata;
    int err, i;
    
    if(!x || !y || !r)
        return ERROR_PTR;
    if(nx < 1 || ny < 1 || nr < 1)
        return ERROR_SIZE;
    
    err = xcorr_fft_size(nx, ny, &nfft, &ndata);
    if(err!= RES_OK)
        goto exit_label;
    
    /* memory allocation */
    px = (complex_t*)malloc(nfft * sizeof(complex_t));
    if(!px)
    {
        err = ERROR_MALLOC;
        goto exit_label;
    }
    
    py = (complex_t*)malloc(nfft * sizeof(complex_t));
    if(!py)
    {
        err = ERROR_MALLOC;
        goto exit_label;
    }
    
    pc = (complex_t*)malloc(nfft * sizeof(complex_t));
    if(!pc)
    {
        err = ERROR_MALLOC;
        goto exit_label;
    }
    
    pX = (complex_t*)malloc(nfft * sizeof(complex_t));
    if(!pX)
    {
        err = ERROR_MALLOC;
        goto exit_label;
    }
    
    pY = (complex_t*)malloc(nfft * sizeof(complex_t));
    if(!pY)
    {
        err = ERROR_MALLOC;
        goto exit_label;
    }
    
    pC = (complex_t*)malloc(nfft * sizeof(complex_t));
    if(!pC)
    {
        err = ERROR_MALLOC;
        goto exit_label;
    }
    
    memset(px, 0, nfft * sizeof(complex_t));
    memset(py, 0, nfft * sizeof(complex_t));
    
    memcpy(px + ndata - 1, x, nx * sizeof(complex_t));
    memcpy(py, y, ny * sizeof(complex_t));
    
    err = fft_cmplx(px, nfft, pfft, pX);
    if(err!= RES_OK)
        goto exit_label;
    
    err = fft_cmplx(py, nfft, pfft, pY);
    if(err!= RES_OK)
        goto exit_label;
      
    for(i = 0; i < nfft; i++)
    {
        RE(pC[i]) = CMCONJRE(pX[i], pY[i]);
        IM(pC[i]) = CMCONJIM(pX[i], pY[i]);
    }
    
    err = ifft_cmplx(pC, nfft, pfft, pc);
    if(err!= RES_OK)
        goto exit_label;
      
    err = xcorr_scale_cmplx(pc, ndata, flag);
    if(err!= RES_OK)
        goto exit_label;
      
    err = xcorr_get_lag_cmplx(pc, ndata, nr, r, t);
    
exit_label:
    if(px)
        free(px);
    if(py)
        free(py);
    if(pc)
        free(pc);
    if(pX)
        free(pX);
    if(pY)
        free(pY);
    if(pC)
        free(pC);
    return err;
}




#ifdef DOXYGEN_ENGLISH
/*******************************************************************************
Return FFT size for autocorrelation or cross correlation vector calculation

Cross-correlation vector size is 
N = 2 * nx - 1,   if   nx >  ny;
N = 2 * ny - 1,   if   nx <= ny.

If cross-correlation size N may not be efficient for FFT 
then we can add zeros to get high-performance FFT size.

For example if N = 1025, then we can add zeros to 2048-points FFT but this way
seems not so good because too much zeros.

If we rewrite  N = 2^L + D, then we can use

               NFFT = 2^L + 2^(L - P), here  P = 0,1,2 or 3.

So NFFT = 2^(L-P) * (2^P + 1). Then 2^(L-P) can use radix-2 FFT, and additional
composite multiplication if P = 0,1,2 or 3 equals 
9, 5, 3 or 2, and we have high-performance FFT algorithms for its points. 
If P = 4 then composite  multiplier is (2^P + 1) = 17, has no good FFT. 
*******************************************************************************/
#endif
#ifdef DOXYGEN_RUSSIAN
/*******************************************************************************
Возвращает размер FFT для расчета полного вектора автокорреляции 
или кросскорреляции.

Размер кросскорреляции равен 
N = 2 * nx - 1,   если   nx >  ny;
N = 2 * ny - 1,   eсли   nx <= ny.

Посколку N может оказаться неудачным размером для FFT, то можно добить нулями
до удобной длины.

Если например N = 1025, то добивать до длины 2048 не очень эффективно, потому
что много лишних нулей.

Если мы рассмотрим N = 2^L + D, то целесообразно использовать 

               NFFT = 2^L + 2^(L - P), где  P = 0,1,2 или 3.

Тогда NFFT = 2^(L-P) * (2^P + 1). Тогда 2^(L-P) реализуем как radix-2, а 
дополнительный составной множитель при P = 0,1,2 или 3 равен соответсвенно 
9, 5, 3 или 2, а для этих длин существуют хорошие процедуры. 
При P = 4 составной множитель будет (2^P + 1) = 17, что не очень хорошо. 
*******************************************************************************/
#endif
int xcorr_fft_size(int nx, int ny, int* pnfft, int* pndata)
{
    int nfft, nfft2, r2, dnfft;
    
    if(nx < 1 || ny < 1)
        return ERROR_SIZE;
    if(!pnfft || !pndata)
        return ERROR_PTR;

    if(nx > ny)
    {
        nfft  = 2*nx - 1;
        *pndata = nx;
    }
    else
    {
        nfft  = 2*ny - 1;
        *pndata = ny;
    }
    nfft2 = nfft;

    r2 = 0;
    while(nfft2 >>= 1)
        r2++;
    
    if(r2 > 3)
    {
        dnfft = 1 << (r2 - 3);
        while(((1 << r2) + dnfft) < nfft)
            dnfft <<= 1;
        nfft = (1 << r2) + dnfft;
    }

    *pnfft = nfft;
    return RES_OK;
}





#ifdef DOXYGEN_ENGLISH

#endif
#ifdef DOXYGEN_RUSSIAN

#endif
int xcorr_scale_cmplx(complex_t* x, int nd, int flag)
{
    int i;
    double w;
    if(!x)
        return ERROR_PTR;
    if(nd < 1)
        return ERROR_SIZE;
      
    switch(flag)
    {
        case DSPL_XCORR_NOSCALE:
            break;
        case  DSPL_XCORR_BIASED:
            for(i = 0; i < 2 * nd - 1; i++)
            {
                w = 1.0 / (double)nd;
                RE(x[i]) *= w;
                IM(x[i]) *= w;
            }
            break;
        case  DSPL_XCORR_UNBIASED:
            for(i = 1; i < 2 * nd - 1; i++)
            {
                w = 1.0 / ((double)nd - fabs((double)(i - nd)));
                RE(x[i-1]) *= w;
                IM(x[i-1]) *= w;
            }
            break;
        default:
            return ERROR_XCORR_FLAG;
    }
    return RES_OK;
}
