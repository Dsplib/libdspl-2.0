/*
* Copyright (c) 2015-2019 Sergey Bakhurin
* Digital Signal Processing Library [http://dsplib.org]
*
* This file is part of libdspl-2.0.
*
* is free software: you can redistribute it and/or modify
* it under the terms of the GNU Lesser    General Public License as published by
* the Free Software Foundation, either version 3 of the License, or
* (at your option) any later version.
*
* DSPL is distributed in the hope that it will be useful,
* but WITHOUT ANY WARRANTY; without even the implied warranty of
* MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.    See the
* GNU General Public License for more details.
*
* You should have received a copy of the GNU Lesser General Public License
* along with Foobar.    If not, see <http://www.gnu.org/licenses/>.
*/



#include <stdio.h>
#include <math.h>
#include "dspl.h"
#include "dspl_internal.h"



/*******************************************************************************
Window function
*******************************************************************************/
int window(double* w, int n, int win_type, double param)
{
    switch(win_type & DSPL_WIN_MASK)
    {
        case DSPL_WIN_BARTLETT:
            return win_bartlett(w, n, win_type);
        case DSPL_WIN_BARTLETT_HANN:
            return win_bartlett_hann(w, n, win_type);
        case DSPL_WIN_BLACKMAN:
            return win_blackman(w, n, win_type);
        case DSPL_WIN_BLACKMAN_HARRIS:
            return win_blackman_harris(w, n, win_type);
        case DSPL_WIN_BLACKMAN_NUTTALL:
            return win_blackman_nuttall(w, n, win_type);
        case DSPL_WIN_CHEBY:
            return win_cheby(w, n, param);
        case DSPL_WIN_FLAT_TOP:
            return win_flat_top(w, n, win_type);
        case DSPL_WIN_GAUSSIAN:
            return win_gaussian(w, n, win_type, param);
        case DSPL_WIN_HAMMING:
            return win_hamming(w, n, win_type);
        case DSPL_WIN_HANN:
            return win_hann(w, n, win_type);
        case DSPL_WIN_KAISER:
            return win_kaiser(w, n, win_type, param);
        case DSPL_WIN_LANCZOS:
            return win_lanczos(w, n, win_type);
        case DSPL_WIN_NUTTALL:
            return win_nuttall(w, n, win_type);
        case DSPL_WIN_RECT:
            return win_rect(w, n);
        case DSPL_WIN_COS:
            return win_cos(w, n, win_type);
        default:
            return ERROR_WIN_TYPE;
    }
    return RES_OK;
}



/******************************************************************************
Barlett window function
*******************************************************************************/
int win_bartlett(double *w, int n, int win_type)
{
    double x = 0.0;
    int i;

    if(!w)
        return ERROR_PTR;
    if(n<2)
        return ERROR_SIZE;

    switch(win_type & DSPL_WIN_SYM_MASK)
    {
        case DSPL_WIN_SYMMETRIC: x = (double)(n-1); break;
        case DSPL_WIN_PERIODIC : x = (double)n;     break;
        default: return ERROR_WIN_SYM;
    }

    for(i = 0; i < n; i++)
    {
        w[i] = 2.0 / x * (x * 0.5-fabs((double)i - x * 0.5));
    }
    return RES_OK;
}





/******************************************************************************
Barlett - Hann    window function
******************************************************************************/
int win_bartlett_hann(double *w, int n, int win_type)
{
    double y;
    double x = 0.0;
    int i;

    if(!w)
    return ERROR_PTR;
    if(n<2)
        return ERROR_SIZE;

    switch(win_type & DSPL_WIN_SYM_MASK)
    {
        case DSPL_WIN_SYMMETRIC: x = 1.0/(double)(n-1); break;
        case DSPL_WIN_PERIODIC : x = 1.0/(double)n; break;
        default: return ERROR_WIN_SYM;
    }

    y = 0.0;
    for(i = 0; i<n; i++)
    {
        w[i] = 0.62 - 0.48 * fabs(y-0.5)-0.38*cos(M_2PI*y);
        y += x;
    }
    return RES_OK;
}





/******************************************************************************
Blackman    window function
******************************************************************************/
int win_blackman(double *w, int n, int win_type)
{
    double y;
    double x = 0.0;
    int i;

    if(!w)
    return ERROR_PTR;
    if(n<2)
        return ERROR_SIZE;

    switch(win_type & DSPL_WIN_SYM_MASK)
    {
        case DSPL_WIN_SYMMETRIC: x = M_2PI/(double)(n-1); break;
        case DSPL_WIN_PERIODIC : x = M_2PI/(double)n; break;
        default: return ERROR_WIN_SYM;
    }

    y = 0.0;
    for(i = 0; i<n; i++)
    {
        w[i] = 0.42 - 0.5* cos(y)+0.08*cos(2.0*y);
        y += x;
    }
    return RES_OK;
}




/******************************************************************************
Blackman - Harris window function
******************************************************************************/
int win_blackman_harris(double *w, int n, int win_type)
{
    double y;
    double x  = 0.0;
    double a0 = 0.35875;
    double a1 = 0.48829;
    double a2 = 0.14128;
    double a3 = 0.01168;
    int i;

    if(!w)
        return ERROR_PTR;
    if(n<2)
        return ERROR_SIZE;

    switch(win_type & DSPL_WIN_SYM_MASK)
    {
        case DSPL_WIN_SYMMETRIC: x = M_2PI/(double)(n-1); break;
        case DSPL_WIN_PERIODIC : x = M_2PI/(double)n; break;
        default: return ERROR_WIN_SYM;
    }

    y = 0.0;
    for(i = 0; i<n; i++)
    {
        w[i] = a0 - a1* cos(y)+a2*cos(2.0*y)-a3*cos(3.0*y);
        y += x;
    }
    return RES_OK;
}




/******************************************************************************
Blackman - Nuttull     window function
******************************************************************************/
int win_blackman_nuttall(double *w, int n, int win_type)
{
    double y;
    double x    = 0.0;
    double a0 = 0.3635819;
    double a1 = 0.4891775;
    double a2 = 0.1365995;
    double a3 = 0.0106411;
    int i;


    if(!w)
        return ERROR_PTR;
    if(n<2)
        return ERROR_SIZE;

    switch(win_type & DSPL_WIN_SYM_MASK)
    {
        case DSPL_WIN_SYMMETRIC:        x = M_2PI/(double)(n-1);        break;
        case DSPL_WIN_PERIODIC :        x = M_2PI/(double)n;                break;
        default: return ERROR_WIN_SYM;
    }

    y = 0.0;
    for(i = 0; i<n; i++)
    {
        w[i] = a0 - a1* cos(y)+a2*cos(2.0*y)-a3*cos(3.0*y);
        y += x;
    }
    return RES_OK;
}






/******************************************************************************
Chebyshev parametric window function
param sets spectrum sidelobes level in dB
ATTENTION! ONLY SYMMETRIC WINDOW
*******************************************************************************/
int win_cheby(double *w, int n, double param)
{
    int k, i, m;
    double z, dz, sum = 0, wmax=0, r1, x0, chx, chy, in;

    if(!w)
        return ERROR_PTR;

    if(n<2)
        return ERROR_SIZE;

    if(param <= 0.0)
        return ERROR_WIN_PARAM;

    r1 = pow(10, param/20);
    x0 = cosh((1.0/(double)(n-1)) * acosh(r1));

    /* check window length even or odd */
    if(n%2==0)
    {
        dz = 0.5;
        m = n/2-1;
    }
    else
    {
        m = (n-1)/2;
        dz = 0.0;
    }

    for(k = 0; k < m+2; k++)
    {
        z = (double)(k - m) - dz;
        sum = 0;

        for(i = 1; i <= m; i++)
        {
            in = (double)i / (double)n;
            chx = x0 * cos(M_PI * in);
            cheby_poly1(&chx, 1, n-1, &chy);
            sum += chy * cos(2.0 * z * M_PI * in);
        }

        w[k] = r1 + 2.0 * sum;
        w[n-1-k] = w[k];

        /* max value calculation */
        if(w[k]>wmax)
            wmax=w[k];
    }

    /* normalization */
    for(k=0; k < n; k++)
        w[k] /= wmax;

    return RES_OK;
}



/******************************************************************************
Cosine window function
******************************************************************************/
int win_cos(double *w, int n, int win_type)
{
    double y;
    double x = 0.0;
    int i;

    if(!w)
        return ERROR_PTR;
    if(n<2)
        return ERROR_SIZE;

    switch(win_type & DSPL_WIN_SYM_MASK)
    {
        case DSPL_WIN_SYMMETRIC: x = M_PI/(double)(n-1); break;
        case DSPL_WIN_PERIODIC : x = M_PI/(double)n; break;
        default: return ERROR_WIN_SYM;
    }

    y = 0.0;
    for(i = 0; i<n; i++)
    {
        w[i] = sin(y);
        y += x;
    }
    return RES_OK;
}






/******************************************************************************
Flat - Top     window function
******************************************************************************/
int win_flat_top(double *w, int n, int win_type)
{
    double y;
    double x  = 0.0;
    double a0 = 1.0;
    double a1 = 1.93;
    double a2 = 1.29;
    double a3 = 0.388;
    double a4 = 0.032;
    int i;

    if(!w)
        return ERROR_PTR;
    if(n<2)
        return ERROR_SIZE;

    switch(win_type & DSPL_WIN_SYM_MASK)
    {
        case DSPL_WIN_SYMMETRIC: x = M_2PI/(double)(n-1); break;
        case DSPL_WIN_PERIODIC : x = M_2PI/(double)n; break;
        default: return ERROR_WIN_SYM;
    }

    y = 0.0;
    for(i = 0; i<n; i++)
    {
        w[i] = a0 - a1* cos(y)+a2*cos(2.0*y)-a3*cos(3.0*y)+a4*cos(4.0*y);
        y += x;
    }
    return RES_OK;
}






/******************************************************************************
Gaussian window function
******************************************************************************/
int win_gaussian(double *w, int n, int win_type, double alpha)
{
    double x = 0.0;
    double y;
    double sigma;
    int i;

    if(!w)
        return ERROR_PTR;
    if(n<2)
        return ERROR_SIZE;

    switch(win_type & DSPL_WIN_SYM_MASK)
    {
        case DSPL_WIN_SYMMETRIC: x = (double)(n-1)*0.5; break;
        case DSPL_WIN_PERIODIC : x = (double)(n)*0.5; break;
        default: return ERROR_WIN_SYM;
    }


    sigma = alpha / x;
    for(i = 0; i<n; i++)
    {
        y = ((double)i - x)*sigma;
        w[i] = exp(-0.5*y*y);
    }
        return RES_OK;
}






/******************************************************************************
Hamming window function
******************************************************************************/
int win_hamming(double *w, int n, int win_type)
{
    double x = 0.0;
    double y;
    int i;

    if(!w)
        return ERROR_PTR;
    if(n<2)
        return ERROR_SIZE;

    switch(win_type & DSPL_WIN_SYM_MASK)
    {
        case DSPL_WIN_SYMMETRIC: x = M_2PI/(double)(n-1); break;
        case DSPL_WIN_PERIODIC : x = M_2PI/(double)n; break;
        default: return ERROR_WIN_SYM;
    }

    y = 0.0;
    for(i = 0; i < n; i++)
    {
        w[i] = 0.54-0.46*cos(y);
        y += x;
    }
    return RES_OK;
}




/******************************************************************************
Hann window function
******************************************************************************/
int win_hann(double *w, int n, int win_type)
{
    double x = 0.0;
    double y;
    int i;

    if(!w)
        return ERROR_PTR;
    if(n<2)
        return ERROR_SIZE;

    switch(win_type & DSPL_WIN_SYM_MASK)
    {
        case DSPL_WIN_SYMMETRIC: x = M_2PI/(double)(n-1); break;
        case DSPL_WIN_PERIODIC : x = M_2PI/(double)n; break;
        default: return ERROR_WIN_SYM;
    }

    y = 0.0;
    for(i = 0; i < n; i++)
    {
        w[i] = 0.5*(1-cos(y));
        y += x;
    }
    return RES_OK;
}


/******************************************************************************
Kaiser window function
******************************************************************************/
int win_kaiser(double* w, int n, int win_type, double param)
{
    double num, den, x, y;
    int i, err;
    if(!w)
        return ERROR_PTR;
    if(n<2)
        return ERROR_SIZE;

    switch(win_type & DSPL_WIN_SYM_MASK)
    {
        case DSPL_WIN_SYMMETRIC: x = 1.0/(double)(n-1); break;
        case DSPL_WIN_PERIODIC : x = 1.0/(double)n; break;
        default: return ERROR_WIN_SYM;
    }

    err = bessel_i0(&param, 1, &den);
    if(err != RES_OK)
        return err;
    for(i = 0; i < n; i++)
    {
        y = (double)(2*i) / x - 1.0;
        y = param * sqrt(1.0 - y*y);
        err = bessel_i0(&y, 1, &num);
        if(err != RES_OK)
            return err;
        w[i] = num / den;
    }
    return err;
}



/******************************************************************************
Lanczos window function
******************************************************************************/
int win_lanczos(double *w, int n, int win_type)
{
    double y;
    double x = 0.0;
    int i;

    if(!w)
        return ERROR_PTR;
    if(n<2)
        return ERROR_SIZE;

    switch(win_type & DSPL_WIN_SYM_MASK)
    {
        case DSPL_WIN_SYMMETRIC: x = M_2PI/(double)(n-1); break;
        case DSPL_WIN_PERIODIC : x = M_2PI/(double)n; break;
        default: return ERROR_WIN_SYM;
    }

    y = 0.0;
    for(i = 0; i < n; i++)
    {
        if((y - M_PI)==0.0)
            w[i] = 1.0;
        else
            w[i] = sin(y - M_PI)/(y - M_PI);
        y += x;
    }
    return RES_OK;

}



/******************************************************************************
Nuttall window function
******************************************************************************/
int win_nuttall(double *w, int n, int win_type)
{
    double y;
    double x  = 0.0;
    double a0 = 0.355768;
    double a1 = 0.487396;
    double a2 = 0.144232;
    double a3 = 0.012604;
    int i;

    if(!w)
        return ERROR_PTR;
    if(n<2)
        return ERROR_SIZE;

    switch(win_type & DSPL_WIN_SYM_MASK)
    {
        case DSPL_WIN_SYMMETRIC: x = M_2PI/(double)(n-1); break;
        case DSPL_WIN_PERIODIC : x = M_2PI/(double)n; break;
        default: return ERROR_WIN_SYM;
    }

    y = 0.0;
    for(i = 0; i < n; i++)
    {
        w[i] = a0 - a1* cos(y)+a2*cos(2.0*y)-a3*cos(3.0*y);
        y += x;
    }
    return RES_OK;
}





/******************************************************************************
Rectangle window function
******************************************************************************/
int win_rect(double *w, int n)
{
    int i;

    if(!w)
        return ERROR_PTR;
    if(n<2)
        return ERROR_SIZE;

    for(i = 0; i < n; i++)
        w[i] = 1.0;
    return RES_OK;
}

