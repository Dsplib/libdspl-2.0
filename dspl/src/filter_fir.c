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
#include <stdio.h>
#include <string.h>
#include "dspl.h"
#include "dspl_internal.h"


/******************************************************************************
 * Linear phase lowpass filter
 ******************************************************************************/
int fir_linphase_lpf(int ord, double wp, int win_type, 
                     double win_param, double* h)
{
  int n, err = RES_OK;
  double *w = NULL;
  
  
  w = (double*)malloc((ord+1)*sizeof(double));
  
  err = linspace(-(double)ord*0.5, (double)ord*0.5, ord+1, DSPL_SYMMETRIC, w);
  
  if(err!=RES_OK)
    goto error_proc;
  
  err = sinc(w, ord+1, M_PI*wp, h);
  
  if(err!=RES_OK)
    goto error_proc;

  err = window(w, ord+1, win_type | DSPL_SYMMETRIC, win_param);
  
  if(err!=RES_OK)
    goto error_proc;
  
  for(n = 0; n < ord+1; n++)
    h[n] *= w[n] * wp;
    
error_proc:
  if(w)
    free(w);
  return err;
}




/******************************************************************************
 * Linear phase FIR filter
 ******************************************************************************/
int DSPL_API  fir_linphase(int ord, double w0, double w1, int filter_type, 
                           int win_type, double win_param, double* h)
{
  int n, err;
  double wc, b, del;
    
  if(ord<1)
    return ERROR_FILTER_ORD;
  if(w0 <= 0.0)
    return ERROR_FILTER_WP;
  if(!h)
    return ERROR_PTR;
  
  switch(filter_type & DSPL_FILTER_TYPE_MASK)
  {
    /* Lowpass FIR coefficients calculation */
    case DSPL_FILTER_LPF:
      err = fir_linphase_lpf(ord, w0, win_type, win_param, h);
      break;
    
    /* Highpass FIR coefficients calculation */
    case DSPL_FILTER_HPF:
      err = fir_linphase_lpf(ord, 1.0-w0, win_type, win_param, h);
      if(err == RES_OK)
      {
        /* LPF filter frequency inversion */
        for(n = 0; n < ord+1; n+=2)
          h[n] = -h[n];
      }
      break;
    
    /* Bandpass FIR coefficients calculation */
    case DSPL_FILTER_BPASS:
      if(w1 < w0)
      {
        err = ERROR_FILTER_WS;
        break;
      }           
      wc = (w0 + w1) * 0.5; /* central frequency  */
      b  =  w1 - w0;        /* bandwidth          */
      err = fir_linphase_lpf(ord, b*0.5, win_type, win_param, h);
      if(err == RES_OK)
      {
        /* LPF frequency shifting to the central frequency */
        del = 0.5 * (double)ord; 
        for(n = 0; n < ord+1; n++)
          h[n] *= 2.0 * cos(M_PI * ((double)n - del) * wc);
      }
      break;
      
    /* BandStop FIR coefficients calculation               */
    /* ATTENTION! Bandstop filter must be even order only! */
    case DSPL_FILTER_BSTOP:
    {  
      double *h0 = NULL;
      
      /* check filter order. Return error if order is odd. */
      if(ord%2)
        return ERROR_FILTER_ORD;
      
      /* check frequency (w1 must be higher than w0) */
      if(w1 < w0)
      {
        err = ERROR_FILTER_WS;
        break;
      }
      /* temp coeff vector */
      h0 = (double*)malloc((ord+1) * sizeof(double));
      
      /* calculate LPF */
      err = fir_linphase(ord, w0, 0.0, DSPL_FILTER_LPF, win_type, win_param, h0);
      if(err!=RES_OK)
      {
        free(h0);
        return err;        
      }      
      /* calculate HPF */      
      err = fir_linphase(ord, w1, 0.0, DSPL_FILTER_HPF, win_type, win_param, h);
      if(err==RES_OK)
      {
        /* Bandstop filter is sum of lowpass and highpass filters */
        for(n = 0; n < ord+1; n++)
          h[n] += h0[n];
      }
      free(h0);
      break;
    }
    default:
      err = ERROR_FILTER_FT;
  }  
  return err;
}
