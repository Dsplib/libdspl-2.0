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




/*******************************************************************************
\ingroup IIR_FILTER_DESIGN_GROUP
\fn int bilinear(double* bs, double* as, int ord, double* bz, double* az)
\brief Analog filter transfer function H(s) bilinear transform to the
       digital filter transfer function H(z).
       
Function calculates digital filter coefficients by rational substitution

\f[
s \leftarrow \frac{1 - z^{-1}}{1 - z^{-1}}.
\f]

Digital filter order still the same as analog prototype order. 
Analog prototype frequency \f$\Omega\f$ related with the  digital filter 
normalized frequency \f$\omega\f$ as:

\f[
\Omega = \tan(\omega / 2).
\f]

\param[in]  bs   Pointer to the numerator coefficients of an 
                 analog prototype transfer function \f$H(s)\f$. <BR>
                 Array size is `[ord+1 x 1]`.<BR>
                 Memory must be allocated.<BR><BR>

\param[in]  as   Pointer to the denominator coefficients of an 
                 analog prototype transfer function \f$H(s)\f$. <BR>
                 Array size is `[ord+1 x 1]`.<BR>
                 Memory must be allocated.<BR><BR>
                 
\param[in]  ord  Filter order.<BR>
                 Number of  coefficients \f$H(s)\f$ and \f$H(z)\f$
                 numerator and denominator equals `ord+1`.<BR><BR>
                 
\param[out]  bz  Pointer to the numerator coefficients of a
                 digital filter transfer function \f$H(z)\f$. <BR>
                 Array size is `[ord+1 x 1]`.<BR>
                 Memory must be allocated.<BR><BR>

\param[out]  az  Pointer to the numerator coefficients of a
                 digital filter transfer function \f$H(z)\f$. <BR>
                 Array size is `[ord+1 x 1]`.<BR>
                 Memory must be allocated.<BR><BR>

\return
  `RES_OK`      if filter is calculated successfully.<BR><BR>
                Else \ref ERROR_CODE_GROUP "code error".<BR>

\author  Sergey Bakhurin www.dsplib.org          
*******************************************************************************/
int DSPL_API bilinear(double* bs, double* as, int ord, double* bz, double* az)
{
  double c[2] = {1.0, -1.0};
  double d[2] = {1.0,  1.0};
  return ratcompos(bs, as, ord, c, d, 1, bz, az);
}




/******************************************************************************
Digital IIR filter coefficients calculation
*******************************************************************************/
int DSPL_API iir(double rp, double rs, int ord,  double  w0, double  w1,
                                       int type, double* b,  double* a)
{
  double *bs = NULL;
  double *as = NULL;
  double *bt = NULL;
  double *at = NULL;
  double wa0, wa1;
  int err, ord_ap = ord;
  
  if(((type & DSPL_FILTER_TYPE_MASK) == DSPL_FILTER_LPF) ||
     ((type & DSPL_FILTER_TYPE_MASK) == DSPL_FILTER_HPF))
  {
    bs = (double*)malloc((ord_ap+1)*sizeof(double));
    as = (double*)malloc((ord_ap+1)*sizeof(double));
    bt = (double*)malloc((ord_ap+1)*sizeof(double));
    at = (double*)malloc((ord_ap+1)*sizeof(double));
  }
  
  
  if(((type & DSPL_FILTER_TYPE_MASK) == DSPL_FILTER_BPASS) ||
     ((type & DSPL_FILTER_TYPE_MASK) == DSPL_FILTER_BSTOP))
  {
    if(ord % 2)
      return ERROR_FILTER_ORD_BP;
    else
    {
      ord_ap = ord / 2;
      bs = (double*)malloc((ord_ap + 1)*sizeof(double));
      as = (double*)malloc((ord_ap + 1)*sizeof(double));
      bt = (double*)malloc((ord    + 1)*sizeof(double));
      at = (double*)malloc((ord    + 1)*sizeof(double));
    }
  }
  err = iir_ap(rp, rs, ord_ap, type, bs, as);
  if(err != RES_OK)
    goto error_proc;
  
  //frequency transformation 
  wa0 = tan(w0 * M_PI * 0.5);
  wa1 = tan(w1 * M_PI * 0.5);
  
  switch(type & DSPL_FILTER_TYPE_MASK)
  {
    
    case DSPL_FILTER_LPF:
      err = low2low(bs, as, ord_ap, 1.0, wa0, bt, at);
      break;
    
    case DSPL_FILTER_HPF:
      err = low2high(bs, as, ord_ap, 1.0, wa0, bt, at);
      break;
      
    case DSPL_FILTER_BPASS:
      err = low2bp(bs, as, ord_ap, 1.0, wa0, wa1, bt, at);
      break;
      
    case DSPL_FILTER_BSTOP:
      // need frequency transform ws ->  1  rad/s  
      err = low2bs(bs, as, ord_ap, 1.0, wa0, wa1, bt, at);
      break;
   
    default:
      err = ERROR_FILTER_TYPE;
      break;
  }
  if(err != RES_OK)
    goto error_proc;
  
  
  err = bilinear(bt, at, ord, b, a);


error_proc:

  if(bs)
    free(bs);
  if(as)
    free(as);
   if(bt)
    free(bt);
  if(at)
    free(at);
  
  return err;
  
}






/******************************************************************************
Analog prototype for IIR 
*******************************************************************************/
int iir_ap(double rp, double rs, int ord, int type, double* b, double* a)
{
  int err;
  switch(type & DSPL_FILTER_APPROX_MASK)
  {
    case DSPL_FILTER_BUTTER:
      err = butter_ap(rp, ord, b, a);
      break;
    case DSPL_FILTER_CHEBY1:
      err = cheby1_ap(rp, ord, b, a);
      break;
    case DSPL_FILTER_CHEBY2:
      err = cheby2_ap_wp1(rp, rs, ord, b, a);
      break;
    case DSPL_FILTER_ELLIP:
      err = ellip_ap(rp, rs, ord, b, a);
      break;
    default:
      err = ERROR_FILTER_APPROX; 
  }
  
  
  return err;
}




