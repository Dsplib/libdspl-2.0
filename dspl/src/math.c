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
module operator for double
*******************************************************************************/
double DSPL_API dmod (double x, double y)
{
  if(y == 0.0)
    return x;
  return  x - floor(x/y) * y;
}





/*******************************************************************************
sinc(x) = sin(pi*x)/(pi*x)
*******************************************************************************/
int DSPL_API sinc(double* x, int n, double a, double* y)
{
  int k;

  if(!x || !y)
    return ERROR_PTR;
  if(n<1)
    return ERROR_SIZE;

  for(k = 0; k < n; k++)
    y[k] = (x[k]==0.0) ? 1.0 : sin(a*x[k])/(a*x[k]);

  return RES_OK;


}



/*******************************************************************************
Sine integral
--------------------------------------------------------------------------------
This function uses Padé approximants of the convergent Taylor series [1]


[1] 
https://www.sciencedirect.com/science/article/pii/S221313371500013X?via%3Dihub

*******************************************************************************/
int DSPL_API sine_int(double* x, int n, double* si)
{
  int k, sgn, p;
  double num, den, y, x2, x22, z, f, g;
  
  double A[8] =   {+1.00000000000000000E0   ,
                   -4.54393409816329991E-2  ,
                   +1.15457225751016682E-3  ,
                   -1.41018536821330254E-5  ,
                   +9.43280809438713025E-8  ,
                   -3.53201978997168357E-10 ,
                   +7.08240282274875911E-13 ,
                   -6.05338212010422477E-16 };
                 
  
  
  double B[7]  =  {+1.0                     ,
                   +1.01162145739225565E-2  ,
                   +4.99175116169755106E-5  ,
                   +1.55654986308745614E-7  ,
                   +3.28067571055789734E-10 ,
                   +4.50490975753865810E-13 ,
                   +3.21107051193712168E-16 };
  
  
                 
  double FA[11] = {+1.000000000000000000000E0,
                   +7.444370681619367006180E2, 
                   +1.963963728951468698010E5, 
                   +2.377503101254318340340E7, 
                   +1.430734038212746368880E9, 
                   +4.33736238870432522765E10, 
                   +6.40533830574022022911E11, 
                   +4.20968180571076940208E12, 
                   +1.00795182980368574617E13, 
                   +4.94816688199951963482E12, 
                   -4.94701168645415959931E11}; 
   
  double FB[10] = {+1.000000000000000000000E0,
                   +7.464370681619276780310E2, 
                   +1.978652470315839514500E5, 
                   +2.415356701651268451440E7, 
                   +1.474789521929854649580E9, 
                   +4.58595115847765779830E10, 
                   +7.08501308149515401563E11, 
                   +5.06084464593475076774E12, 
                   +1.43468549171581016479E13, 
                   +1.11535493509914254097E13};  
                   
                   
    
  double GA[11] = {+1.000000000000000000E0,
                   +8.135952011516861500E2, 
                   +2.352391816264782000E5, 
                   +3.125575707957787310E7, 
                   +2.062975951467633540E9, 
                   +6.83052205423625007E10, 
                   +1.09049528450362786E12, 
                   +7.57664583257834349E12, 
                   +1.81004487464664575E13, 
                   +6.43291613143049485E12, 
                   -1.36517137670871689E12}; 
    
    
  double GB[10] = {+1.000000000000000000E0,
                   +8.195952011514515640E2, 
                   +2.400367528355787770E5, 
                   +3.260266616470908220E7, 
                   +2.233555432780993600E9, 
                   +7.87465017341829930E10, 
                   +1.39866710696414565E12, 
                   +1.17164723371736605E13, 
                   +4.01839087307656620E13, 
                   +3.99653257887490811E13};
                 
  if(!x || !si)
    return ERROR_PTR;
  if(n<1)
    return ERROR_SIZE;
  
  
  for(p = 0; p < n; p++)
  {    
    sgn = x[p] > 0.0 ?  0 : 1;
    y   = x[p] < 0.0 ? -x[p] : x[p];  
    
    if(y < 4)
    {
      x2 = y * y;
      z = 1.0;
      num = 0.0;
      for(k = 0; k < 8; k++)
      {
        num += A[k] * z;
        z*=x2;       
      }
      z = 1.0;
      den = 0.0;
      for(k = 0; k < 7; k++)
      {
        den += B[k]*z;
        z*=x2;
      }      
      si[p] = x[p] * num/den; 
    }
    else
    {
      
      x2 = 1.0/y;
      x22 = x2*x2;
      z = 1.0;
      num = 0.0;
      for(k = 0; k < 11; k++)
      {      
        num += FA[k] * z;
        z*=x22;       
      }
      z = 1.0;
      den = 0.0;
      for(k = 0; k < 10; k++)
      {
        den += FB[k]*z;
        z*=x22;
      }   
    
      f = x2 * num / den;     
      
      z = 1.0;
      num = 0.0;
      for(k = 0; k < 11; k++)
      {      
        num += GA[k] * z;
        z*=x22;       
      }
      z = 1.0;
      den = 0.0;
      for(k = 0; k < 10; k++)
      {
        den += GB[k]*z;
        z*=x22;
      }   
    
      g = x22 * num / den;
      
      si[p] = sgn ? f * cos(y) + g * sin(y) - M_PI * 0.5 :
                    M_PI * 0.5 - f * cos(y) - g * sin(y);
    }  
  }
  return RES_OK; 
}
