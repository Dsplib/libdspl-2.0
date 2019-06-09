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
2 points DFT
*******************************************************************************/
void dft2 (complex_t *x,  complex_t* y)
{
  RE(y[0]) = RE(x[0]) + RE(x[1]);
  IM(y[0]) = IM(x[0]) + IM(x[1]);

  RE(y[1]) = RE(x[0]) - RE(x[1]);
  IM(y[1]) = IM(x[0]) - IM(x[1]);
}



/*******************************************************************************
3 points DFT (Winograd algorithm)
*******************************************************************************/
void dft3 (complex_t *x,  complex_t* y)
{
  double a, b, c, d;
  
  a = RE(x[0]) - 0.5 * (RE(x[1]) + RE(x[2]));
  b =         DFT3_W * (IM(x[1]) - IM(x[2]));
  c = IM(x[0]) - 0.5 * (IM(x[1]) + IM(x[2]));
  d =         DFT3_W * (RE(x[2]) - RE(x[1]));
   
  RE(y[0]) = RE(x[0]) + RE(x[1]) + RE(x[2]);
  IM(y[0]) = IM(x[0]) + IM(x[1]) + IM(x[2]);
  
  RE(y[1]) = a + b;
  IM(y[1]) = c + d;

  RE(y[2]) = a - b;
  IM(y[2]) = c - d;
}



/*******************************************************************************
4 points DFT (Winograd algorithm)
*******************************************************************************/
void dft4 (complex_t *x,  complex_t* y)
{
  double a0, b0, c0, d0, a1, b1, c1, d1;
  
  a0 = RE(x[0]) + RE(x[2]);
  b0 = RE(x[1]) + RE(x[3]);
  c0 = IM(x[0]) + IM(x[2]);
  d0 = IM(x[1]) + IM(x[3]);
  
  a1 = RE(x[0]) - RE(x[2]);
  b1 = RE(x[1]) - RE(x[3]);
  c1 = IM(x[0]) - IM(x[2]);
  d1 = IM(x[1]) - IM(x[3]);
  
  RE(y[0]) = a0 + b0;
  IM(y[0]) = c0 + d0;

  RE(y[1]) = a1 + d1;
  IM(y[1]) = c1 - b1;
  
  RE(y[2]) = a0 - b0;
  IM(y[2]) = c0 - d0;

  RE(y[3]) = a1 - d1;
  IM(y[3]) = c1 + b1;  
}




/*******************************************************************************
5 points DFT (Winograd algorithm)
*******************************************************************************/
void dft5 (complex_t *x,  complex_t* y)
{
  complex_t sum[14];
  complex_t mul[6];

  RE(sum[1]) = RE(x[1]) + RE(x[4]);
  IM(sum[1]) = IM(x[1]) + IM(x[4]);
  
  RE(sum[2]) = RE(x[1]) - RE(x[4]);
  IM(sum[2]) = IM(x[1]) - IM(x[4]);
  
  RE(sum[3]) = RE(x[3]) + RE(x[2]);
  IM(sum[3]) = IM(x[3]) + IM(x[2]);
  
  RE(sum[4]) = RE(x[3]) - RE(x[2]);
  IM(sum[4]) = IM(x[3]) - IM(x[2]);  
  
  RE(sum[5]) = RE(sum[1]) + RE(sum[3]);
  IM(sum[5]) = IM(sum[1]) + IM(sum[3]);  
  
  RE(sum[6]) = RE(sum[1]) - RE(sum[3]);
  IM(sum[6]) = IM(sum[1]) - IM(sum[3]);  
  
  RE(sum[7]) = RE(sum[2]) + RE(sum[4]);
  IM(sum[7]) = IM(sum[2]) + IM(sum[4]);  
  
  RE(y[0]) = RE(sum[5]) + RE(x[0]);
  IM(y[0]) = IM(sum[5]) + IM(x[0]);
  
  RE(mul[1]) = RE(sum[5]) * DFT5_W1;
  IM(mul[1]) = IM(sum[5]) * DFT5_W1;
  
  RE(mul[2]) = RE(sum[6]) * DFT5_W2;
  IM(mul[2]) = IM(sum[6]) * DFT5_W2;  
  
  RE(mul[3]) = -IM(sum[2]) * DFT5_W3;
  IM(mul[3]) =  RE(sum[2]) * DFT5_W3;
  
  RE(mul[4]) = -IM(sum[7]) * DFT5_W4;
  IM(mul[4]) =  RE(sum[7]) * DFT5_W4;
  
  RE(mul[5]) = -IM(sum[4]) * DFT5_W5;
  IM(mul[5]) =  RE(sum[4]) * DFT5_W5;
  
  RE(sum[9]) = RE(y[0]) + RE(mul[1]);
  IM(sum[9]) = IM(y[0]) + IM(mul[1]);
  
  RE(sum[10]) = RE(sum[9]) + RE(mul[2]);
  IM(sum[10]) = IM(sum[9]) + IM(mul[2]);
  
  RE(sum[11]) = RE(sum[9]) - RE(mul[2]);
  IM(sum[11]) = IM(sum[9]) - IM(mul[2]);
  
  RE(sum[12]) = RE(mul[3]) - RE(mul[4]);
  IM(sum[12]) = IM(mul[3]) - IM(mul[4]);
  
  RE(sum[13]) = RE(mul[4]) + RE(mul[5]);
  IM(sum[13]) = IM(mul[4]) + IM(mul[5]);
  
  RE(y[4]) = RE(sum[10]) + RE(sum[12]);
  IM(y[4]) = IM(sum[10]) + IM(sum[12]);
  
  RE(y[3]) = RE(sum[11]) + RE(sum[13]);
  IM(y[3]) = IM(sum[11]) + IM(sum[13]);
  
  RE(y[2]) = RE(sum[11]) - RE(sum[13]);
  IM(y[2]) = IM(sum[11]) - IM(sum[13]);  
  
  RE(y[1]) = RE(sum[10]) - RE(sum[12]);
  IM(y[1]) = IM(sum[10]) - IM(sum[12]);
}


void dft7 (complex_t *x,  complex_t* y)
{
  complex_t sum[31];
  complex_t mul[9];
  
  RE(sum[1]) = RE(x[1]) + RE(x[6]);
  IM(sum[1]) = IM(x[1]) + IM(x[6]);
  
  RE(sum[2]) = RE(x[1]) - RE(x[6]);
  IM(sum[2]) = IM(x[1]) - IM(x[6]);
  
  RE(sum[3]) = RE(x[4]) + RE(x[3]);
  IM(sum[3]) = IM(x[4]) + IM(x[3]);
  
  RE(sum[4]) = RE(x[4]) - RE(x[3]);
  IM(sum[4]) = IM(x[4]) - IM(x[3]);
  
  // Winograd paper mistake?!
  RE(sum[5]) = RE(x[2]) + RE(x[5]);
  IM(sum[5]) = IM(x[2]) + IM(x[5]);
  // Winograd paper mistake?!
  RE(sum[6]) = RE(x[2]) - RE(x[5]);
  IM(sum[6]) = IM(x[2]) - IM(x[5]);
  
  RE(sum[7]) = RE(sum[1]) + RE(sum[3]);
  IM(sum[7]) = IM(sum[1]) + IM(sum[3]);
  
  RE(sum[8]) = RE(sum[7]) + RE(sum[5]);
  IM(sum[8]) = IM(sum[7]) + IM(sum[5]);
  
  RE(y[0]) = RE(sum[9]) = RE(sum[8]) + RE(x[0]);
  IM(y[0]) = IM(sum[9]) = IM(sum[8]) + IM(x[0]);
  
  RE(sum[10]) = RE(sum[1]) - RE(sum[3]);
  IM(sum[10]) = IM(sum[1]) - IM(sum[3]);
  
  RE(sum[11]) = RE(sum[3]) - RE(sum[5]);
  IM(sum[11]) = IM(sum[3]) - IM(sum[5]);

  RE(sum[12]) = RE(sum[5]) - RE(sum[1]);
  IM(sum[12]) = IM(sum[5]) - IM(sum[1]);
    
  RE(sum[13]) = RE(sum[2]) + RE(sum[4]);
  IM(sum[13]) = IM(sum[2]) + IM(sum[4]);
  
  RE(sum[14]) = RE(sum[13]) + RE(sum[6]);
  IM(sum[14]) = IM(sum[13]) + IM(sum[6]);
  
  RE(sum[15]) = RE(sum[2]) - RE(sum[4]);
  IM(sum[15]) = IM(sum[2]) - IM(sum[4]);
  
  RE(sum[16]) = RE(sum[4]) - RE(sum[6]);
  IM(sum[16]) = IM(sum[4]) - IM(sum[6]);
  
  RE(sum[17]) = RE(sum[6]) - RE(sum[2]);
  IM(sum[17]) = IM(sum[6]) - IM(sum[2]);
  
  RE(mul[1]) = DFT7_W1 * RE(sum[8]); 
  IM(mul[1]) = DFT7_W1 * IM(sum[8]);
  
  RE(mul[2]) = DFT7_W2 * RE(sum[10]); 
  IM(mul[2]) = DFT7_W2 * IM(sum[10]);
  
  RE(mul[3]) = DFT7_W3 * RE(sum[11]); 
  IM(mul[3]) = DFT7_W3 * IM(sum[11]);
  
  RE(mul[4]) = DFT7_W4 * RE(sum[12]); 
  IM(mul[4]) = DFT7_W4 * IM(sum[12]);
  
  RE(mul[5]) = -DFT7_W5 * IM(sum[14]); 
  IM(mul[5]) =  DFT7_W5 * RE(sum[14]);
  
  RE(mul[6]) = -DFT7_W6 * IM(sum[15]); 
  IM(mul[6]) =  DFT7_W6 * RE(sum[15]);
  
  RE(mul[7]) = -DFT7_W7 * IM(sum[16]); 
  IM(mul[7]) =  DFT7_W7 * RE(sum[16]);
  
  RE(mul[8]) = -DFT7_W8 * IM(sum[17]); 
  IM(mul[8]) =  DFT7_W8 * RE(sum[17]);
  
  RE(sum[18]) = RE(y[0]) + RE(mul[1]);
  IM(sum[18]) = IM(y[0]) + IM(mul[1]);
  
  RE(sum[19]) = RE(sum[18]) + RE(mul[2]);
  IM(sum[19]) = IM(sum[18]) + IM(mul[2]);
  
  RE(sum[20]) = RE(sum[19]) + RE(mul[3]);
  IM(sum[20]) = IM(sum[19]) + IM(mul[3]);
  
  RE(sum[21]) = RE(sum[18]) - RE(mul[2]);
  IM(sum[21]) = IM(sum[18]) - IM(mul[2]);
  
  RE(sum[22]) = RE(sum[21]) - RE(mul[4]);
  IM(sum[22]) = IM(sum[21]) - IM(mul[4]);
  
  RE(sum[23]) = RE(sum[18]) - RE(mul[3]);
  IM(sum[23]) = IM(sum[18]) - IM(mul[3]);
  
  RE(sum[24]) = RE(sum[23]) + RE(mul[4]);
  IM(sum[24]) = IM(sum[23]) + IM(mul[4]);
  
  RE(sum[25]) = RE(mul[5]) + RE(mul[6]);
  IM(sum[25]) = IM(mul[5]) + IM(mul[6]);
  
  RE(sum[26]) = RE(sum[25]) + RE(mul[7]);
  IM(sum[26]) = IM(sum[25]) + IM(mul[7]);
  
  RE(sum[27]) = RE(mul[5]) - RE(mul[6]);
  IM(sum[27]) = IM(mul[5]) - IM(mul[6]);
  
  RE(sum[28]) = RE(sum[27]) - RE(mul[8]);
  IM(sum[28]) = IM(sum[27]) - IM(mul[8]);
  
  RE(sum[29]) = RE(mul[5]) - RE(mul[7]);
  IM(sum[29]) = IM(mul[5]) - IM(mul[7]);
  
  RE(sum[30]) = RE(sum[29]) + RE(mul[8]);
  IM(sum[30]) = IM(sum[29]) + IM(mul[8]);
  
  RE(y[1]) = RE(sum[20]) + RE(sum[26]);
  IM(y[1]) = IM(sum[20]) + IM(sum[26]);
  
  RE(y[2]) = RE(sum[22]) + RE(sum[28]);
  IM(y[2]) = IM(sum[22]) + IM(sum[28]);
  
  RE(y[3]) = RE(sum[24]) - RE(sum[30]);
  IM(y[3]) = IM(sum[24]) - IM(sum[30]);
  
  RE(y[4]) = RE(sum[24]) + RE(sum[30]);
  IM(y[4]) = IM(sum[24]) + IM(sum[30]);
  
  RE(y[5]) = RE(sum[22]) - RE(sum[28]);
  IM(y[5]) = IM(sum[22]) - IM(sum[28]);  
  
  RE(y[6]) = RE(sum[20]) - RE(sum[26]);
  IM(y[6]) = IM(sum[20]) - IM(sum[26]);
  
}


/*******************************************************************************
16 points DFT (Winograd algorithm)
*******************************************************************************/
void dft16 (complex_t *x,  complex_t* y)
{
  complex_t t0[16];
  complex_t t1[16];
  double tmp;
  
  transpose4x4(x, t0);
  
  dft4(t0,    t1);
  dft4(t0+4,  t1+4);
  dft4(t0+8,  t1+8);
  dft4(t0+12, t1+12);
  
  //#define DFT16_W1       0.923879532511287
  //#define DFT16_W2       0.382683432365090
  //#define DFT16_W3       0.707106781186548
  
  // 0.923879532511287 - 0.382683432365090i
  tmp       =  RE(t1[5]) * DFT16_W1 + IM(t1[5]) * DFT16_W2;
  IM(t1[5]) = -RE(t1[5]) * DFT16_W2 + IM(t1[5]) * DFT16_W1;
  RE(t1[5]) = tmp;
  
  // 0.707106781186548 - 0.707106781186547i
  tmp       =  ( RE(t1[6]) + IM(t1[6])) * DFT16_W3;
  IM(t1[6]) =  (-RE(t1[6]) + IM(t1[6])) * DFT16_W3;
  RE(t1[6]) = tmp;
    
  // 0.382683432365090 - 0.923879532511287i
  tmp       =  RE(t1[7]) * DFT16_W2 + IM(t1[7]) * DFT16_W1;
  IM(t1[7]) = -RE(t1[7]) * DFT16_W1 + IM(t1[7]) * DFT16_W2;
  RE(t1[7]) = tmp;
  
  // 0.707106781186548 - 0.707106781186547i
  tmp       =  ( RE(t1[9]) + IM(t1[9])) * DFT16_W3;
  IM(t1[9]) =  (-RE(t1[9]) + IM(t1[9])) * DFT16_W3;
  RE(t1[9]) = tmp;
  
  // 0.000000000000000 - 1.000000000000000i
  tmp       =   RE(t1[10]);
  RE(t1[10]) =  IM(t1[10]);
  IM(t1[10]) = -tmp;
  
  //-0.707106781186547 - 0.707106781186548i
  tmp        =  (-RE(t1[11]) + IM(t1[11])) * DFT16_W3;
  IM(t1[11]) =  (-RE(t1[11]) - IM(t1[11])) * DFT16_W3;
  RE(t1[11]) = tmp;
  
  // 0.382683432365090 - 0.923879532511287i
  tmp       =  RE(t1[13]) * DFT16_W2 + IM(t1[13]) * DFT16_W1;
  IM(t1[13]) = -RE(t1[13]) * DFT16_W1 + IM(t1[13]) * DFT16_W2;
  RE(t1[13]) = tmp;
  
  //-0.707106781186547 - 0.707106781186548i
  tmp        =  (-RE(t1[14]) + IM(t1[14])) * DFT16_W3;
  IM(t1[14]) =  (-RE(t1[14]) - IM(t1[14])) * DFT16_W3;
  RE(t1[14]) = tmp;
  
  
  //-0.923879532511287 + 0.382683432365090i
  tmp       =  -RE(t1[15]) * DFT16_W1 - IM(t1[15]) * DFT16_W2;
  IM(t1[15]) =  RE(t1[15]) * DFT16_W2 - IM(t1[15]) * DFT16_W1;
  RE(t1[15]) = tmp;
  
  transpose4x4(t1, t0);
  
  dft4(t0,    t1);
  dft4(t0+4,  t1+4);
  dft4(t0+8,  t1+8);
  dft4(t0+12, t1+12);
  
  transpose4x4(t1, y);

}

/*******************************************************************************
4 x 4 matrix transpose
*******************************************************************************/
void transpose4x4(complex_t *x,  complex_t* y)
{
  RE(y[ 0]) = RE(x[ 0]);  IM(y[ 0]) = IM(x[ 0]);   
  RE(y[ 1]) = RE(x[ 4]);  IM(y[ 1]) = IM(x[ 4]);   
  RE(y[ 2]) = RE(x[ 8]);  IM(y[ 2]) = IM(x[ 8]);   
  RE(y[ 3]) = RE(x[12]);  IM(y[ 3]) = IM(x[12]);   
  RE(y[ 4]) = RE(x[ 1]);  IM(y[ 4]) = IM(x[ 1]);   
  RE(y[ 5]) = RE(x[ 5]);  IM(y[ 5]) = IM(x[ 5]);   
  RE(y[ 6]) = RE(x[ 9]);  IM(y[ 6]) = IM(x[ 9]);   
  RE(y[ 7]) = RE(x[13]);  IM(y[ 7]) = IM(x[13]);   
  RE(y[ 8]) = RE(x[ 2]);  IM(y[ 8]) = IM(x[ 2]);  
  RE(y[ 9]) = RE(x[ 6]);  IM(y[ 9]) = IM(x[ 6]); 
  RE(y[10]) = RE(x[10]);  IM(y[10]) = IM(x[10]); 
  RE(y[11]) = RE(x[14]);  IM(y[11]) = IM(x[14]); 
  RE(y[12]) = RE(x[ 3]);  IM(y[12]) = IM(x[ 3]); 
  RE(y[13]) = RE(x[ 7]);  IM(y[13]) = IM(x[ 7]); 
  RE(y[14]) = RE(x[11]);  IM(y[14]) = IM(x[11]); 
  RE(y[15]) = RE(x[15]);  IM(y[15]) = IM(x[15]); 
}


