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

#include <stdlib.h>
#include <stdio.h>
#include <string.h>

#include "dspl.h"
#include "dspl_internal.h"



/*******************************************************************************
2 points DFT
*******************************************************************************/
void dft2 (complex_t *x, complex_t* y)
{
    RE(y[0]) = RE(x[0]) + RE(x[1]);
    IM(y[0]) = IM(x[0]) + IM(x[1]);

    RE(y[1]) = RE(x[0]) - RE(x[1]);
    IM(y[1]) = IM(x[0]) - IM(x[1]);
}



/*******************************************************************************
3 points DFT (Winograd algorithm)
*******************************************************************************/
void dft3 (complex_t *x, complex_t* y)
{
    double a, b, c, d;

    a = RE(x[0]) - 0.5 * (RE(x[1]) + RE(x[2]));
    b = DFT3_W * (IM(x[1]) - IM(x[2]));
    c = IM(x[0]) - 0.5 * (IM(x[1]) + IM(x[2]));
    d = DFT3_W * (RE(x[2]) - RE(x[1]));

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
void dft4 (complex_t *x, complex_t* y)
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
void dft5 (complex_t *x, complex_t* y)
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



/*******************************************************************************
7 points DFT (Winograd algorithm)
*******************************************************************************/
void dft7 (complex_t *x, complex_t* y)
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

    /* Winograd paper mistake?! */
    RE(sum[5]) = RE(x[2]) + RE(x[5]);
    IM(sum[5]) = IM(x[2]) + IM(x[5]);

    /* Winograd paper mistake?! */
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
8 points DFT
*******************************************************************************/
void dft8 (complex_t *x, complex_t* y)
{
    complex_t t0[8];
    complex_t t1[8];
    double tmp;

    transpose2x4(x, t0);

    dft4(t0,   t1);
    dft4(t0+4, t1+4);

    /*    0.707106781186548 - 707106781186548i */
    tmp       = (RE(t1[5]) + IM(t1[5])) * DFT8_W;
    IM(t1[5]) = (IM(t1[5]) - RE(t1[5])) * DFT8_W;
    RE(t1[5]) = tmp;

     /*    0.000000000000000 - 1.000000000000000i */
    tmp       =  RE(t1[6]);
    RE(t1[6]) =  IM(t1[6]);
    IM(t1[6]) = -tmp;

    /* -0.707106781186548 - 707106781186548i */
    tmp       =  (IM(t1[7]) - RE(t1[7])) * DFT8_W;
    IM(t1[7]) = -(IM(t1[7]) + RE(t1[7])) * DFT8_W;
    RE(t1[7]) = tmp;

    transpose4x2(t1, t0);

    RE(t1[0]) = RE(t0[0]) + RE(t0[1]);
    IM(t1[0]) = IM(t0[0]) + IM(t0[1]);
    RE(t1[1]) = RE(t0[0]) - RE(t0[1]);
    IM(t1[1]) = IM(t0[0]) - IM(t0[1]);

    RE(t1[2]) = RE(t0[2]) + RE(t0[3]);
    IM(t1[2]) = IM(t0[2]) + IM(t0[3]);
    RE(t1[3]) = RE(t0[2]) - RE(t0[3]);
    IM(t1[3]) = IM(t0[2]) - IM(t0[3]);

    RE(t1[4]) = RE(t0[4]) + RE(t0[5]);
    IM(t1[4]) = IM(t0[4]) + IM(t0[5]);
    RE(t1[5]) = RE(t0[4]) - RE(t0[5]);
    IM(t1[5]) = IM(t0[4]) - IM(t0[5]);

    RE(t1[6]) = RE(t0[6]) + RE(t0[7]);
    IM(t1[6]) = IM(t0[6]) + IM(t0[7]);
    RE(t1[7]) = RE(t0[6]) - RE(t0[7]);
    IM(t1[7]) = IM(t0[6]) - IM(t0[7]);

    transpose2x4(t1, y);
}

/*******************************************************************************
16 points DFT (Winograd algorithm)
*******************************************************************************/
void dft16(complex_t *x, complex_t* y)
{
    complex_t t0[16];
    complex_t t1[16];
    double tmp;

    transpose4x4(x, t0);

    dft4(t0,    t1);
    dft4(t0+4,  t1+4);
    dft4(t0+8,  t1+8);
    dft4(t0+12, t1+12);

    /* #define DFT16_W1             0.923879532511287 */
    /* #define DFT16_W2             0.382683432365090 */
    /* #define DFT16_W3             0.707106781186548 */

    /* 0.923879532511287 - 0.382683432365090i */
    tmp       =  RE(t1[5]) * DFT16_W1 + IM(t1[5]) * DFT16_W2;
    IM(t1[5]) = -RE(t1[5]) * DFT16_W2 + IM(t1[5]) * DFT16_W1;
    RE(t1[5]) = tmp;

    /* 0.707106781186548 - 0.707106781186547i */
    tmp       = ( RE(t1[6]) + IM(t1[6])) * DFT16_W3;
    IM(t1[6]) = (-RE(t1[6]) + IM(t1[6])) * DFT16_W3;
    RE(t1[6]) = tmp;

    /* 0.382683432365090 - 0.923879532511287i */
    tmp       =  RE(t1[7]) * DFT16_W2 + IM(t1[7]) * DFT16_W1;
    IM(t1[7]) = -RE(t1[7]) * DFT16_W1 + IM(t1[7]) * DFT16_W2;
    RE(t1[7]) = tmp;

    /* 0.707106781186548 - 0.707106781186547i */
    tmp       = ( RE(t1[9]) + IM(t1[9])) * DFT16_W3;
    IM(t1[9]) = (-RE(t1[9]) + IM(t1[9])) * DFT16_W3;
    RE(t1[9]) = tmp;

    /* 0.000000000000000 - 1.000000000000000i */
    tmp        = RE(t1[10]);
    RE(t1[10]) = IM(t1[10]);
    IM(t1[10]) = -tmp;

    /* -0.707106781186547 - 0.707106781186548i */
    tmp        = (-RE(t1[11]) + IM(t1[11])) * DFT16_W3;
    IM(t1[11]) = (-RE(t1[11]) - IM(t1[11])) * DFT16_W3;
    RE(t1[11]) = tmp;

    /* 0.382683432365090 - 0.923879532511287i */
    tmp        =  RE(t1[13]) * DFT16_W2 + IM(t1[13]) * DFT16_W1;
    IM(t1[13]) = -RE(t1[13]) * DFT16_W1 + IM(t1[13]) * DFT16_W2;
    RE(t1[13]) = tmp;

    /* -0.707106781186547 - 0.707106781186548i */
    tmp        = (-RE(t1[14]) + IM(t1[14])) * DFT16_W3;
    IM(t1[14]) = (-RE(t1[14]) - IM(t1[14])) * DFT16_W3;
    RE(t1[14]) = tmp;


    /* -0.923879532511287 + 0.382683432365090i */
    tmp        = -RE(t1[15]) * DFT16_W1 - IM(t1[15]) * DFT16_W2;
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
32 points DFT (Winograd algorithm)
*******************************************************************************/
void dft32(complex_t *x, complex_t* y, complex_t* w)
{
    complex_t t0[32];
    complex_t t1[32];

    int i;

    transpose4x8(x, t0);
    dft8(t0,    t1);
    dft8(t0+8,  t1+8);
    dft8(t0+16, t1+16);
    dft8(t0+24, t1+24);
    

    for(i = 0; i < 32; i++)
    {
        RE(t0[i]) = CMRE(t1[i], w[i]); 
        IM(t0[i]) = CMIM(t1[i], w[i]); 
    }

    transpose8x4(t0, t1);

    for(i = 0; i < 8; i++)
        dft4(t1 + i*4,    t0 + i*4);

    transpose4x8(t0, y);
}


/*******************************************************************************
64 points DFT (Winograd algorithm)
*******************************************************************************/
void dft64(complex_t *x, complex_t* y, complex_t* w)
{
    complex_t t0[64];
    complex_t t1[64];

    int i;

    transpose8x8(x, t0);

    for(i = 0; i < 8; i++)
        dft8(t0 + i*8,    t1 + i*8);

    for(i = 0; i < 64; i++)
    {
        RE(t0[i]) = CMRE(t1[i], w[i]); 
        IM(t0[i]) = CMIM(t1[i], w[i]); 
    }
    
    transpose8x8(t0, t1);
    
    for(i = 0; i < 8; i++)
        dft8(t1 + i*8,    t0 + i*8);
    transpose8x8(t0, y);
}



/*******************************************************************************
256 points DFT (Winograd algorithm)
*******************************************************************************/
void dft128(complex_t *x, complex_t* y, complex_t* w)
{
    complex_t t0[128];
    complex_t t1[128];

    int i;

    matrix_transpose_cmplx(x,8,16,t0);

    for(i = 0; i < 8; i++)
        dft16(t0 + i*16,    t1 + i*16);

    for(i = 0; i < 128; i++)
    {
        RE(t0[i]) = CMRE(t1[i], w[i]); 
        IM(t0[i]) = CMIM(t1[i], w[i]); 
    }

    matrix_transpose_cmplx(t0, 16, 8, t1);

    for(i = 0; i < 16; i++)
        dft8(t1 + i*8,    t0 + i*8);

    matrix_transpose_cmplx(t0, 8, 16, y);
}



/*******************************************************************************
256 points DFT (Winograd algorithm)
*******************************************************************************/
void dft256(complex_t *x, complex_t* y, complex_t* w)
{
    complex_t t0[256];
    complex_t t1[256];

    int i;

    transpose16x16(x, t0);

    for(i = 0; i < 16; i++)
        dft16(t0 + i*16,    t1 + i*16);

    for(i = 0; i < 256; i++)
    {
        RE(t0[i]) = CMRE(t1[i], w[i]); 
        IM(t0[i]) = CMIM(t1[i], w[i]); 
    }

    transpose16x16(t0, t1);

    for(i = 0; i < 16; i++)
        dft16(t1 + i*16,    t0 + i*16);

    transpose16x16(t0, y);
}


/*******************************************************************************
512 points DFT (Winograd algorithm)
*******************************************************************************/
void dft512(complex_t *x, complex_t* y, complex_t* w, complex_t* w32)
{
    complex_t t0[512];
    complex_t t1[512];

    int i;

    matrix_transpose_cmplx(x,16,32,t0);

    for(i = 0; i < 16; i++)
        dft32(t0 + i*32,    t1 + i*32, w32);

    for(i = 0; i < 512; i++)
    {
        RE(t0[i]) = CMRE(t1[i], w[i]); 
        IM(t0[i]) = CMIM(t1[i], w[i]); 
    }

    matrix_transpose_cmplx(t0, 32, 16, t1);

    for(i = 0; i < 32; i++)
        dft16(t1 + i*16,    t0 + i*16);

    matrix_transpose_cmplx(t0, 16, 32, y);
}


/*******************************************************************************
1024 points DFT (Winograd algorithm)
*******************************************************************************/
void dft1024(complex_t *x, complex_t* y, complex_t* w, complex_t* w32)
{
    complex_t t0[1024];
    complex_t t1[1024];

    int i;

    matrix_transpose_cmplx(x,32,32,t0);

    for(i = 0; i < 32; i++)
        dft32(t0 + i*32,    t1 + i*32, w32);

    for(i = 0; i < 1024; i++)
    {
        RE(t0[i]) = CMRE(t1[i], w[i]); 
        IM(t0[i]) = CMIM(t1[i], w[i]); 
    }

    matrix_transpose_cmplx(t0, 32, 32, t1);

    for(i = 0; i < 32; i++)
        dft32(t1 + i*32,    t0 + i*32, w32);

    matrix_transpose_cmplx(t0, 32, 32, y);
}

/*******************************************************************************
2048 points DFT (Winograd algorithm)
*******************************************************************************/
void dft2048(complex_t *x, complex_t* y, complex_t* w, 
             complex_t* w32, complex_t* w64)
{
    complex_t *t0 = NULL;
    complex_t *t1 = NULL;

    int i;
    
    t0 = (complex_t*)malloc(2048*sizeof(complex_t));
    t1 = (complex_t*)malloc(2048*sizeof(complex_t));
    
    matrix_transpose_cmplx(x,32,64,t0);

    for(i = 0; i < 32; i++)
        dft64(t0 + i*64,    t1 + i*64, w64);

    for(i = 0; i < 2048; i++)
    {
        RE(t0[i]) = CMRE(t1[i], w[i]); 
        IM(t0[i]) = CMIM(t1[i], w[i]); 
    }

    matrix_transpose_cmplx(t0, 64, 32, t1);

    for(i = 0; i < 64; i++)
        dft32(t1 + i*32,    t0 + i*32, w32);

    matrix_transpose_cmplx(t0, 32, 64, y);
    
    free(t0);
    free(t1);
}



/*******************************************************************************
4096 points DFT (Winograd algorithm)
*******************************************************************************/
void dft4096(complex_t *x, complex_t* y, complex_t* w,  complex_t* w256)
{
    complex_t *t0 = NULL;
    complex_t *t1 = NULL;

    int i;
    
    t0 = (complex_t*)malloc(4096*sizeof(complex_t));
    t1 = (complex_t*)malloc(4096*sizeof(complex_t));
    
    matrix_transpose_cmplx(x,16,256,t0);

    for(i = 0; i < 16; i++)
        dft256(t0 + i*256,    t1 + i*256, w256);

    for(i = 0; i < 4096; i++)
    {
        RE(t0[i]) = CMRE(t1[i], w[i]); 
        IM(t0[i]) = CMIM(t1[i], w[i]); 
    }

    matrix_transpose_cmplx(t0, 256, 16, t1);

    for(i = 0; i < 256; i++)
        dft16(t1 + i*16,    t0 + i*16);

    matrix_transpose_cmplx(t0, 16, 256, y);
    
    free(t0);
    free(t1);
}



/*******************************************************************************
4 x 2 matrix transpose
*******************************************************************************/
void transpose4x2(complex_t *x, complex_t* y)
{
    RE(y[ 0]) = RE(x[ 0]);    IM(y[ 0]) = IM(x[ 0]);
    RE(y[ 1]) = RE(x[ 4]);    IM(y[ 1]) = IM(x[ 4]);
    RE(y[ 2]) = RE(x[ 1]);    IM(y[ 2]) = IM(x[ 1]);
    RE(y[ 3]) = RE(x[ 5]);    IM(y[ 3]) = IM(x[ 5]);
    RE(y[ 4]) = RE(x[ 2]);    IM(y[ 4]) = IM(x[ 2]);
    RE(y[ 5]) = RE(x[ 6]);    IM(y[ 5]) = IM(x[ 6]);
    RE(y[ 6]) = RE(x[ 3]);    IM(y[ 6]) = IM(x[ 3]);
    RE(y[ 7]) = RE(x[ 7]);    IM(y[ 7]) = IM(x[ 7]);
}


/*******************************************************************************
2 x 4 matrix transpose
*******************************************************************************/
void transpose2x4(complex_t *x, complex_t* y)
{
    RE(y[ 0]) = RE(x[ 0]);    IM(y[ 0]) = IM(x[ 0]);
    RE(y[ 1]) = RE(x[ 2]);    IM(y[ 1]) = IM(x[ 2]);
    RE(y[ 2]) = RE(x[ 4]);    IM(y[ 2]) = IM(x[ 4]);
    RE(y[ 3]) = RE(x[ 6]);    IM(y[ 3]) = IM(x[ 6]);
    RE(y[ 4]) = RE(x[ 1]);    IM(y[ 4]) = IM(x[ 1]);
    RE(y[ 5]) = RE(x[ 3]);    IM(y[ 5]) = IM(x[ 3]);
    RE(y[ 6]) = RE(x[ 5]);    IM(y[ 6]) = IM(x[ 5]);
    RE(y[ 7]) = RE(x[ 7]);    IM(y[ 7]) = IM(x[ 7]);
}



/*******************************************************************************
4 x 4 matrix transpose
*******************************************************************************/
void transpose4x4(complex_t *x, complex_t* y)
{
    RE(y[ 0]) = RE(x[ 0]);    IM(y[ 0]) = IM(x[ 0]);
    RE(y[ 1]) = RE(x[ 4]);    IM(y[ 1]) = IM(x[ 4]);
    RE(y[ 2]) = RE(x[ 8]);    IM(y[ 2]) = IM(x[ 8]);
    RE(y[ 3]) = RE(x[12]);    IM(y[ 3]) = IM(x[12]);
    RE(y[ 4]) = RE(x[ 1]);    IM(y[ 4]) = IM(x[ 1]);
    RE(y[ 5]) = RE(x[ 5]);    IM(y[ 5]) = IM(x[ 5]);
    RE(y[ 6]) = RE(x[ 9]);    IM(y[ 6]) = IM(x[ 9]);
    RE(y[ 7]) = RE(x[13]);    IM(y[ 7]) = IM(x[13]);
    RE(y[ 8]) = RE(x[ 2]);    IM(y[ 8]) = IM(x[ 2]);
    RE(y[ 9]) = RE(x[ 6]);    IM(y[ 9]) = IM(x[ 6]);
    RE(y[10]) = RE(x[10]);    IM(y[10]) = IM(x[10]);
    RE(y[11]) = RE(x[14]);    IM(y[11]) = IM(x[14]);
    RE(y[12]) = RE(x[ 3]);    IM(y[12]) = IM(x[ 3]);
    RE(y[13]) = RE(x[ 7]);    IM(y[13]) = IM(x[ 7]);
    RE(y[14]) = RE(x[11]);    IM(y[14]) = IM(x[11]);
    RE(y[15]) = RE(x[15]);    IM(y[15]) = IM(x[15]);
}


/*******************************************************************************
8 x 4 matrix transpose
*******************************************************************************/
void transpose8x4(complex_t *x, complex_t* y)
{
    RE(y[ 0]) = RE(x[  0]);    IM(y[ 0]) = IM(x[  0]);
    RE(y[ 1]) = RE(x[  8]);    IM(y[ 1]) = IM(x[  8]);
    RE(y[ 2]) = RE(x[ 16]);    IM(y[ 2]) = IM(x[ 16]);
    RE(y[ 3]) = RE(x[ 24]);    IM(y[ 3]) = IM(x[ 24]);
    RE(y[ 4]) = RE(x[  1]);    IM(y[ 4]) = IM(x[  1]);
    RE(y[ 5]) = RE(x[  9]);    IM(y[ 5]) = IM(x[  9]);
    RE(y[ 6]) = RE(x[ 17]);    IM(y[ 6]) = IM(x[ 17]);
    RE(y[ 7]) = RE(x[ 25]);    IM(y[ 7]) = IM(x[ 25]);
    RE(y[ 8]) = RE(x[  2]);    IM(y[ 8]) = IM(x[  2]);
    RE(y[ 9]) = RE(x[ 10]);    IM(y[ 9]) = IM(x[ 10]);
    RE(y[10]) = RE(x[ 18]);    IM(y[10]) = IM(x[ 18]);
    RE(y[11]) = RE(x[ 26]);    IM(y[11]) = IM(x[ 26]);
    RE(y[12]) = RE(x[  3]);    IM(y[12]) = IM(x[  3]);
    RE(y[13]) = RE(x[ 11]);    IM(y[13]) = IM(x[ 11]);
    RE(y[14]) = RE(x[ 19]);    IM(y[14]) = IM(x[ 19]);
    RE(y[15]) = RE(x[ 27]);    IM(y[15]) = IM(x[ 27]);
    RE(y[16]) = RE(x[  4]);    IM(y[16]) = IM(x[  4]);
    RE(y[17]) = RE(x[ 12]);    IM(y[17]) = IM(x[ 12]);
    RE(y[18]) = RE(x[ 20]);    IM(y[18]) = IM(x[ 20]);
    RE(y[19]) = RE(x[ 28]);    IM(y[19]) = IM(x[ 28]);
    RE(y[20]) = RE(x[  5]);    IM(y[20]) = IM(x[  5]);
    RE(y[21]) = RE(x[ 13]);    IM(y[21]) = IM(x[ 13]);
    RE(y[22]) = RE(x[ 21]);    IM(y[22]) = IM(x[ 21]);
    RE(y[23]) = RE(x[ 29]);    IM(y[23]) = IM(x[ 29]);
    RE(y[24]) = RE(x[  6]);    IM(y[24]) = IM(x[  6]);
    RE(y[25]) = RE(x[ 14]);    IM(y[25]) = IM(x[ 14]);
    RE(y[26]) = RE(x[ 22]);    IM(y[26]) = IM(x[ 22]);
    RE(y[27]) = RE(x[ 30]);    IM(y[27]) = IM(x[ 30]);
    RE(y[28]) = RE(x[  7]);    IM(y[28]) = IM(x[  7]);
    RE(y[29]) = RE(x[ 15]);    IM(y[29]) = IM(x[ 15]);
    RE(y[30]) = RE(x[ 23]);    IM(y[30]) = IM(x[ 23]);
    RE(y[31]) = RE(x[ 31]);    IM(y[31]) = IM(x[ 31]);
}


/*******************************************************************************
4 x 8 matrix transpose
*******************************************************************************/
void transpose4x8(complex_t *x, complex_t* y)
{
    RE(y[ 0]) = RE(x[  0]);    IM(y[ 0]) = IM(x[  0]);
    RE(y[ 1]) = RE(x[  4]);    IM(y[ 1]) = IM(x[  4]);
    RE(y[ 2]) = RE(x[  8]);    IM(y[ 2]) = IM(x[  8]);
    RE(y[ 3]) = RE(x[ 12]);    IM(y[ 3]) = IM(x[ 12]);
    RE(y[ 4]) = RE(x[ 16]);    IM(y[ 4]) = IM(x[ 16]);
    RE(y[ 5]) = RE(x[ 20]);    IM(y[ 5]) = IM(x[ 20]);
    RE(y[ 6]) = RE(x[ 24]);    IM(y[ 6]) = IM(x[ 24]);
    RE(y[ 7]) = RE(x[ 28]);    IM(y[ 7]) = IM(x[ 28]);
    RE(y[ 8]) = RE(x[  1]);    IM(y[ 8]) = IM(x[  1]);
    RE(y[ 9]) = RE(x[  5]);    IM(y[ 9]) = IM(x[  5]);
    RE(y[10]) = RE(x[  9]);    IM(y[10]) = IM(x[  9]);
    RE(y[11]) = RE(x[ 13]);    IM(y[11]) = IM(x[ 13]);
    RE(y[12]) = RE(x[ 17]);    IM(y[12]) = IM(x[ 17]);
    RE(y[13]) = RE(x[ 21]);    IM(y[13]) = IM(x[ 21]);
    RE(y[14]) = RE(x[ 25]);    IM(y[14]) = IM(x[ 25]);
    RE(y[15]) = RE(x[ 29]);    IM(y[15]) = IM(x[ 29]);
    RE(y[16]) = RE(x[  2]);    IM(y[16]) = IM(x[  2]);
    RE(y[17]) = RE(x[  6]);    IM(y[17]) = IM(x[  6]);
    RE(y[18]) = RE(x[ 10]);    IM(y[18]) = IM(x[ 10]);
    RE(y[19]) = RE(x[ 14]);    IM(y[19]) = IM(x[ 14]);
    RE(y[20]) = RE(x[ 18]);    IM(y[20]) = IM(x[ 18]);
    RE(y[21]) = RE(x[ 22]);    IM(y[21]) = IM(x[ 22]);
    RE(y[22]) = RE(x[ 26]);    IM(y[22]) = IM(x[ 26]);
    RE(y[23]) = RE(x[ 30]);    IM(y[23]) = IM(x[ 30]);
    RE(y[24]) = RE(x[  3]);    IM(y[24]) = IM(x[  3]);
    RE(y[25]) = RE(x[  7]);    IM(y[25]) = IM(x[  7]);
    RE(y[26]) = RE(x[ 11]);    IM(y[26]) = IM(x[ 11]);
    RE(y[27]) = RE(x[ 15]);    IM(y[27]) = IM(x[ 15]);
    RE(y[28]) = RE(x[ 19]);    IM(y[28]) = IM(x[ 19]);
    RE(y[29]) = RE(x[ 23]);    IM(y[29]) = IM(x[ 23]);
    RE(y[30]) = RE(x[ 27]);    IM(y[30]) = IM(x[ 27]);
    RE(y[31]) = RE(x[ 31]);    IM(y[31]) = IM(x[ 31]);
}


/*******************************************************************************
4 x 8 matrix transpose
*******************************************************************************/
void transpose8x8(complex_t *x, complex_t* y)
{
    RE(y[  0]) = RE(x[ 0]);    IM(y[  0]) = IM(x[ 0]);
    RE(y[  1]) = RE(x[ 8]);    IM(y[  1]) = IM(x[ 8]);
    RE(y[  2]) = RE(x[16]);    IM(y[  2]) = IM(x[16]);
    RE(y[  3]) = RE(x[24]);    IM(y[  3]) = IM(x[24]);
    RE(y[  4]) = RE(x[32]);    IM(y[  4]) = IM(x[32]);
    RE(y[  5]) = RE(x[40]);    IM(y[  5]) = IM(x[40]);
    RE(y[  6]) = RE(x[48]);    IM(y[  6]) = IM(x[48]);
    RE(y[  7]) = RE(x[56]);    IM(y[  7]) = IM(x[56]);
    RE(y[  8]) = RE(x[ 1]);    IM(y[  8]) = IM(x[ 1]);
    RE(y[  9]) = RE(x[ 9]);    IM(y[  9]) = IM(x[ 9]);
    RE(y[ 10]) = RE(x[17]);    IM(y[ 10]) = IM(x[17]);
    RE(y[ 11]) = RE(x[25]);    IM(y[ 11]) = IM(x[25]);
    RE(y[ 12]) = RE(x[33]);    IM(y[ 12]) = IM(x[33]);
    RE(y[ 13]) = RE(x[41]);    IM(y[ 13]) = IM(x[41]);
    RE(y[ 14]) = RE(x[49]);    IM(y[ 14]) = IM(x[49]);
    RE(y[ 15]) = RE(x[57]);    IM(y[ 15]) = IM(x[57]);
    RE(y[ 16]) = RE(x[ 2]);    IM(y[ 16]) = IM(x[ 2]);
    RE(y[ 17]) = RE(x[10]);    IM(y[ 17]) = IM(x[10]);
    RE(y[ 18]) = RE(x[18]);    IM(y[ 18]) = IM(x[18]);
    RE(y[ 19]) = RE(x[26]);    IM(y[ 19]) = IM(x[26]);
    RE(y[ 20]) = RE(x[34]);    IM(y[ 20]) = IM(x[34]);
    RE(y[ 21]) = RE(x[42]);    IM(y[ 21]) = IM(x[42]);
    RE(y[ 22]) = RE(x[50]);    IM(y[ 22]) = IM(x[50]);
    RE(y[ 23]) = RE(x[58]);    IM(y[ 23]) = IM(x[58]);
    RE(y[ 24]) = RE(x[ 3]);    IM(y[ 24]) = IM(x[ 3]);
    RE(y[ 25]) = RE(x[11]);    IM(y[ 25]) = IM(x[11]);
    RE(y[ 26]) = RE(x[19]);    IM(y[ 26]) = IM(x[19]);
    RE(y[ 27]) = RE(x[27]);    IM(y[ 27]) = IM(x[27]);
    RE(y[ 28]) = RE(x[35]);    IM(y[ 28]) = IM(x[35]);
    RE(y[ 29]) = RE(x[43]);    IM(y[ 29]) = IM(x[43]);
    RE(y[ 30]) = RE(x[51]);    IM(y[ 30]) = IM(x[51]);
    RE(y[ 31]) = RE(x[59]);    IM(y[ 31]) = IM(x[59]);
    RE(y[ 32]) = RE(x[ 4]);    IM(y[ 32]) = IM(x[ 4]);
    RE(y[ 33]) = RE(x[12]);    IM(y[ 33]) = IM(x[12]);
    RE(y[ 34]) = RE(x[20]);    IM(y[ 34]) = IM(x[20]);
    RE(y[ 35]) = RE(x[28]);    IM(y[ 35]) = IM(x[28]);
    RE(y[ 36]) = RE(x[36]);    IM(y[ 36]) = IM(x[36]);
    RE(y[ 37]) = RE(x[44]);    IM(y[ 37]) = IM(x[44]);
    RE(y[ 38]) = RE(x[52]);    IM(y[ 38]) = IM(x[52]);
    RE(y[ 39]) = RE(x[60]);    IM(y[ 39]) = IM(x[60]);
    RE(y[ 40]) = RE(x[ 5]);    IM(y[ 40]) = IM(x[ 5]);
    RE(y[ 41]) = RE(x[13]);    IM(y[ 41]) = IM(x[13]);
    RE(y[ 42]) = RE(x[21]);    IM(y[ 42]) = IM(x[21]);
    RE(y[ 43]) = RE(x[29]);    IM(y[ 43]) = IM(x[29]);
    RE(y[ 44]) = RE(x[37]);    IM(y[ 44]) = IM(x[37]);
    RE(y[ 45]) = RE(x[45]);    IM(y[ 45]) = IM(x[45]);
    RE(y[ 46]) = RE(x[53]);    IM(y[ 46]) = IM(x[53]);
    RE(y[ 47]) = RE(x[61]);    IM(y[ 47]) = IM(x[61]);
    RE(y[ 48]) = RE(x[ 6]);    IM(y[ 48]) = IM(x[ 6]);
    RE(y[ 49]) = RE(x[14]);    IM(y[ 49]) = IM(x[14]);
    RE(y[ 50]) = RE(x[22]);    IM(y[ 50]) = IM(x[22]);
    RE(y[ 51]) = RE(x[30]);    IM(y[ 51]) = IM(x[30]);
    RE(y[ 52]) = RE(x[38]);    IM(y[ 52]) = IM(x[38]);
    RE(y[ 53]) = RE(x[46]);    IM(y[ 53]) = IM(x[46]);
    RE(y[ 54]) = RE(x[54]);    IM(y[ 54]) = IM(x[54]);
    RE(y[ 55]) = RE(x[62]);    IM(y[ 55]) = IM(x[62]);
    RE(y[ 56]) = RE(x[ 7]);    IM(y[ 56]) = IM(x[ 7]);
    RE(y[ 57]) = RE(x[15]);    IM(y[ 57]) = IM(x[15]);
    RE(y[ 58]) = RE(x[23]);    IM(y[ 58]) = IM(x[23]);
    RE(y[ 59]) = RE(x[31]);    IM(y[ 59]) = IM(x[31]);
    RE(y[ 60]) = RE(x[39]);    IM(y[ 60]) = IM(x[39]);
    RE(y[ 61]) = RE(x[47]);    IM(y[ 61]) = IM(x[47]);
    RE(y[ 62]) = RE(x[55]);    IM(y[ 62]) = IM(x[55]);
    RE(y[ 63]) = RE(x[63]);    IM(y[ 63]) = IM(x[63]);
}



/*******************************************************************************
16 x 16 matrix transpose
*******************************************************************************/
void transpose16x16(complex_t* x, complex_t* y)
{
    RE(y[  0]) = RE(x[  0]);    IM(y[  0]) = IM(x[  0]);
    RE(y[  1]) = RE(x[ 16]);    IM(y[  1]) = IM(x[ 16]);
    RE(y[  2]) = RE(x[ 32]);    IM(y[  2]) = IM(x[ 32]);
    RE(y[  3]) = RE(x[ 48]);    IM(y[  3]) = IM(x[ 48]);
    RE(y[  4]) = RE(x[ 64]);    IM(y[  4]) = IM(x[ 64]);
    RE(y[  5]) = RE(x[ 80]);    IM(y[  5]) = IM(x[ 80]);
    RE(y[  6]) = RE(x[ 96]);    IM(y[  6]) = IM(x[ 96]);
    RE(y[  7]) = RE(x[112]);    IM(y[  7]) = IM(x[112]);
    RE(y[  8]) = RE(x[128]);    IM(y[  8]) = IM(x[128]);
    RE(y[  9]) = RE(x[144]);    IM(y[  9]) = IM(x[144]);
    RE(y[ 10]) = RE(x[160]);    IM(y[ 10]) = IM(x[160]);
    RE(y[ 11]) = RE(x[176]);    IM(y[ 11]) = IM(x[176]);
    RE(y[ 12]) = RE(x[192]);    IM(y[ 12]) = IM(x[192]);
    RE(y[ 13]) = RE(x[208]);    IM(y[ 13]) = IM(x[208]);
    RE(y[ 14]) = RE(x[224]);    IM(y[ 14]) = IM(x[224]);
    RE(y[ 15]) = RE(x[240]);    IM(y[ 15]) = IM(x[240]);
    RE(y[ 16]) = RE(x[  1]);    IM(y[ 16]) = IM(x[  1]);
    RE(y[ 17]) = RE(x[ 17]);    IM(y[ 17]) = IM(x[ 17]);
    RE(y[ 18]) = RE(x[ 33]);    IM(y[ 18]) = IM(x[ 33]);
    RE(y[ 19]) = RE(x[ 49]);    IM(y[ 19]) = IM(x[ 49]);
    RE(y[ 20]) = RE(x[ 65]);    IM(y[ 20]) = IM(x[ 65]);
    RE(y[ 21]) = RE(x[ 81]);    IM(y[ 21]) = IM(x[ 81]);
    RE(y[ 22]) = RE(x[ 97]);    IM(y[ 22]) = IM(x[ 97]);
    RE(y[ 23]) = RE(x[113]);    IM(y[ 23]) = IM(x[113]);
    RE(y[ 24]) = RE(x[129]);    IM(y[ 24]) = IM(x[129]);
    RE(y[ 25]) = RE(x[145]);    IM(y[ 25]) = IM(x[145]);
    RE(y[ 26]) = RE(x[161]);    IM(y[ 26]) = IM(x[161]);
    RE(y[ 27]) = RE(x[177]);    IM(y[ 27]) = IM(x[177]);
    RE(y[ 28]) = RE(x[193]);    IM(y[ 28]) = IM(x[193]);
    RE(y[ 29]) = RE(x[209]);    IM(y[ 29]) = IM(x[209]);
    RE(y[ 30]) = RE(x[225]);    IM(y[ 30]) = IM(x[225]);
    RE(y[ 31]) = RE(x[241]);    IM(y[ 31]) = IM(x[241]);
    RE(y[ 32]) = RE(x[  2]);    IM(y[ 32]) = IM(x[  2]);
    RE(y[ 33]) = RE(x[ 18]);    IM(y[ 33]) = IM(x[ 18]);
    RE(y[ 34]) = RE(x[ 34]);    IM(y[ 34]) = IM(x[ 34]);
    RE(y[ 35]) = RE(x[ 50]);    IM(y[ 35]) = IM(x[ 50]);
    RE(y[ 36]) = RE(x[ 66]);    IM(y[ 36]) = IM(x[ 66]);
    RE(y[ 37]) = RE(x[ 82]);    IM(y[ 37]) = IM(x[ 82]);
    RE(y[ 38]) = RE(x[ 98]);    IM(y[ 38]) = IM(x[ 98]);
    RE(y[ 39]) = RE(x[114]);    IM(y[ 39]) = IM(x[114]);
    RE(y[ 40]) = RE(x[130]);    IM(y[ 40]) = IM(x[130]);
    RE(y[ 41]) = RE(x[146]);    IM(y[ 41]) = IM(x[146]);
    RE(y[ 42]) = RE(x[162]);    IM(y[ 42]) = IM(x[162]);
    RE(y[ 43]) = RE(x[178]);    IM(y[ 43]) = IM(x[178]);
    RE(y[ 44]) = RE(x[194]);    IM(y[ 44]) = IM(x[194]);
    RE(y[ 45]) = RE(x[210]);    IM(y[ 45]) = IM(x[210]);
    RE(y[ 46]) = RE(x[226]);    IM(y[ 46]) = IM(x[226]);
    RE(y[ 47]) = RE(x[242]);    IM(y[ 47]) = IM(x[242]);
    RE(y[ 48]) = RE(x[  3]);    IM(y[ 48]) = IM(x[  3]);
    RE(y[ 49]) = RE(x[ 19]);    IM(y[ 49]) = IM(x[ 19]);
    RE(y[ 50]) = RE(x[ 35]);    IM(y[ 50]) = IM(x[ 35]);
    RE(y[ 51]) = RE(x[ 51]);    IM(y[ 51]) = IM(x[ 51]);
    RE(y[ 52]) = RE(x[ 67]);    IM(y[ 52]) = IM(x[ 67]);
    RE(y[ 53]) = RE(x[ 83]);    IM(y[ 53]) = IM(x[ 83]);
    RE(y[ 54]) = RE(x[ 99]);    IM(y[ 54]) = IM(x[ 99]);
    RE(y[ 55]) = RE(x[115]);    IM(y[ 55]) = IM(x[115]);
    RE(y[ 56]) = RE(x[131]);    IM(y[ 56]) = IM(x[131]);
    RE(y[ 57]) = RE(x[147]);    IM(y[ 57]) = IM(x[147]);
    RE(y[ 58]) = RE(x[163]);    IM(y[ 58]) = IM(x[163]);
    RE(y[ 59]) = RE(x[179]);    IM(y[ 59]) = IM(x[179]);
    RE(y[ 60]) = RE(x[195]);    IM(y[ 60]) = IM(x[195]);
    RE(y[ 61]) = RE(x[211]);    IM(y[ 61]) = IM(x[211]);
    RE(y[ 62]) = RE(x[227]);    IM(y[ 62]) = IM(x[227]);
    RE(y[ 63]) = RE(x[243]);    IM(y[ 63]) = IM(x[243]);
    RE(y[ 64]) = RE(x[  4]);    IM(y[ 64]) = IM(x[  4]);
    RE(y[ 65]) = RE(x[ 20]);    IM(y[ 65]) = IM(x[ 20]);
    RE(y[ 66]) = RE(x[ 36]);    IM(y[ 66]) = IM(x[ 36]);
    RE(y[ 67]) = RE(x[ 52]);    IM(y[ 67]) = IM(x[ 52]);
    RE(y[ 68]) = RE(x[ 68]);    IM(y[ 68]) = IM(x[ 68]);
    RE(y[ 69]) = RE(x[ 84]);    IM(y[ 69]) = IM(x[ 84]);
    RE(y[ 70]) = RE(x[100]);    IM(y[ 70]) = IM(x[100]);
    RE(y[ 71]) = RE(x[116]);    IM(y[ 71]) = IM(x[116]);
    RE(y[ 72]) = RE(x[132]);    IM(y[ 72]) = IM(x[132]);
    RE(y[ 73]) = RE(x[148]);    IM(y[ 73]) = IM(x[148]);
    RE(y[ 74]) = RE(x[164]);    IM(y[ 74]) = IM(x[164]);
    RE(y[ 75]) = RE(x[180]);    IM(y[ 75]) = IM(x[180]);
    RE(y[ 76]) = RE(x[196]);    IM(y[ 76]) = IM(x[196]);
    RE(y[ 77]) = RE(x[212]);    IM(y[ 77]) = IM(x[212]);
    RE(y[ 78]) = RE(x[228]);    IM(y[ 78]) = IM(x[228]);
    RE(y[ 79]) = RE(x[244]);    IM(y[ 79]) = IM(x[244]);
    RE(y[ 80]) = RE(x[  5]);    IM(y[ 80]) = IM(x[  5]);
    RE(y[ 81]) = RE(x[ 21]);    IM(y[ 81]) = IM(x[ 21]);
    RE(y[ 82]) = RE(x[ 37]);    IM(y[ 82]) = IM(x[ 37]);
    RE(y[ 83]) = RE(x[ 53]);    IM(y[ 83]) = IM(x[ 53]);
    RE(y[ 84]) = RE(x[ 69]);    IM(y[ 84]) = IM(x[ 69]);
    RE(y[ 85]) = RE(x[ 85]);    IM(y[ 85]) = IM(x[ 85]);
    RE(y[ 86]) = RE(x[101]);    IM(y[ 86]) = IM(x[101]);
    RE(y[ 87]) = RE(x[117]);    IM(y[ 87]) = IM(x[117]);
    RE(y[ 88]) = RE(x[133]);    IM(y[ 88]) = IM(x[133]);
    RE(y[ 89]) = RE(x[149]);    IM(y[ 89]) = IM(x[149]);
    RE(y[ 90]) = RE(x[165]);    IM(y[ 90]) = IM(x[165]);
    RE(y[ 91]) = RE(x[181]);    IM(y[ 91]) = IM(x[181]);
    RE(y[ 92]) = RE(x[197]);    IM(y[ 92]) = IM(x[197]);
    RE(y[ 93]) = RE(x[213]);    IM(y[ 93]) = IM(x[213]);
    RE(y[ 94]) = RE(x[229]);    IM(y[ 94]) = IM(x[229]);
    RE(y[ 95]) = RE(x[245]);    IM(y[ 95]) = IM(x[245]);
    RE(y[ 96]) = RE(x[  6]);    IM(y[ 96]) = IM(x[  6]);
    RE(y[ 97]) = RE(x[ 22]);    IM(y[ 97]) = IM(x[ 22]);
    RE(y[ 98]) = RE(x[ 38]);    IM(y[ 98]) = IM(x[ 38]);
    RE(y[ 99]) = RE(x[ 54]);    IM(y[ 99]) = IM(x[ 54]);
    RE(y[100]) = RE(x[ 70]);    IM(y[100]) = IM(x[ 70]);
    RE(y[101]) = RE(x[ 86]);    IM(y[101]) = IM(x[ 86]);
    RE(y[102]) = RE(x[102]);    IM(y[102]) = IM(x[102]);
    RE(y[103]) = RE(x[118]);    IM(y[103]) = IM(x[118]);
    RE(y[104]) = RE(x[134]);    IM(y[104]) = IM(x[134]);
    RE(y[105]) = RE(x[150]);    IM(y[105]) = IM(x[150]);
    RE(y[106]) = RE(x[166]);    IM(y[106]) = IM(x[166]);
    RE(y[107]) = RE(x[182]);    IM(y[107]) = IM(x[182]);
    RE(y[108]) = RE(x[198]);    IM(y[108]) = IM(x[198]);
    RE(y[109]) = RE(x[214]);    IM(y[109]) = IM(x[214]);
    RE(y[110]) = RE(x[230]);    IM(y[110]) = IM(x[230]);
    RE(y[111]) = RE(x[246]);    IM(y[111]) = IM(x[246]);
    RE(y[112]) = RE(x[  7]);    IM(y[112]) = IM(x[  7]);
    RE(y[113]) = RE(x[ 23]);    IM(y[113]) = IM(x[ 23]);
    RE(y[114]) = RE(x[ 39]);    IM(y[114]) = IM(x[ 39]);
    RE(y[115]) = RE(x[ 55]);    IM(y[115]) = IM(x[ 55]);
    RE(y[116]) = RE(x[ 71]);    IM(y[116]) = IM(x[ 71]);
    RE(y[117]) = RE(x[ 87]);    IM(y[117]) = IM(x[ 87]);
    RE(y[118]) = RE(x[103]);    IM(y[118]) = IM(x[103]);
    RE(y[119]) = RE(x[119]);    IM(y[119]) = IM(x[119]);
    RE(y[120]) = RE(x[135]);    IM(y[120]) = IM(x[135]);
    RE(y[121]) = RE(x[151]);    IM(y[121]) = IM(x[151]);
    RE(y[122]) = RE(x[167]);    IM(y[122]) = IM(x[167]);
    RE(y[123]) = RE(x[183]);    IM(y[123]) = IM(x[183]);
    RE(y[124]) = RE(x[199]);    IM(y[124]) = IM(x[199]);
    RE(y[125]) = RE(x[215]);    IM(y[125]) = IM(x[215]);
    RE(y[126]) = RE(x[231]);    IM(y[126]) = IM(x[231]);
    RE(y[127]) = RE(x[247]);    IM(y[127]) = IM(x[247]);
    RE(y[128]) = RE(x[  8]);    IM(y[128]) = IM(x[  8]);
    RE(y[129]) = RE(x[ 24]);    IM(y[129]) = IM(x[ 24]);
    RE(y[130]) = RE(x[ 40]);    IM(y[130]) = IM(x[ 40]);
    RE(y[131]) = RE(x[ 56]);    IM(y[131]) = IM(x[ 56]);
    RE(y[132]) = RE(x[ 72]);    IM(y[132]) = IM(x[ 72]);
    RE(y[133]) = RE(x[ 88]);    IM(y[133]) = IM(x[ 88]);
    RE(y[134]) = RE(x[104]);    IM(y[134]) = IM(x[104]);
    RE(y[135]) = RE(x[120]);    IM(y[135]) = IM(x[120]);
    RE(y[136]) = RE(x[136]);    IM(y[136]) = IM(x[136]);
    RE(y[137]) = RE(x[152]);    IM(y[137]) = IM(x[152]);
    RE(y[138]) = RE(x[168]);    IM(y[138]) = IM(x[168]);
    RE(y[139]) = RE(x[184]);    IM(y[139]) = IM(x[184]);
    RE(y[140]) = RE(x[200]);    IM(y[140]) = IM(x[200]);
    RE(y[141]) = RE(x[216]);    IM(y[141]) = IM(x[216]);
    RE(y[142]) = RE(x[232]);    IM(y[142]) = IM(x[232]);
    RE(y[143]) = RE(x[248]);    IM(y[143]) = IM(x[248]);
    RE(y[144]) = RE(x[  9]);    IM(y[144]) = IM(x[  9]);
    RE(y[145]) = RE(x[ 25]);    IM(y[145]) = IM(x[ 25]);
    RE(y[146]) = RE(x[ 41]);    IM(y[146]) = IM(x[ 41]);
    RE(y[147]) = RE(x[ 57]);    IM(y[147]) = IM(x[ 57]);
    RE(y[148]) = RE(x[ 73]);    IM(y[148]) = IM(x[ 73]);
    RE(y[149]) = RE(x[ 89]);    IM(y[149]) = IM(x[ 89]);
    RE(y[150]) = RE(x[105]);    IM(y[150]) = IM(x[105]);
    RE(y[151]) = RE(x[121]);    IM(y[151]) = IM(x[121]);
    RE(y[152]) = RE(x[137]);    IM(y[152]) = IM(x[137]);
    RE(y[153]) = RE(x[153]);    IM(y[153]) = IM(x[153]);
    RE(y[154]) = RE(x[169]);    IM(y[154]) = IM(x[169]);
    RE(y[155]) = RE(x[185]);    IM(y[155]) = IM(x[185]);
    RE(y[156]) = RE(x[201]);    IM(y[156]) = IM(x[201]);
    RE(y[157]) = RE(x[217]);    IM(y[157]) = IM(x[217]);
    RE(y[158]) = RE(x[233]);    IM(y[158]) = IM(x[233]);
    RE(y[159]) = RE(x[249]);    IM(y[159]) = IM(x[249]);
    RE(y[160]) = RE(x[ 10]);    IM(y[160]) = IM(x[ 10]);
    RE(y[161]) = RE(x[ 26]);    IM(y[161]) = IM(x[ 26]);
    RE(y[162]) = RE(x[ 42]);    IM(y[162]) = IM(x[ 42]);
    RE(y[163]) = RE(x[ 58]);    IM(y[163]) = IM(x[ 58]);
    RE(y[164]) = RE(x[ 74]);    IM(y[164]) = IM(x[ 74]);
    RE(y[165]) = RE(x[ 90]);    IM(y[165]) = IM(x[ 90]);
    RE(y[166]) = RE(x[106]);    IM(y[166]) = IM(x[106]);
    RE(y[167]) = RE(x[122]);    IM(y[167]) = IM(x[122]);
    RE(y[168]) = RE(x[138]);    IM(y[168]) = IM(x[138]);
    RE(y[169]) = RE(x[154]);    IM(y[169]) = IM(x[154]);
    RE(y[170]) = RE(x[170]);    IM(y[170]) = IM(x[170]);
    RE(y[171]) = RE(x[186]);    IM(y[171]) = IM(x[186]);
    RE(y[172]) = RE(x[202]);    IM(y[172]) = IM(x[202]);
    RE(y[173]) = RE(x[218]);    IM(y[173]) = IM(x[218]);
    RE(y[174]) = RE(x[234]);    IM(y[174]) = IM(x[234]);
    RE(y[175]) = RE(x[250]);    IM(y[175]) = IM(x[250]);
    RE(y[176]) = RE(x[ 11]);    IM(y[176]) = IM(x[ 11]);
    RE(y[177]) = RE(x[ 27]);    IM(y[177]) = IM(x[ 27]);
    RE(y[178]) = RE(x[ 43]);    IM(y[178]) = IM(x[ 43]);
    RE(y[179]) = RE(x[ 59]);    IM(y[179]) = IM(x[ 59]);
    RE(y[180]) = RE(x[ 75]);    IM(y[180]) = IM(x[ 75]);
    RE(y[181]) = RE(x[ 91]);    IM(y[181]) = IM(x[ 91]);
    RE(y[182]) = RE(x[107]);    IM(y[182]) = IM(x[107]);
    RE(y[183]) = RE(x[123]);    IM(y[183]) = IM(x[123]);
    RE(y[184]) = RE(x[139]);    IM(y[184]) = IM(x[139]);
    RE(y[185]) = RE(x[155]);    IM(y[185]) = IM(x[155]);
    RE(y[186]) = RE(x[171]);    IM(y[186]) = IM(x[171]);
    RE(y[187]) = RE(x[187]);    IM(y[187]) = IM(x[187]);
    RE(y[188]) = RE(x[203]);    IM(y[188]) = IM(x[203]);
    RE(y[189]) = RE(x[219]);    IM(y[189]) = IM(x[219]);
    RE(y[190]) = RE(x[235]);    IM(y[190]) = IM(x[235]);
    RE(y[191]) = RE(x[251]);    IM(y[191]) = IM(x[251]);
    RE(y[192]) = RE(x[ 12]);    IM(y[192]) = IM(x[ 12]);
    RE(y[193]) = RE(x[ 28]);    IM(y[193]) = IM(x[ 28]);
    RE(y[194]) = RE(x[ 44]);    IM(y[194]) = IM(x[ 44]);
    RE(y[195]) = RE(x[ 60]);    IM(y[195]) = IM(x[ 60]);
    RE(y[196]) = RE(x[ 76]);    IM(y[196]) = IM(x[ 76]);
    RE(y[197]) = RE(x[ 92]);    IM(y[197]) = IM(x[ 92]);
    RE(y[198]) = RE(x[108]);    IM(y[198]) = IM(x[108]);
    RE(y[199]) = RE(x[124]);    IM(y[199]) = IM(x[124]);
    RE(y[200]) = RE(x[140]);    IM(y[200]) = IM(x[140]);
    RE(y[201]) = RE(x[156]);    IM(y[201]) = IM(x[156]);
    RE(y[202]) = RE(x[172]);    IM(y[202]) = IM(x[172]);
    RE(y[203]) = RE(x[188]);    IM(y[203]) = IM(x[188]);
    RE(y[204]) = RE(x[204]);    IM(y[204]) = IM(x[204]);
    RE(y[205]) = RE(x[220]);    IM(y[205]) = IM(x[220]);
    RE(y[206]) = RE(x[236]);    IM(y[206]) = IM(x[236]);
    RE(y[207]) = RE(x[252]);    IM(y[207]) = IM(x[252]);
    RE(y[208]) = RE(x[ 13]);    IM(y[208]) = IM(x[ 13]);
    RE(y[209]) = RE(x[ 29]);    IM(y[209]) = IM(x[ 29]);
    RE(y[210]) = RE(x[ 45]);    IM(y[210]) = IM(x[ 45]);
    RE(y[211]) = RE(x[ 61]);    IM(y[211]) = IM(x[ 61]);
    RE(y[212]) = RE(x[ 77]);    IM(y[212]) = IM(x[ 77]);
    RE(y[213]) = RE(x[ 93]);    IM(y[213]) = IM(x[ 93]);
    RE(y[214]) = RE(x[109]);    IM(y[214]) = IM(x[109]);
    RE(y[215]) = RE(x[125]);    IM(y[215]) = IM(x[125]);
    RE(y[216]) = RE(x[141]);    IM(y[216]) = IM(x[141]);
    RE(y[217]) = RE(x[157]);    IM(y[217]) = IM(x[157]);
    RE(y[218]) = RE(x[173]);    IM(y[218]) = IM(x[173]);
    RE(y[219]) = RE(x[189]);    IM(y[219]) = IM(x[189]);
    RE(y[220]) = RE(x[205]);    IM(y[220]) = IM(x[205]);
    RE(y[221]) = RE(x[221]);    IM(y[221]) = IM(x[221]);
    RE(y[222]) = RE(x[237]);    IM(y[222]) = IM(x[237]);
    RE(y[223]) = RE(x[253]);    IM(y[223]) = IM(x[253]);
    RE(y[224]) = RE(x[ 14]);    IM(y[224]) = IM(x[ 14]);
    RE(y[225]) = RE(x[ 30]);    IM(y[225]) = IM(x[ 30]);
    RE(y[226]) = RE(x[ 46]);    IM(y[226]) = IM(x[ 46]);
    RE(y[227]) = RE(x[ 62]);    IM(y[227]) = IM(x[ 62]);
    RE(y[228]) = RE(x[ 78]);    IM(y[228]) = IM(x[ 78]);
    RE(y[229]) = RE(x[ 94]);    IM(y[229]) = IM(x[ 94]);
    RE(y[230]) = RE(x[110]);    IM(y[230]) = IM(x[110]);
    RE(y[231]) = RE(x[126]);    IM(y[231]) = IM(x[126]);
    RE(y[232]) = RE(x[142]);    IM(y[232]) = IM(x[142]);
    RE(y[233]) = RE(x[158]);    IM(y[233]) = IM(x[158]);
    RE(y[234]) = RE(x[174]);    IM(y[234]) = IM(x[174]);
    RE(y[235]) = RE(x[190]);    IM(y[235]) = IM(x[190]);
    RE(y[236]) = RE(x[206]);    IM(y[236]) = IM(x[206]);
    RE(y[237]) = RE(x[222]);    IM(y[237]) = IM(x[222]);
    RE(y[238]) = RE(x[238]);    IM(y[238]) = IM(x[238]);
    RE(y[239]) = RE(x[254]);    IM(y[239]) = IM(x[254]);
    RE(y[240]) = RE(x[ 15]);    IM(y[240]) = IM(x[ 15]);
    RE(y[241]) = RE(x[ 31]);    IM(y[241]) = IM(x[ 31]);
    RE(y[242]) = RE(x[ 47]);    IM(y[242]) = IM(x[ 47]);
    RE(y[243]) = RE(x[ 63]);    IM(y[243]) = IM(x[ 63]);
    RE(y[244]) = RE(x[ 79]);    IM(y[244]) = IM(x[ 79]);
    RE(y[245]) = RE(x[ 95]);    IM(y[245]) = IM(x[ 95]);
    RE(y[246]) = RE(x[111]);    IM(y[246]) = IM(x[111]);
    RE(y[247]) = RE(x[127]);    IM(y[247]) = IM(x[127]);
    RE(y[248]) = RE(x[143]);    IM(y[248]) = IM(x[143]);
    RE(y[249]) = RE(x[159]);    IM(y[249]) = IM(x[159]);
    RE(y[250]) = RE(x[175]);    IM(y[250]) = IM(x[175]);
    RE(y[251]) = RE(x[191]);    IM(y[251]) = IM(x[191]);
    RE(y[252]) = RE(x[207]);    IM(y[252]) = IM(x[207]);
    RE(y[253]) = RE(x[223]);    IM(y[253]) = IM(x[223]);
    RE(y[254]) = RE(x[239]);    IM(y[254]) = IM(x[239]);
    RE(y[255]) = RE(x[255]);    IM(y[255]) = IM(x[255]);  
}


