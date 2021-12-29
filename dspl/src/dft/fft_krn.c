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
#include <float.h>

#include "dspl.h"
#include "dspl_internal.h"




#ifdef DOXYGEN_ENGLISH

#endif
#ifdef DOXYGEN_RUSSIAN

#endif
int fft_krn(complex_t* t0, complex_t* t1, fft_t* p, int n, int addr)
{
    int n1, n2, k, m, i;
    complex_t *pw = p->w+addr;
    complex_t tmp;
    
    n1 = 1;
    if(n % 4096 == 0) { n1 = 4096; goto label_size; }
    if(n % 2048 == 0) { n1 = 2048; goto label_size; }
    if(n % 1024 == 0) { n1 = 1024; goto label_size; }
    if(n %  512 == 0) { n1 =  512; goto label_size; }
    if(n %  256 == 0) { n1 =  256; goto label_size; }
    if(n %  128 == 0) { n1 =  128; goto label_size; }
    if(n %   64 == 0) { n1 =   64; goto label_size; }
    if(n %   32 == 0) { n1 =   32; goto label_size; }
    if(n %   16 == 0) { n1 =   16; goto label_size; }
    if(n %    7 == 0) { n1 =    7; goto label_size; }
    if(n %    8 == 0) { n1 =    8; goto label_size; }
    if(n %    5 == 0) { n1 =    5; goto label_size; }
    if(n %    4 == 0) { n1 =    4; goto label_size; }
    if(n %    3 == 0) { n1 =    3; goto label_size; }
    if(n %    2 == 0) { n1 =    2; goto label_size; }

label_size:
    if(n1 == 1)
    {
        for(k = 0; k < n; k++)
        {
            RE(t1[k]) = IM(t1[k]) = 0.0;
            for(m = 0; m < n; m++)
            {
                i = (k*m) % n;
                RE(tmp) = CMRE(t0[m], pw[i]);
                IM(tmp) = CMIM(t0[m], pw[i]);
                RE(t1[k]) += RE(tmp);
                IM(t1[k]) += IM(tmp);
            }
        }
    }
    else
    {
        n2 = n / n1;
        
        if(n2>1)
        {
            memcpy(t1, t0, n*sizeof(complex_t));
            matrix_transpose_cmplx(t1, n2, n1, t0);
        }
        
        if(n1 == 4096)
            for(k = 0; k < n2; k++)
                dft4096(t0+4096*k, t1+4096*k, p->w4096, p->w256);
              
        if(n1 == 2048)
            for(k = 0; k < n2; k++)
                dft2048(t0+2048*k, t1+2048*k, p->w2048, p->w32, p->w64);
        
        if(n1 == 1024)
            for(k = 0; k < n2; k++)
                dft1024(t0+1024*k, t1+1024*k, p->w1024, p->w32);
              
        if(n1 == 512)
            for(k = 0; k < n2; k++)
                dft512(t0+512*k, t1+512*k, p->w512, p->w32);
              
        if(n1 == 256)
            for(k = 0; k < n2; k++)
                dft256(t0+256*k, t1+256*k, p->w256);
              
        if(n1 == 128)
            for(k = 0; k < n2; k++)
                dft128(t0+128*k, t1+128*k, p->w128);
              
        if(n1 == 64)
            for(k = 0; k < n2; k++)
                dft64(t0+64*k, t1+64*k, p->w64);

        if(n1 == 32)
            for(k = 0; k < n2; k++)
                dft32(t0+32*k, t1+32*k, p->w32);
        
        if(n1 == 16)
            for(k = 0; k < n2; k++)
                dft16(t0+16*k, t1+16*k);
                
        if(n1 == 7)
            for(k = 0; k < n2; k++)
                dft7(t0+7*k, t1+7*k);
            
        if(n1 == 8)
            for(k = 0; k < n2; k++)
                dft8(t0+8*k, t1+8*k); 
                
        if(n1 == 5)
            for(k = 0; k < n2; k++)
                dft5(t0+5*k, t1+5*k);
     
        if(n1 == 4)
            for(k = 0; k < n2; k++)
                dft4(t0+4*k, t1+4*k);
        
        if(n1 == 3)
            for(k = 0; k < n2; k++)
                dft3(t0+3*k, t1+3*k);
        
        if(n1 == 2)
            for(k = 0; k < n2; k++)
                dft2(t0+2*k, t1+2*k);

        if(n2 > 1)
        {

            for(k =0; k < n; k++)
            {
                RE(t0[k]) = CMRE(t1[k], pw[k]);
                IM(t0[k]) = CMIM(t1[k], pw[k]);
            }

            matrix_transpose_cmplx(t0, n1, n2, t1);
            
            
            for(k = 0; k < n1; k++)
            {
                fft_krn(t1+k*n2, t0+k*n2, p, n2, addr+n);
            }
            
            matrix_transpose_cmplx(t0, n2, n1, t1);
        }
    }
    return RES_OK;
}


