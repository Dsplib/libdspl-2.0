#ifndef RANDOMGEN_H
#define RANDOMGEN_H


#include "dspl.h"


#define DSPL_RAND_MOD_X1            2147483647
#define DSPL_RAND_MOD_X2            2145483479

/* random MRG32K3A algorithm constants */
#define MRG32K3A_NORM       2.328306549295728e-10
#define MRG32K3A_M1         4294967087.0
#define MRG32K3A_M2         4294944443.0
#define MRG32K3A_A12           1403580.0
#define MRG32K3A_A13            810728.0
#define MRG32K3A_A21            527612.0
#define MRG32K3A_A23           1370589.0
#define RAND_BUFSIZE        512


int randu_mrg32k3a (double* u, int n, random_t* prnd);

/* 
   A C-program for MT19937-64 (2004/9/29 version).
   Coded by Takuji Nishimura and Makoto Matsumoto.

   This is a 64-bit version of Mersenne Twister pseudorandom number
   generator.

   Before using, initialize the state by using init_genrand64(seed)  
   or init_by_array64(init_key, key_length).

   Copyright (C) 2004, Makoto Matsumoto and Takuji Nishimura,
   All rights reserved.                          

   Redistribution and use in source and binary forms, with or without
   modification, are permitted provided that the following conditions
   are met:

     1. Redistributions of source code must retain the above copyright
        notice, this list of conditions and the following disclaimer.

     2. Redistributions in binary form must reproduce the above copyright
        notice, this list of conditions and the following disclaimer in the
        documentation and/or other materials provided with the distribution.

     3. The names of its contributors may not be used to endorse or promote 
        products derived from this software without specific prior written 
        permission.

   THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
   "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
   LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR
   A PARTICULAR PURPOSE ARE DISCLAIMED.  IN NO EVENT SHALL THE COPYRIGHT OWNER OR
   CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
   EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
   PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR
   PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF
   LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
   NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
   SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

   References:
   T. Nishimura, ``Tables of 64-bit Mersenne Twisters''
     ACM Transactions on Modeling and 
     Computer Simulation 10. (2000) 348--357.
   M. Matsumoto and T. Nishimura,
     ``Mersenne Twister: a 623-dimensionally equidistributed
       uniform pseudorandom number generator''
     ACM Transactions on Modeling and 
     Computer Simulation 8. (Jan. 1998) 3--30.

   Any feedback is very welcome.
   http://www.math.hiroshima-u.ac.jp/~m-mat/MT/emt.html
   email: m-mat @ math.sci.hiroshima-u.ac.jp (remove spaces)
*/


/* initializes mt[NN] with a seed */
void mt19937_init_genrand64(unsigned long long seed, random_t* prnd);

/* initialize by an array with array-length */
/* init_key is the array for initializing keys */
/* key_length is its length */
void mt19937_init_by_array64(unsigned long long init_key[], 
                             unsigned long long key_length,
                             random_t* prnd);

/* generates a random number on [0, 2^64-1]-interval */
unsigned long long mt19937_genrand64_int64(random_t* prnd);

/* generates a random number on [0, 2^63-1]-interval */
long long mt19937_genrand64_int63(random_t* prnd);

/* generates a random number on [0,1]-real-interval */
double mt19937_genrand64_real1(random_t* prnd);

/* generates a random number on [0,1)-real-interval */
double mt19937_genrand64_real2(random_t* prnd);

/* generates a random number on (0,1)-real-interval */
double mt19937_genrand64_real3(random_t* prnd);



#endif