/*
* Copyright (c) 2015-2018 Sergey Bakhurin
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


#ifndef DSPL_H
#define DSPL_H



#include <math.h>
  

/* math const definition */
#ifndef M_PI
	#define M_PI	3.1415926535897932384626433832795
#endif
#ifndef M_2PI
	#define M_2PI	6.283185307179586476925286766559
#endif
     


typedef double complex_t[2];	 



typedef struct
{
    void      *pfftw;
    complex_t *in;
    complex_t *out;
    size_t     size;
} fft_t;


#define RE(x)           (x[0])
#define IM(x)           (x[1]) 


#define SQR(x)          ((x) * (x)) 
#define ABSSQR(x)       ((SQR(RE(x))) + (SQR(IM(x))))
#define ABS(x)          sqrt((ABSSQR(x)))

#define CMRE(a,b)       ((RE(a)) * (RE(b)) - (IM(a)) * (IM(b)))
#define CMIM(a,b)       ((RE(a)) * (IM(b)) + (IM(a)) * (RE(b)))

#define CMCONJRE(a, b)  ((RE(a)) * (RE(b)) + (IM(a)) * (IM(b)))
#define CMCONJIM(a, b)  ((IM(a)) * (RE(b)) - (RE(a)) * (IM(b))) 





#define RES_OK						    0

/* Error codes 									  */
/* A									0x01xxxxxx*/
/* B									0x02xxxxxx*/
/* C									0x03xxxxxx*/
/* D									0x04xxxxxx*/
#define ERROR_DAT_TYPE                  0x04012020
#define ERROR_DIV_ZERO                  0x04102226
/* E									0x05xxxxxx*/
#define ERROR_ELLIP_MODULE              0x05121315
/* F									0x06xxxxxx*/
#define ERROR_FFT_CREATE                0x06060318
#define ERROR_FILTER_A0                 0x06100100
#define ERROR_FILTER_ORD                0x06101518
#define ERROR_FILTER_RP                 0x06101816
#define ERROR_FILTER_RS                 0x06101819
#define ERROR_FNAME                     0x06140113
#define ERROR_FOPEN                     0x06151605
#define ERROR_FWRITE_SIZE               0x06231820
/* G									0x07xxxxxx*/
/* H									0x08xxxxxx*/
/* I									0x09xxxxxx*/
#define ERROR_INC_SIZE                  0x09140319
/* J									0x10xxxxxx*/
/* K									0x11xxxxxx*/
/* L									0x12xxxxxx*/
/* M									0x13xxxxxx*/
/* N									0x14xxxxxx*/
#define ERROR_NEGATIVE                  0x14050701
/* O									0x15xxxxxx*/
/* P									0x16xxxxxx*/
#define ERROR_POLY_ORD                  0x16151518
#define ERROR_PTR					    0x16201800
/* Q									0x17xxxxxx*/
/* R									0x18xxxxxx*/
#define ERROR_RAND_SIGMA                0x18011909
#define ERROR_RESAMPLE_RATIO            0x18051801
#define ERROR_RESAMPLE_FRAC_DELAY       0x18050604
/* S									0x19xxxxxx*/
#define ERROR_SIZE					    0x19092605
#define ERROR_SYM_TYPE                  0x19251320
/* T									0x20xxxxxx*/
/* U									0x21xxxxxx*/
/* V									0x22xxxxxx*/
/* W									0x23xxxxxx*/
#define ERROR_WIN_TYPE                  0x23092025
#define ERROR_WIN_SYM                   0x23091925
/* X									0x24xxxxxx*/
/* Y									0x25xxxxxx*/
/* Z									0x26xxxxxx*/


#define DAT_DOUBLE      0
#define DAT_COMPLEX     1




#define DSPL_SYMMETRIC                  0x00000000
#define DSPL_PERIODIC                   0x00000001

#define DSPL_WIN_SYM_MASK               0x00000001
#define DSPL_WIN_MASK                   0x000FFFFE

#define DSPL_WIN_SYMMETRIC              DSPL_SYMMETRIC
#define DSPL_WIN_PERIODIC               DSPL_PERIODIC


#define DSPL_WIN_BARTLETT               0x00000004
#define DSPL_WIN_BARTLETT_HANN          0x00000008	
#define DSPL_WIN_BLACKMAN               0x00000010			
#define DSPL_WIN_BLACKMAN_HARRIS        0x00000040	 
#define DSPL_WIN_BLACKMAN_NUTTALL       0x00000080    	
#define DSPL_WIN_FLAT_TOP               0x00000100			
#define DSPL_WIN_GAUSSIAN               0x00000400			
#define DSPL_WIN_HAMMING                0x00000800			
#define DSPL_WIN_HANN                   0x00001000				
#define DSPL_WIN_LANCZOS                0x00004000			
#define DSPL_WIN_NUTTALL                0x00008000			
#define DSPL_WIN_RECT                   0x00010000				
#define DSPL_WIN_COS                    0x00040000


#define ELLIP_ITER       16
#define ELLIP_MAX_ORD    24



#ifdef __cplusplus
extern "C" {
#endif


// Declare DSPL_API for Windows OS
#ifdef BUILD_LIB

#ifdef WIN_OS
#define DSPL_API __declspec(dllexport)
#endif // WIN_OS

#ifdef LINUX_OS
#define DSPL_API
#endif //LINUX_OS


int DSPL_API blas_dscal(int n, double a, double* x, int incx);

int DSPL_API cheby_poly1(double* x, int n, int ord, double* y);
int DSPL_API cheby_poly2(double* x, int n, int ord, double* y);

int	DSPL_API	dft			(double* x,		int n,	complex_t *y);
int	DSPL_API	dft_cmplx	(complex_t* x,	int n,	complex_t *y);


int DSPL_API	fft_create	( fft_t *pfft,	int n);
void DSPL_API	fft_destroy	(fft_t *pfft);
int DSPL_API	fft_cmplx	(complex_t *x,	int n, fft_t* pfft, complex_t* y);
int DSPL_API	fft_shift	(double* x,		int n, double* y);


#else //BUILD_LIB 



typedef int (*p_blas_dscal  )   (int n, double a, double* x, int incx);

typedef int (*p_cheby_poly1 )   (double* x, int n, int ord, double* y);
typedef int (*p_cheby_poly2 )   (double* x, int n, int ord, double* y);

typedef int (*p_dft			)	(double* x,		int n,	complex_t *y);
typedef int (*p_dft_cmplx	)	(complex_t* x,	int n,	complex_t *y);
typedef int (*p_fft_create	)	(fft_t *pfft,	int n);
typedef void(*p_fft_destroy	)	(fft_t *pfft);
typedef int (*p_fft_cmplx	)	(complex_t *x,	int n, fft_t* pfft, complex_t* y);
typedef int (*p_fft_shift	)	(double* x,		int n, double* y);



extern p_blas_dscal         blas_dscal      ;
      
extern p_cheby_poly1        cheby_poly1     ;
extern p_cheby_poly2        cheby_poly2     ;

extern p_dft				dft				;
extern p_dft_cmplx			dft_cmplx		;

extern p_fft_create			fft_create		;	
extern p_fft_destroy	    fft_destroy		;
extern p_fft_cmplx	        fft_cmplx		;
extern p_fft_shift	        fft_shift		;

#endif //BUILD_DLL


#ifdef __cplusplus
}
#endif


void*	dspl_load();
void	dspl_free(void* handle);



#endif //DSPL_H

