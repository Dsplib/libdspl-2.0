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



#ifdef WIN_OS
#include <windows.h>
#endif  //WIN_OS

#ifdef LINUX_OS
#include <dlfcn.h>
#endif //LINUX_OS


#include <stdio.h>
#include "dspl.h"


#ifndef BUILD_LIB



p_acos_cmplx                acos_cmplx          ;
p_asin_cmplx                asin_cmplx          ;
p_cheby_poly1               cheby_poly1		    ;
p_cheby_poly2               cheby_poly2		    ;
p_cmplx2re                  cmplx2re            ;
p_concat                    concat              ;
p_conv                      conv                ;
p_conv_cmplx                conv_cmplx          ;
p_cos_cmplx                 cos_cmplx           ;
p_dft                      	dft				    ;
p_dft_cmplx               	dft_cmplx		    ;
p_dmod                      dmod                ;
p_farrow_lagrange           farrow_lagrange     ;
p_farrow_spline             farrow_spline       ;
p_filter_iir                filter_iir          ;
p_fft                       fft                 ;
p_fft_cmplx                 fft_cmplx           ;   
p_fft_create                fft_create          ; 
p_fft_free                  fft_free            ; 
p_fft_shift                 fft_shift           ;
p_flipip                    flipip              ;
p_flipip_cmplx              flipip_cmplx        ;
p_fourier_series_dec        fourier_series_dec  ;
p_fourier_series_rec        fourier_series_rec  ;
p_goertzel                  goertzel            ;
p_goertzel_cmplx            goertzel_cmplx      ;
p_linspace                  linspace            ;
p_log_cmplx                 log_cmplx           ;
p_logspace                  logspace            ;
p_polyval                   polyval             ;
p_polyval_cmplx             polyval_cmplx       ;
p_randn                     randn               ;
p_randu                     randu               ;
p_re2cmplx                  re2cmplx            ;
p_signal_pimp               signal_pimp         ;
p_sin_cmplx                 sin_cmplx           ;
p_sqrt_cmplx                sqrt_cmplx          ;
p_trapint                   trapint             ;
p_trapint_cmplx             trapint_cmplx       ;
p_writebin                  writebin            ;
p_writetxt                  writetxt            ;

#endif //BUILD_LIB




#ifdef WIN_OS
#define LOAD_FUNC(fn) \
        fname = #fn;\
        fn = (p_##fn)GetProcAddress(handle, fname);\
        if(! fn) goto exit_label;
#endif



#ifdef LINUX_OS
#define LOAD_FUNC(fn) \
        fname = #fn;\
        fn = (p_##fn)dlsym(handle, fname);\
        if ((error = dlerror()) != NULL) goto exit_label
#endif




void* dspl_load()
{
    char* fname;
	#ifdef WIN_OS
		HINSTANCE handle;		
		handle = LoadLibrary(TEXT("libdspl.dll"));
		if (!handle)
		{
			printf("libdspl.dll loading ERROR!\n");
			return NULL;
		}
    #endif //WIN_OS		


    #ifdef LINUX_OS
		char* error;
		void *handle;
		// open the *.so
		handle = dlopen ("./libdspl.so", RTLD_LAZY);
		if (!handle) 
		{
			printf("libdspl.so loading ERROR!\n");
			return NULL;
		}
    #endif	//LINUX_OS

    LOAD_FUNC(acos_cmplx);
    LOAD_FUNC(asin_cmplx);    
    LOAD_FUNC(cheby_poly1);
    LOAD_FUNC(cheby_poly2);
    LOAD_FUNC(cmplx2re);
    LOAD_FUNC(concat);
    LOAD_FUNC(conv);
    LOAD_FUNC(conv_cmplx);
    LOAD_FUNC(cos_cmplx);
    LOAD_FUNC(dft);
    LOAD_FUNC(dft_cmplx);
    LOAD_FUNC(dmod);
    LOAD_FUNC(farrow_lagrange);
    LOAD_FUNC(farrow_spline);
    LOAD_FUNC(filter_iir);
    LOAD_FUNC(fft);
    LOAD_FUNC(fft_cmplx);
    LOAD_FUNC(fft_create);
    LOAD_FUNC(fft_free);
    LOAD_FUNC(fft_shift);
    LOAD_FUNC(flipip);
    LOAD_FUNC(flipip_cmplx);
    LOAD_FUNC(fourier_series_dec);
    LOAD_FUNC(fourier_series_rec);
    LOAD_FUNC(goertzel);
    LOAD_FUNC(goertzel_cmplx);
    LOAD_FUNC(linspace);
    LOAD_FUNC(log_cmplx);
    LOAD_FUNC(logspace);
    LOAD_FUNC(polyval);
    LOAD_FUNC(polyval_cmplx);
    LOAD_FUNC(randn);
    LOAD_FUNC(randu);
    LOAD_FUNC(re2cmplx);
    LOAD_FUNC(signal_pimp);
    LOAD_FUNC(sin_cmplx);
    LOAD_FUNC(sqrt_cmplx);
    LOAD_FUNC(trapint);
    LOAD_FUNC(trapint_cmplx);
    LOAD_FUNC(writebin);
    LOAD_FUNC(writetxt);


    #ifdef WIN_OS
        return (void*)handle;   
	exit_label:
		printf("function %s loading ERROR!\n", fname);
		if(handle)
			FreeLibrary(handle);
		return NULL;		
	#endif //WIN_OS


     #ifdef LINUX_OS
		return handle;		
	exit_label:
		printf("function %s loading ERROR!\n", fname);
		if(handle)
			dlclose(handle);
		return NULL;		
	#endif //LINUX_OS	
}







void dspl_free(void* handle)
{
	#ifdef WIN_OS
		FreeLibrary((HINSTANCE)handle);
	#endif //WIN_OS


	
	#ifdef LINUX_OS
		dlclose(handle);
	#endif //LINUX_OS	
	
}




