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
p_blas_dscal            blas_dscal      ; 

p_cheby_poly1           cheby_poly1     ;
p_cheby_poly2           cheby_poly2     ;

p_dft                 	dft				;
p_dft_cmplx           	dft_cmplx		;

p_fft_create			fft_create		;
p_fft_destroy	   		fft_destroy		;
p_fft_cmplx	       		fft_cmplx		;
p_fft_shift	       		fft_shift		;


#endif //BUILD_LIB







void* dspl_load()
{
	#ifdef WIN_OS
		HINSTANCE handle;
		char* fname;
		handle = LoadLibrary(TEXT("libdspl.dll"));
		if (!handle)
		{
			printf("libdspl.dll loading ERROR!\n");
			return NULL;
		}
			

        fname = "blas_dscal";
		blas_dscal = (p_blas_dscal)GetProcAddress(handle, fname);
		if(!blas_dscal) goto exit_label;


        fname = "cheby_poly1";
		cheby_poly1 = (p_cheby_poly1)GetProcAddress(handle, fname);
		if(!cheby_poly1) goto exit_label;

        fname = "cheby_poly2";
		cheby_poly2 = (p_cheby_poly2)GetProcAddress(handle, fname);
		if(!cheby_poly2) goto exit_label;
	
		fname = "dft";
		dft = (p_dft)GetProcAddress(handle, fname);
		if(!dft) goto exit_label;
		
		fname = "dft_cmplx";
		dft_cmplx = (p_dft_cmplx)GetProcAddress(handle, fname);
		if(!dft_cmplx) goto exit_label;
		
		
		fname = "fft_create";
		fft_create = (p_fft_create)GetProcAddress(handle, fname);
		if(!fft_create) goto exit_label;
		
		fname = "fft_destroy";
		fft_destroy = (p_fft_destroy)GetProcAddress(handle, fname);
		if(!fft_destroy) goto exit_label;
		
		
		fname = "fft_cmplx";
		fft_cmplx = (p_fft_cmplx)GetProcAddress(handle, fname);
		if(!fft_cmplx) goto exit_label;
		
		fname = "fft_shift";
		fft_shift = (p_fft_shift)GetProcAddress(handle, fname);
		if(!fft_shift) goto exit_label;
		
		
		return (void*)handle;
	exit_label:
		printf("function %s loading ERROR!\n", fname);
		if(handle)
			FreeLibrary(handle);
		return NULL;		
	#endif //WIN_OS



	


	#ifdef LINUX_OS
		char* error;
		void *handle;
		char* fname;
		
		// open the *.so
		handle = dlopen ("./libdspl.so", RTLD_LAZY);
		if (!handle) 
		{
			printf("libdspl.so loading ERROR!\n");
			return NULL;
		}
		

        fname = "blas_dscal";
		blas_dscal = (p_blas_dscal)dlsym(handle, fname);
		if ((error = dlerror()) != NULL) goto exit_label;

        fname = "cheby_poly1";
		cheby_poly1 = (p_cheby_poly1)dlsym(handle, fname);
		if ((error = dlerror()) != NULL) goto exit_label;

        fname = "cheby_poly2";
		cheby_poly2 = (p_cheby_poly2)dlsym(handle, fname);
   		if ((error = dlerror()) != NULL) goto exit_label;
		
		
		fname = "dft";
		dft = (p_dft)dlsym(handle, fname);
		if ((error = dlerror()) != NULL) goto exit_label;
		
		
		fname = "dft_cmplx";
		dft_cmplx = (p_dft_cmplx)dlsym(handle, fname);
		if ((error = dlerror()) != NULL) goto exit_label;

        fname = "fft_create";
		fft_create = (p_fft_create)dlsym(handle, fname);
		if ((error = dlerror()) != NULL) goto exit_label;


        fname = "fft_destroy";
		fft_destroy = (p_fft_destroy)dlsym(handle, fname);
		if ((error = dlerror()) != NULL) goto exit_label;


        fname = "fft_cmplx";
		fft_cmplx = (p_fft_cmplx)dlsym(handle, fname);
		if ((error = dlerror()) != NULL) goto exit_label;


        fname = "fft_shift";
		fft_shift = (p_fft_shift)dlsym(handle, fname);
		if ((error = dlerror()) != NULL) goto exit_label;

		
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




