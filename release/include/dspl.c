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

p_cheby_poly1           cheby_poly1		;
p_cheby_poly2           cheby_poly2		;
p_conv                  conv            ;
p_conv_cmplx            conv_cmplx      ;
p_dft                 	dft				;
p_dft_cmplx           	dft_cmplx		;
p_filter_iir            filter_iir      ;
p_goertzel              goertzel        ;
p_goertzel_cmplx        goertzel_cmplx  ;
p_linspace              linspace        ;
p_logspace              logspace        ;
p_polyval               polyval         ;
p_polyval_cmplx         polyval_cmplx   ;


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


    
    LOAD_FUNC(cheby_poly1);
    LOAD_FUNC(cheby_poly2);
    LOAD_FUNC(conv);
    LOAD_FUNC(conv_cmplx);
    LOAD_FUNC(dft);
    LOAD_FUNC(dft_cmplx);
    LOAD_FUNC(filter_iir);
    LOAD_FUNC(goertzel);
    LOAD_FUNC(goertzel_cmplx);
    LOAD_FUNC(linspace);
    LOAD_FUNC(logspace);
    LOAD_FUNC(polyval);
    LOAD_FUNC(polyval_cmplx);




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




