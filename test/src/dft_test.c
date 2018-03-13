#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "dspl.h"

#define N 16
int main()
{
	void* handle;
	handle = dspl_load();
	
    complex_t x[N];
    complex_t y[N];
	//complex_t z[N];
	
	
    for(int k = 0; k < N; k++)
    {
        RE(x[k]) = (double)k;
		IM(x[k]) = 0.0;
    }

    dft_cmplx(x,N,y);

	for(int k = 0; k < N; k++)
        printf("y[%2d] = %9.3f%9.3f\n", k, RE(y[k]), IM(y[k]));
	
    // remember to free the resource
    dspl_free(handle);
    return 0;
}


