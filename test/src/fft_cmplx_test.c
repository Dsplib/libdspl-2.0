#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "dspl.h"

#define N 16
int main()
{
    void* handle;           // DSPL handle
    handle = dspl_load();   // Load DSPL function

    complex_t x[N];         // complex input signal
    complex_t y[N];         // DFT
    fft_t pfft;             // FFT object
    
    //clear fft object
    memset(&pfft, 0, sizeof(fft_t));
    
    // Create FFT object   
    fft_create(&pfft, N);       
    
    for(int k = 0; k < N; k++)
    {
        RE(x[k]) = (double)k;
        IM(x[k]) = 0.0;
    }
    
    //FFT
    fft_cmplx(x, N, &pfft, y);

    for(int k = 0; k < N; k++)
        printf("y[%2d] = %9.3f%9.3f\n", k, RE(y[k]), IM(y[k]));

    fft_free(&pfft);        // clear FFT object
    dspl_free(handle);      // free dspl handle
    return 0;
}


