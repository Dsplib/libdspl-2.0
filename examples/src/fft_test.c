#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "dspl.h"

/* FFT size */
#define N 14

int main()
{
    void* handle;           /* DSPL handle              */
    handle = dspl_load();   /* Load libdspl             */
    double    x[N];         /* Input signal array       */
    complex_t y[N];         /* Output signal array      */
    fft_t pfft = {0};       /* FFT object (fill zeros)  */
    int k;


    /* Fill FFT structure                               */
    fft_create(&pfft, N);

    /* Fill input signal x[k] = k                        */
    for(k = 0; k < N; k++)
        x[k] = (double)k;

    /* FFT                                              */
    fft(x, N, &pfft, y);

    /* print result                                     */
    for(k = 0; k < N; k++)
        printf("y[%2d] = %9.3f%9.3f\n", k, RE(y[k]), IM(y[k]));

    fft_free(&pfft);        /* Clear fft_t object       */
    dspl_free(handle);      /* Clear DSPL handle        */
    return 0;
}



