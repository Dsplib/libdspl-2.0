#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "dspl.h"

/* DFT size */
#define N 16

int main()
{

    void* handle;           /* DSPL handle         */
    handle = dspl_load();   /* Load DSPL function  */
    double    x[N];         /* real input signal   */
    complex_t y[N];         /* DFT vector          */
    int k;

    /* fill input signal */
    for(k = 0; k < N; k++)
        x[k] = (double)k;

    /* DFT calculation */
    dft(x, N, y);

    /* Print result   */
    for(k = 0; k < N; k++)
        printf("y[%2d] = %9.3f%9.3f\n", k, RE(y[k]), IM(y[k]));

    /* remember to free the resource */
    dspl_free(handle);
    return 0;
}


