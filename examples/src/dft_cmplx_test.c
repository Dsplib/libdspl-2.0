#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "dspl.h"

#define N 16
int main()
{

    void* handle;           /* DSPL handle           */
    handle = dspl_load();   /* Load DSPL function    */
    complex_t y[N];         /* DFT                   */

    complex_t x[N];         /* complex input signal  */
    int k;

    for(k = 0; k < N; k++)
    {
        RE(x[k]) = (double)k;
        IM(x[k]) = 0.0;
    }

    dft_cmplx(x,N,y);

    for(k = 0; k < N; k++)
        printf("y[%2d] = %9.3f%9.3f\n", k, RE(y[k]), IM(y[k]));


    dspl_free(handle);  /* remember to free the resource */
    return 0;
}


