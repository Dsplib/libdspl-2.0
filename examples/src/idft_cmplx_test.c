#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "dspl.h"

#define N 16
int main()
{

    void* handle;           /* DSPL handle            */
    handle = dspl_load();   /* Load DSPL function     */
    complex_t x[N];         /* complex input signal   */
    complex_t y[N];         /* DFT                    */
    complex_t z[N];         /* IDFT                   */
    int k;

    /* Fill input signal */
    for(k = 0; k < N; k++)
    {
        RE(x[k]) = (double)k;
        IM(x[k]) = 0.0;
    }

    /* DFT */
    dft_cmplx(x,N,y);

    /* IDFT */
    idft_cmplx(y,N,z);

    /* print result                                     */
    for(k = 0; k < N; k++)
        printf("x[%2d] = %9.3f%+9.3fj,    z[%2d] = %9.3f%+9.3f\n",
        k, RE(x[k]), IM(x[k]), k, RE(z[k]), IM(z[k]));


    dspl_free(handle);  /* remember to free the resource */
    return 0;
}


