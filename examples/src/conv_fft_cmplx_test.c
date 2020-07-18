#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "dspl.h"

#define N 15
#define M 5
int main()
{
    void* handle;           /* DSPL handle        */
    handle = dspl_load();   /* Load DSPL function */
    complex_t a[N], b[M], c[N+M-1], d[N+M-1];
    fft_t pfft;
    int n;

    linspace(0, 2*N, 2*N, DSPL_PERIODIC, (double*)a);
    linspace(0, 2*M, 2*M, DSPL_PERIODIC, (double*)b);
    memset(&pfft, 0, sizeof(fft_t));

    conv_fft_cmplx(a, N, b, M, &pfft, 8, c);
    conv_cmplx(a, N, b, M, d);

    /* print result */
    for(n = 0; n < N+M-1; n++)
    {
        printf("c[%3d] = %9.2f%+9.2fj    ", n, RE(c[n]), IM(c[n]));
        printf("d[%3d] = %9.2f%+9.2fj  \n", n, RE(d[n]), IM(d[n]));
    }

    fft_free(&pfft);        /* free fft structure memory */
    dspl_free(handle);      /* free dspl handle          */
    return 0;
}





