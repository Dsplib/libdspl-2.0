#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "dspl.h"

#define N 13
#define M 7
int main()
{
    void* handle;           /* DSPL handle        */
    handle = dspl_load();   /* Load DSPL function */
    double a[N], b[M], c[N+M-1], d[N+M-1];
    fft_t pfft;
    int n, err;

    linspace(0, N, N, DSPL_PERIODIC, a);
    linspace(0, M, M, DSPL_PERIODIC, b);
    memset(&pfft, 0, sizeof(fft_t));

    err = conv_fft(a, N, b, M, &pfft, 16, c);
    printf("conv_fft error: 0x%.8x\n", err);

    err = conv(a, N, b, M, d);
    printf("conv error:     0x%.8x\n", err);

    /* print result */
    for(n = 0; n < N+M-1; n++)
        printf("c[%3d] = %9.2f    d[%3d] = %9.2f\n", n, c[n], n, d[n]);

    fft_free(&pfft);        /* free fft structure memory */
    dspl_free(handle);      /* free dspl handle          */
    return 0;
}





