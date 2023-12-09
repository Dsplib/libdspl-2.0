#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "dspl.h"

#define N 6
int main(int argc, char* argv[])
{
    void* hdspl;            /* DSPL handle        */
    random_t rnd = {0};     /* random structure   */
    hdspl = dspl_load();    /* Load DSPL function */
    
    /* input vector */
    double x[N];
    
    double *y = NULL;
    int err, n, m, t;

    /* Generate random input vector */
    err = random_init(&rnd, RAND_TYPE_MRG32K3A, NULL);
    if(err != RES_OK)
        goto exit_label;

    err = randn(x, N, 0.0, 1.0, &rnd);
    if(err != RES_OK)
        goto exit_label;
    
    /* Write input vector x to the bin file x.bin */
    err = writebin(x, N, 1, DAT_DOUBLE, "x.bin");
    if(err != RES_OK)
        goto exit_label;
    
    /* Read x.bin to the y vector */
    err = readbin("x.bin", (void**)&y, &n, &m, &t);
    if(err != RES_OK)
        goto exit_label;

    /*Print x and y vectors */
    printf("n = %d, m = %d, t = %d\n", n, m, t);
    for(t = 0; t < n; t++)
        printf("%+8.4f, %+8.4f\n", x[t], y[t]);


exit_label:
    if(y)
        free(y);
    /* free dspl handle */
    dspl_free(hdspl);
    return 0;
}

