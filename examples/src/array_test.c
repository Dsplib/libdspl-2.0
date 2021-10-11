#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "dspl.h"


  
int main()
{
    void* handle;             /* DSPL handle */
    handle = dspl_load();     /* Load DSPL function */

    double a[5] = {0.0, 1.0, 2.0, 3.0, 4.0};
    double b[5] = {5.0, 6.0, 7.0, 8.0, 9.0};
    double c[10], d[5], r;
    int err, k, n;


    err = vector_dot(a, b, 5, &r);
    printf("\n\ndot product result: %d a^T * b = %f ", err, r);

    /* Concatenate arrays a and b. Result keeps to the array c */
    err = concat((void*)a, 5*sizeof(double),
                 (void*)b, 5*sizeof(double), (void*)c);
    printf("\n\nconcatenation result: %d\n\narray c = ", err);
    for(k = 0; k < 10; k++)
        printf("%6.1f", c[k]);

    /* Decimate array c 2 times. Result keeps to the array d */
    err = decimate(c, 10, 2, d, &n);
    printf("\n\ndecimation result: %d\n\narray d = ", err);
    for (k = 0; k < n; k++)
        printf("%6.1f", d[k]);

    /* find max abs value */
    double am[5] = {0.0, 2.0, -5.0, 4.0, 2.0};
    double m;
    int ind;
    err = find_max_abs(am, 5, &m, &ind);
    printf("\n\nmax absolute value: %8.1f (index %d)", m, ind);

    /* Flip in place */
    printf("\n\nFlipip function test:\n\n");
    double x[5] = {0.0, 1.0, 2.0, 3.0, 4.0};
    int i;
    for(i = 0; i < 5; i++)
        printf("%6.1f  ", x[i]);
    flipip(x, 5);
    printf("\n");
    for(i = 0; i < 5; i++)
        printf("%6.1f  ", x[i]);


    printf("\n\nflipip_cmplx function test:\n\n");
    complex_t y[5] = {{0.0, 0.0}, {1.0, 1.0},
                      {2.0, 2.0}, {3.0, 3.0}, {4.0, 4.0}};
    for(i = 0; i < 5; i++)
        printf("%6.1f%+.1fj  ", RE(y[i]), IM(y[i]));
    flipip_cmplx(y, 5);
    printf("\n");
    for(i = 0; i < 5; i++)
        printf("%6.1f%+.1fj  ", RE(y[i]), IM(y[i]));

    /* free dspl handle */
    dspl_free(handle);

    return err;
}


