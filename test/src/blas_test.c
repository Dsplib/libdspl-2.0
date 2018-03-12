#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "dspl.h"

#define N 16
int main()
{
    void* handle;
    handle = dspl_load();

    double x[N];
    double a = 2.0;

    for(int k = 0; k < N; k++)
        x[k] = (double)k;

    blas_dscal(N, a, x, 1);

    for(int k = 0; k < N; k++)
        printf("x[%2d] = %9.3f \n", k, x[k]);

    // remember to free the resource
    dspl_free(handle);
    return 0;
}


