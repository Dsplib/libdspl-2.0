#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "dspl.h"

#define N 4
#define M 3
#define K 5

int main()
{
    void* handle;           /* DSPL handle         */
    handle = dspl_load();   /* Load DSPL function  */

    double a[N*K], b[K*M], c[N*M];

    linspace(0, N*K, N*K, DSPL_PERIODIC, a);
    matrix_print(a, N, K, "A", "%8.2f");

    linspace(0, K*M, K*M, DSPL_PERIODIC, b);
    matrix_print(b, K, M, "B", "%8.2f");
    
    matrix_mul(a, N, K, b, K, M, c);
    matrix_print(c, N, M, "C", "%8.2f");

    dspl_free(handle);      /* free dspl handle  */
    return 0; 
}

