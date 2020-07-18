#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "dspl.h"

#define N 4
#define M 3

int main()
{

    void* handle;           /* DSPL handle         */
    handle = dspl_load();   /* Load DSPL function  */

    double    r[N*M], p[N*M];
    complex_t c[N*M], d[N*M];

    linspace(0, N*M, N*M, DSPL_PERIODIC, r);
    matrix_print(r, N, M, "R", "%8.2f");

    matrix_transpose(r, N, M, p);
    matrix_print(p, M, N, "P", "%8.2f");

    linspace(0, N*M*2, N*M*2, DSPL_PERIODIC, (double*)c);
    matrix_print_cmplx(c, N, M, "C", "%8.2f%+8.2fi");

    matrix_transpose_cmplx(c, N, M, d);
    matrix_print_cmplx(d, M, N, "D", "%8.2f%+8.2fi");

    matrix_transpose_hermite(c, N, M, d);
    matrix_print_cmplx(d, M, N, "D", "%8.2f%+8.2fi");


    dspl_free(handle);      /* free dspl handle  */

    return 0;

}

