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

    double r[N*M];
    complex_t c[N*M];

    linspace(0, N*M, N*M, DSPL_PERIODIC, r);
    matrix_print(r, N, M, "R", "%8.2f");

    linspace(0, N*M*2, N*M*2, DSPL_PERIODIC, (double*)c);
    matrix_print_cmplx(c, N, M, "C", "%8.2f%+8.2fi");

    matrix_eye(r, N, M);
    matrix_print(r, N, M, "I", "%8.2f");

    matrix_eye_cmplx(c, N, M);
    matrix_print_cmplx(c, N, M, "I", "%8.2f%+8.2fi");

    dspl_free(handle);      /* free dspl handle  */
    return 0;
}

