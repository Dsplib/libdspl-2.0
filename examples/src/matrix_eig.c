#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "dspl.h"

#define N  3

int main()
{
    void* handle;           /* DSPL handle         */
    handle = dspl_load();   /* Load DSPL function  */

    complex_t a[N*N] = {{1.0, 0.0}, {1.0, 0.0}, {0.0, 0.0},
                        {2.0, 0.0}, {0.0, 0.0}, {1.0, 0.0},
                        {3.0, 0.0}, {0.0, 0.0}, {0.0, 0.0}};

    complex_t v[N];
    int err, info;
    matrix_print_cmplx(a, N, N, "A", "%8.2f%+8.2fi");

    err = matrix_eig_cmplx(a, N, v, &info);
    if(err!=RES_OK)
    {
        printf("ERROR CODE: 0x%.8x, info = %d\n", err, info);
    }

    matrix_print_cmplx(v, N, 1, "v", "%10.6f%+10.6fi");

    dspl_free(handle);      /* free dspl handle  */
    return 0;
}

