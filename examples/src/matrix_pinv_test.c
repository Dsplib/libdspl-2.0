#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "dspl.h"

#define N  2
#define M  6



int main(int argc, char* argv[])
{
    void* handle;           /* DSPL handle         */
    handle = dspl_load();   /* Load DSPL function  */

    /* Matrix 
    A = [ 1   2   2;
          2   4   5;
          2   0   1;
          2  -1   0;];
    in array a by columns
    */
    double   a[N*M] = { 1, 2, 2, 5,    2, 4, 4, -1,   0, 0, 3, -2};
    double inv[M*N]; /* left orthogonal matrix U              */

    
    int err, info, mn;
    
    /* print input matrix */
    matrix_print(a, N, M, "A", "%8.2f");

    /* SVD decomposition A = U*S*V^T */
    /*-----------------------------------------------------*/
    err = matrix_pinv(a, N, M, NULL, inv, &info);
    if(err != RES_OK)
        printf("err = %.8x  info = %d\n", err, info);
    
    /* Print SVD decomposition */
    matrix_print(inv,  M, N, "inv(A)",   "%8.8f");




    dspl_free(handle);      /* free dspl handle  */
    return 0;
}
