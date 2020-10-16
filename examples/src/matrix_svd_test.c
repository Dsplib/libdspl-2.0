#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "dspl.h"

#define N  4
#define M  3


int main(int argc, char* argv[])
{
    void* handle;           /* DSPL handle         */
    handle = dspl_load();   /* Load DSPL function  */

    /* Matrix 
    A = [ 1   2   3;
          3   4   5;
          2   0   1;
          2  -1   0;];
    in array a by columns
    */
    double a[N*M] = { 1, 3, 2, 2,   2, 4, 0, -1,   3, 5, 1, 0};
    double u[N*N]; /* left orthogonal matrix U              */
    double s[M];   /* Singular values diagonal              */
    double vt[M*M];/* transposed left orthogonal matrix V^T */

    double ur[N*M] = {0}; /* matrix UR = U*S */
    double ar[N*M];       /* AR = UR * V^T   */ 
    
    int err, info, i, j, mn;
    
    /* print input matrix */
    matrix_print(a, N, M, "A", "%8.2f");

    /* SVD decomposition A = U*S*V^T */
    /*-----------------------------------------------------*/
    err = matrix_svd(a, N, M, u, s, vt, &info);
    if(err != RES_OK)
        printf("err = %.8x  info = %d\n", err, info);
    
    /* Print SVD decomposition */
    matrix_print(u,  N, N, "U",   "%8.4f");
    matrix_print(vt, M, M, "V^T", "%8.4f");
    matrix_print(s,  M, 1, "S",   "%8.6f");


    /* SVD reconstruction AR = U*S*V^T */
    /*-----------------------------------------------------*/

    /* step 1: UR = U*S. 
       We can process this by columns 
       because S is diagonal 
    */
    mn = (M > N) ? N : M;
    for(i = 0; i < mn; i++)
        for(j = 0; j < N; j++)
            ur[j + i*N] = u[j + i*N] * s[i];

    /* step 2: AR = U * S * V^T */
    matrix_mul(ur, N, M, vt, M, M, ar);
    matrix_print(ar, N, M, "AR = U*S*V^T",   "%8.2f");

    dspl_free(handle);      /* free dspl handle  */
    return 0;
}
