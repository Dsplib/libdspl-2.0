#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "dspl.h"

#define N  3
#define M  4

int main(int argc, char* argv[])
{
    /* Matrix 
    z = [ 0   0   1   1;
          0   1   0   0;
          0   1   0   0];
    in array a by columns
    */
    double z[N*M] = { 1, 0, 1,   0, 0, 0,   0, 1, 0,   1, 0, 1};
    double x[N] = {0.0, 1.0, 2.0};
    double y[M] = {0.0, 1.0, 2.0, 3.0};
    void* handle;           /* DSPL handle         */
    handle = dspl_load();   /* Load DSPL function  */

    contour2d_t contour = {0};

    contour2d(z, x, y, N, M, 0.5,&contour);

    dspl_free(handle);      /* free dspl handle  */
    return 0;
}
