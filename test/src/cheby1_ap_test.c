#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "dspl.h"

#define ORD 4
#define N   1000


int main()
{
    void* handle;           // DSPL handle
    handle = dspl_load();   // Load DSPL function

    double    a[ORD+1], b[ORD+1];
    double Rp = 3.0;
   
    int k;
    int res = cheby1_ap(Rp, ORD, b, a);
    if(res != RES_OK)
        printf("error code = 0x%8x\n", res);
    
    for(k = 0; k < ORD+1; k++)
        printf("b[%2d] = %9.3f     a[%2d] = %9.3f\n", k, b[k], k, a[k]);


    double w[N], hdb[N];
    complex_t h[N];

    logspace(-2.0, 2.0, N , DSPL_SYMMETRIC, w);
    freqs(b, a, ORD, w, N, h);
    for(k = 0; k < N; k++)
        hdb[k] = 10.0 * log10(ABSSQR(h[k]));

    writetxt(w, hdb, N, "dat/cheby1_ap_test_mag.txt");    
    
    dspl_free(handle);      // free dspl handle

    system("gnuplot -p  gnuplot/cheby1_ap_test.plt");

    return 0;
}


