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

    double a[ORD+1], b[ORD+1];
    double Rp = 3.0;
    double w[N], mag[N], phi[N], tau[N];

   
    int k;
    int res = butter_ap(Rp, ORD, b, a);
    if(res != RES_OK)
        printf("error code = 0x%8x\n", res);
    
    for(k = 0; k < ORD+1; k++)
        printf("b[%2d] = %9.3f     a[%2d] = %9.3f\n", k, b[k], k, a[k]);


    logspace(-2.0, 2.0, N , DSPL_SYMMETRIC, w);
    freqs_resp(b, a, ORD, w, N, DSPL_FLAG_LOG|DSPL_FLAG_UNWRAP, mag, phi, tau);

    writetxt(w, mag, N, "dat/butter_ap_test_mag.txt");    
    writetxt(w, phi, N, "dat/butter_ap_test_phi.txt");
    writetxt(w, tau, N, "dat/butter_ap_test_tau.txt");

    dspl_free(handle);      // free dspl handle

    res = system("gnuplot gnuplot/butter_ap_test.plt");

    return 0;
}


