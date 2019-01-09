#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "dspl.h"

#define ORD 5
#define N   1000
#define K   1024

int main()
{
    void* handle;           // DSPL handle
    handle = dspl_load();   // Load DSPL function

    double a[ORD+1], b[ORD+1];
    double Rs = 50;
    double w[N], mag[N], phi[N], tau[N];
    double t[K], h[K];
    fft_t pfft;
   
    int k;
    int res = cheby2_ap(Rs, ORD, b, a);
    if(res != RES_OK)
        printf("error code = 0x%8x\n", res);
    
    for(k = 0; k < ORD+1; k++)
        printf("b[%2d] = %9.5f     a[%2d] = %9.5f\n", k, b[k], k, a[k]);


    logspace(-2.0, 2.0, N , DSPL_SYMMETRIC, w);
    freqs_resp(b, a, ORD, w, N, DSPL_FLAG_LOG|DSPL_FLAG_UNWRAP, mag, phi, tau);

    writetxt(w, mag, N, "dat/cheby2_ap_test_mag.txt");    
    writetxt(w, phi, N, "dat/cheby2_ap_test_phi.txt");
    writetxt(w, tau, N, "dat/cheby2_ap_test_tau.txt");

   memset(&pfft, 0, sizeof(fft_t));
   freqs2time(b, a, ORD, 10.0, K, &pfft, t,h);
   writetxt(t, h, K, "dat/cheby2_ap_test_time.txt");

    dspl_free(handle);      // free dspl handle

    return 0;
}


