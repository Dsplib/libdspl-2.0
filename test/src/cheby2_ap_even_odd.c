#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "dspl.h"

#define ORD 9
#define N   1000
#define K   1024

int main()
{
	void* handle;           // DSPL handle
	handle = dspl_load();   // Load DSPL function
	
	double a[ORD+1], b[ORD+1];
	double Rs = 40.0;
	double w[N], mag[N], phi[N], tau[N];
	double t[K], h[K];
	fft_t pfft;
	
	int k;
	int res = cheby2_ap(Rs, ORD-1, b, a);
	if(res != RES_OK)
		printf("error code = 0x%8x\n", res);
	
	for(k = 0; k < ORD; k++)
		printf("b[%2d] = %9.3f     a[%2d] = %9.3f\n", k, b[k], k, a[k]);
	
	
	logspace(-2.0, 2.0, N , DSPL_SYMMETRIC, w);
	freqs_resp(b, a, ORD-1, w, N, 
		   DSPL_FLAG_LOG|DSPL_FLAG_UNWRAP, mag, phi, tau);
	
	writetxt(w, mag, N, "dat/cheby2_ap_ord6_mag.txt");    
	writetxt(w, phi, N, "dat/cheby2_ap_ord6_phi.txt");
	writetxt(w, tau, N, "dat/cheby2_ap_ord6_tau.txt");
	
	memset(&pfft, 0, sizeof(fft_t));
	freqs2time(b, a, ORD-1, 50.0, K, &pfft, t,h);
	writetxt(t, h, K, "dat/cheby2_ap_ord6_time.txt");
	
	
	
	res = cheby2_ap(Rs, ORD, b, a);
	if(res != RES_OK)
		printf("error code = 0x%8x\n", res);
	
	for(k = 0; k < ORD+1; k++)
		printf("b[%2d] = %9.3f     a[%2d] = %9.3f\n", k, b[k], k, a[k]);
	
	freqs_resp(b, a, ORD, w, N, 
		   DSPL_FLAG_LOG|DSPL_FLAG_UNWRAP, mag, phi, tau);
	
	writetxt(w, mag, N, "dat/cheby2_ap_ord7_mag.txt");    
	writetxt(w, phi, N, "dat/cheby2_ap_ord7_phi.txt");
	writetxt(w, tau, N, "dat/cheby2_ap_ord7_tau.txt");
	
	
	memset(&pfft, 0, sizeof(fft_t));
	freqs2time(b, a, ORD, 50.0, K, &pfft, t,h);
	writetxt(t, h, K, "dat/cheby2_ap_ord7_time.txt");
	
	
	
	dspl_free(handle);      // free dspl handle
	
	return 0;
}


