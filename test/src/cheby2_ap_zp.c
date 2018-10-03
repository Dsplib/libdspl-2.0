#include <stdio.h>
#include <stdlib.h>
#include "dspl.h"

#define ORD 	9
#define N	1024
#define SCALE 1.0

int main()
{
	void* handle;           // DSPL handle
	handle = dspl_load();   // Load DSPL function
	
	complex_t z[ORD], p[ORD];
	
	int nz, np, k;
	double Rs = 40.0;
	
	printf("\nORD = 9\n");
	int res = cheby2_ap_zp(ORD, Rs, z, &nz, p, &np);
	if(res != RES_OK)
		printf("error code = 0x%8x\n", res);
	
	printf("\nChebyshev type 2 zeros:\n");
	for(k = 0; k < nz; k++)
		printf("(%9.3f,  %9.3f)\n", k, SCALE*RE(z[k]), SCALE*IM(z[k]));
	
	printf("\nChebyshev type 2 poles:\n");
	for(k = 0; k < np; k++)
		printf("(%9.3f,  %9.3f)\n", k, SCALE*RE(p[k]), SCALE*IM(p[k]));
	
	
	double alpha[N], sigma[N], omega[N];
	double rs, es, beta, shb, chb, ca, sa, den;
	char fname[64];
	int m;
	
	linspace(0, M_2PI, N, DSPL_SYMMETRIC, alpha);
	
	
	for(m = 1; m < 12; m++)
	{
		rs = (double)m*10.0;
		es = sqrt(pow(10.0, rs*0.1) - 1.0);
		beta   = asinh(es)/(double)ORD;
	
		chb = cosh(beta);
		shb = sinh(beta);
	
		for(k = 0; k < N; k++)
		{
			ca = cos(alpha[k]);
			sa = sin(alpha[k]);
			den =  (ca*ca*chb*chb + sa*sa*shb*shb);		
			sigma[k] = -SCALE * sa * shb/den;
			omega[k] =  SCALE * ca * chb/den;
		}
		
		sprintf(fname, "dat/cheby2_ap_poles_ord9_rs%d.txt", m*10);		
		writetxt(sigma, omega, N, fname);
	}
	
	printf("\nORD = 8\n");
	res = cheby2_ap_zp(ORD-1, Rs, z, &nz, p, &np);
	if(res != RES_OK)
		printf("error code = 0x%8x\n", res);
	
	printf("\nChebyshev type 2 zeros:\n");
	for(k = 0; k < nz; k++)
		printf("(%9.3f,  %9.3f)\n", k, SCALE*RE(z[k]), SCALE*IM(z[k]));
	
	printf("\nChebyshev type 2 poles:\n");
	for(k = 0; k < np; k++)
		printf("(%9.3f,  %9.3f)\n", k, SCALE*RE(p[k]), SCALE*IM(p[k]));
	
	for(m = 1; m < 12; m++)
	{
		rs = (double)m*10.0;
		es = sqrt(pow(10.0, rs*0.1) - 1.0);
		beta   = asinh(es)/(double)(ORD-1);
		
		chb = cosh(beta);
		shb = sinh(beta);
		
		for(k = 0; k < N; k++)
		{
			ca = cos(alpha[k]);
			sa = sin(alpha[k]);
			den =  (ca*ca*chb*chb + sa*sa*shb*shb);		
			sigma[k] = -SCALE * sa * shb/den;
			omega[k] =  SCALE * ca * chb/den;
		}
		
		sprintf(fname, "dat/cheby2_ap_poles_ord8_rs%d.txt", m*10);		
		writetxt(sigma, omega, N, fname);
	}
	
	
	dspl_free(handle);      // free dspl handle
	return 0;
}
