#include <stdio.h>
#include <string.h>
#include <math.h>
#include "dspl.h"

#define NW  	60
#define NS	30
#define ORD	2

#define R	1.0
#define C	0.5
#define L	2.0

int main(int argc, char* argv[]) 
{
	
	double w[NW], sigma[NS], habs[NW*NS];
	void* handle;
	double b[ORD+1] = {1, 0, 0};
	double a[ORD+1] = {1, R*C, L*C};
	complex_t hs[NW*NS];
	complex_t  s[NW*NS];
	int k, n;
	
	handle = dspl_load();
	if(!handle)
	{
		printf("cannot to load libdspl!\n");
		return 0;
	}
	
	
	linspace(-2, 0,  NS, DSPL_SYMMETRIC, sigma);
	linspace(-2, 2,  NW, DSPL_SYMMETRIC, w);
	
	
	for(k = 0; k < NW; k++)
	{
		for(n = 0; n < NS; n++)
		{
			RE(s[k*NS + n]) = sigma[n];
			IM(s[k*NS + n]) = w[k];
		}
	}
	
	
	freqs_cmplx(b, a, ORD, s, NW*NS, hs);
	//cmplx2re(hs, N, hr, hi);

	for(k = 0; k < NW*NS; k++)
	{
		habs[k] = ABS(hs[k]) > 16.0 ? 16.0 : ABS(hs[k]) ;
		
	}
	writetxt_3d(sigma, NS, w, NW,  habs,  "dat/laplace_tf_3d_abs.txt");
	
	freqs(b,a,ORD, w,NW, hs);
	
	for(k = 0; k<NW; k++)
		habs[k] =  ABS(hs[k]); 	
	
	sigma[0] = 0.0;
	writetxt_3d(sigma, 1, w, NW,  habs,  "dat/laplace_tf_3d_hjw.txt");

	
	
	
	dspl_free(handle);
	
	return 0;
}
