#include <stdio.h>
#include <string.h>
#include <math.h>
#include "dspl.h"

#define N  	200
#define ORD	5
#define SCALEH 2.0
#define SCALEP 0.35
#define SCALES 2.0

int main(int argc, char* argv[]) 
{
	
	double w[N], sigma[N], hr[N], hi[N];
	void* handle;
	double b[ORD+1], a[ORD+1];
	//complex_t p[ORD], z[ORD];
	complex_t hs[N];
	int k;
	
	handle = dspl_load();
	if(!handle)
	{
		printf("cannot to load libdspl!\n");
		return 0;
	}
	
	
	linspace(-2.6, 2.6,  N, DSPL_SYMMETRIC, w);
	
	butter_ap(2, ORD, b, a);
	
	
	freqs(b,a, ORD, w, N, hs);
	//cmplx2re(hs, N, hr, hi);

	for(k = 0; k < N; k++)
	{
		hr[k] = RE(hs[k]) * SCALES;
		hi[k] = IM(hs[k]) * SCALES;
	}
	
	memset(sigma, 0, N*sizeof(double));
	
	writetxt_3dline(sigma,w, hr, N, "dat/laplace_tf_planes_3d_re.txt");
	writetxt_3dline(sigma,w, hi, N, "dat/laplace_tf_planes_3d_im.txt");
	
	
	writetxt(w, hr, N, "dat/laplace_tf_planes_2d_re.txt");
	writetxt(w, hi, N, "dat/laplace_tf_planes_2d_im.txt");
	
	freqs_resp(b,a, ORD, w, N, DSPL_FLAG_UNWRAP, hr, hi, NULL);
	for(k = 0; k < N; k++)
	{
		hr[k] *= SCALEH;
		hi[k] = (hi[k] + M_2PI) * SCALEP;
	}
	
	writetxt(w, hr, N, "dat/laplace_tf_planes_2d_abs.txt");
	writetxt(w, hi, N, "dat/laplace_tf_planes_2d_phi.txt");
	
	
	
	dspl_free(handle);
	
	return 0;
}
