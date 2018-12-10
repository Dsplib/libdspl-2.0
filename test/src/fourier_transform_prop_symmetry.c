#include <stdio.h>
#include <string.h>
#include <math.h>
#include "dspl.h"

#define N  1000

int main(int argc, char* argv[])
{
	double t[N],  s[N], S[N], PHI[N], w[N];
	complex_t SC[N], sc[N] ;
	void* handle;
  int k;

	handle = dspl_load();
	if(!handle)
	{
		printf("cannot to load libdspl!\n");
		return 0;
	}


	linspace(-1, 10, N, DSPL_SYMMETRIC, t);
	linspace(-2.0*M_2PI, 2.0*M_2PI, N, DSPL_SYMMETRIC, w);

  for(k =0; k < N; k++)
  {
    IM(sc[k]) = 0.0;
    if(t[k]<0)
    {
      s[k] = RE(sc[k]) = 0.0;
    }
    else
    {
      s[k] = RE(sc[k]) = exp(-t[k]);
    }
  }

	writetxt(t, s, N/2, "dat/fourier_transform_prop_sym_time.txt");

  fourier_integral_cmplx(t, sc, N, N, w, SC);

  for(k = 0; k < N; k++)
  {
    S[k] = ABS(SC[k]);
    PHI[k] = ARG(SC[k]);
  }
  writetxt(w, S,   N, "dat/fourier_transform_prop_sym_mag.txt");
  writetxt(w, PHI, N, "dat/fourier_transform_prop_sym_phi.txt");

	// remember to free the resource
	dspl_free(handle);

	return 0;
}
