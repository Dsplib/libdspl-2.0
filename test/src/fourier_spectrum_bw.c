#include <stdio.h>
#include <string.h>
#include <math.h>
#include "dspl.h"

#define N  1024
#define K  8
#define DT 0.01
#define M 20

int main(int argc, char* argv[])
{
  double t[N*K],  s[N*K], z[N*K];
  double w[N*K],  S[N*K], Z[N*K];
  
  double amp[M] = {0.098, -0.160, -0.134, 0.259, 0.248,
                  -0.525,  0.124, -0.080, 0.179, 0.160,
                   0.071, -0.055, -0.119, 0.061, 0.297,
                   0.023, -0.347,  0.085, 0.046, 0.022};  
  
  double mu[M] = {478.0, 552.0, 672.0, 189.0, 268.0, 
                  504.0, 824.0, 482.0, 471.0, 643.0, 
                  685.0, 163.0, 870.0, 424.0, 662.0, 
                  611.0, 377.0, 549.0, 824.0, 731.0};
  
  double sigma[M] = {527.0, 1399.0, 1499.0,  543.0,  523.0,  
                    1465.0, 1175.0, 1535.0,  604.0, 1074.0,
                     890.0,  616.0,  563.0, 1201.0, 1884.0,
                    1745.0, 1465.0, 1885.0, 1767.0, 1635.0};
  
  
  
  complex_t SC[N*K];
  int k, m;
  double ms;
	void* handle;
  fft_t pfft;

	handle = dspl_load();
	if(!handle)
	{
		printf("cannot to load libdspl!\n");
		return 0;
	}


	linspace(0, (double)(N*K), N*K, DSPL_PERIODIC, t);
  memset(z, 0, N*K*sizeof(double));
  ms = 0.0;
  for(k = 0; k < N; k++)
  {
    for(m =0; m < M; m++)
    {
      s[k] += amp[m] * exp(-(t[k]-mu[m])*(t[k]-mu[m]) / sigma[m]); 
    }
		s[k] += 2.0 * exp(-(t[k]-(double)N * 0.5)*(t[k]-(double)N * 0.5) / 150.0); 

  }
  
  for(k = 503; k < 523; k++)
    z[k] = 1.8;
  
  for(k = 0; k < N; k++)
    t[k] *= DT;
  
  
  writetxt(t, s, N, "dat/dat_s.txt");
  writetxt(t, z, N, "dat/dat_z.txt");
  
  memset(&pfft, 0, sizeof(fft_t));
  fft_create(&pfft, N*K);
  fft(s, N*K, &pfft, SC);
   
  for(k = 0; k < N*K; k++)
  {    
    S[k] = ABS(SC[k])*DT;
    w[k] = -M_PI/DT + M_2PI/DT/(double)(N*K) * (double)k;
  }
  fft_shift(S, N*K, S);
  
  fft(z, N*K, &pfft, SC);
   
  for(k = 0; k < N*K; k++)
  {    
    Z[k] = ABS(SC[k])*DT;
  }
  fft_shift(Z, N*K, Z);
  
  
  writetxt(w, S, N*K, "dat/spec_s.txt");
  writetxt(w, Z, N*K, "dat/spec_z.txt");
  
	dspl_free(handle);

	return 0;
}
