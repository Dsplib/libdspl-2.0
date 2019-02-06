#include <stdio.h>
#include <string.h>
#include "dspl.h"

#define N  800
#define K  19
void scale_vector(double* x, int n, double a)
{
  int k;
  for(k = 0; k < n; k++)
    x[k] *= a;
}


void magphi(double* x, int n, double* mag, double* phi)
{
  int k;
  for(k = 0; k < n/2; k++)
  {
    mag[k] = fabs(x[k]);
    phi[k] = x[k] < 0.0 ? 1.0 : 0.0;
  }
  
  for(k = n/2; k < n; k++)
  {
    mag[k] = fabs(x[k]);
    phi[k] = x[k] < 0.0 ? -1.0 : 0.0;
  }
}




int main(int argc, char* argv[])
{

  double w[N];
  double S[N];
  
  double mag[N];
  double phi[N];
  
  
  double wn[K];
  double Sn[K];

  int n;
  
	void* handle;
  handle = dspl_load();
	if(!handle)
	{
		printf("cannot to load libdspl!\n");
		return 0;
	}
  
  linspace(-8.0, 8.0, N, DSPL_SYMMETRIC, w);
  linspace(-8.0, 8.0, K, DSPL_SYMMETRIC, wn);
  
  sinc(w, N, M_PI*0.5, S);
  scale_vector(S, N, 2);
  
  sinc(wn, K, M_PI*0.5, Sn);
  scale_vector(Sn, K, 2);
  
  writetxt(w+10,  S+10,  N - 20, "dat/pimp_env.txt");
  writetxt(wn+1,  Sn+1,  K - 2,  "dat/pimp_spec.txt");
  
  magphi(S, N, mag, phi);
  writetxt(w+10,  mag+10,  N - 20, "dat/pimp_mag_env.txt");
  writetxt(w+10,  mag+10,  N - 20, "dat/pimp_mag_env_shift.txt");  
  writetxt(w+10,  phi+10,  N - 20, "dat/pimp_phi_env.txt");
  
  magphi(Sn, K, mag, phi);
  
  writetxt(wn+1,  mag+1,  K-2,  "dat/pimp_mag_spec.txt");
  writetxt(wn+1,  mag+1,  K-2,  "dat/pimp_mag_spec_shift.txt");
  writetxt(wn+1,  phi+1,  K-2,  "dat/pimp_phi_spec.txt");
  
  
  magphi(S, N, mag, phi);
  for(n = 0; n < N; n++)
    phi[n] -= w[n]/M_PI;  
  writetxt(w+10,  phi+10,  N - 20, "dat/pimp_phi_env_shift.txt");
  
  
  magphi(Sn, K, mag, phi);
  for(n = 0; n < K; n++)
    phi[n] -= wn[n]/M_PI;  
  writetxt(wn+1,  phi+1,  K - 2, "dat/pimp_phi_spec_shift.txt");
    
	// remember to free the resource
	dspl_free(handle);
  return 0;
}
