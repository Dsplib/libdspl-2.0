#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "dspl.h"


#define N  2000
#define A  1
#define TAU  1
int main(int argc, char* argv[])
{
	double    t[N], w[N], Sr[N], St[N], Se[N], Sg[N], sigma;
	void*	  handle;
    int k;
	handle = dspl_load();
	if(!handle)
	{
		printf("cannot to load libdspl!\n");
		return 0;
	}

	linspace(-10, 10,  N, DSPL_PERIODIC, t);
	linspace(-20*M_PI, 20*M_PI,  N, DSPL_PERIODIC, w);

    
    for(k = 0; k < N; k++)
    {
        Sr[k] = sin(w[k]*TAU*0.5+1E-9) / (w[k]*TAU*0.5+1E-9);
        Sr[k] *= Sr[k];
        
        St[k] = sin(w[k]*TAU*0.5+1E-9)/(w[k]*TAU*0.5+1E-9);
        St[k] *= St[k];
        St[k] *= St[k];
        sigma = 0.5;
        Sg[k] = w[k]*sigma;
        Sg[k] *= w[k]*sigma / 4.0;
        Sg[k] = exp(-Sg[k]);     
        Sg[k] *= Sg[k];
        sigma = 2;
        Se[k] = sigma*sigma/(w[k]*w[k] + sigma*sigma);
        
        Se[k] *= Se[k];       
  
    }     
	  writetxt(w, Sr, N, "dat/fourier_transform_energy_r_lin.txt");
    writetxt(w, St, N, "dat/fourier_transform_energy_t_lin.txt");
    writetxt(w, Se, N, "dat/fourier_transform_energy_e_lin.txt");
    writetxt(w, Sg, N, "dat/fourier_transform_energy_g_lin.txt");
    
    
    for(k = 0; k < N; k++)
    {
        Sr[k] = 10.0*log10(Sr[k]);
        Sg[k] = 10.0*log10(Sg[k]);  
        St[k] = 10.0*log10(St[k]);
        Se[k] = 10.0*log10(Se[k]);           
           
  
    }
    
    
    writetxt(w, Sr, N, "dat/fourier_transform_energy_r_log.txt");
    writetxt(w, St, N, "dat/fourier_transform_energy_t_log.txt");
    writetxt(w, Se, N, "dat/fourier_transform_energy_e_log.txt");
    writetxt(w, Sg, N, "dat/fourier_transform_energy_g_log.txt");
    


    // remember to free the resource
	dspl_free(handle);

}
