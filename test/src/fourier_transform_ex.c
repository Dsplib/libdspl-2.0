#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "dspl.h"


#define N  2000
#define A  1
int main(int argc, char* argv[])
{
	double    t[N], w[N], s[N], S[N], PHI[N], sigma;
	void*	  handle;
    int k;
	handle = dspl_load();
	if(!handle)
	{
		printf("cannot to load libdspl!\n");
		return 0;
	}

	linspace(-10, 10,  N, DSPL_PERIODIC, t);
	linspace(-4*M_PI, 4*M_PI,  N, DSPL_PERIODIC, w);

    sigma = 0.5;
    for(k = 0; k < N; k++)
    {
        s[k] = A * exp(-t[k]*t[k]/(sigma*sigma));
        S[k] = A * sigma * sqrt(M_PI) * exp(-w[k]*w[k]*sigma*sigma*0.25);   
    }     
	writetxt(t, s, N, "dat/fourier_transform_ex_gauss_time_0.5.txt");
    writetxt(w, S, N, "dat/fourier_transform_ex_gauss_freq_0.5.txt");    

    for(k = 0; k < N; k++)
    {
        s[k] = A * exp(-fabs(t[k]) * sigma);
        S[k] = 2 * A * sigma / (sigma*sigma + w[k]*w[k]);   
    }     
	writetxt(t, s, N, "dat/fourier_transform_ex_exp2_time_0.5.txt");
    writetxt(w, S, N, "dat/fourier_transform_ex_exp2_freq_0.5.txt"); 

    for(k = 0; k < N; k++)
    {
        s[k] = t[k] >=0 ? A * exp(-t[k] * sigma) : 0.0;
        S[k] = A  / sqrt(sigma*sigma + w[k]*w[k]);
        PHI[k] = -atan(w[k] / sigma);   
    }     
	writetxt(t, s, N, "dat/fourier_transform_ex_exp1_time_0.5.txt");
    writetxt(w, S, N, "dat/fourier_transform_ex_exp_freq_mag_0.5.txt");      
    writetxt(w, PHI, N, "dat/fourier_transform_ex_exp_freq_phi_0.5.txt"); 


    sigma = 1;
    for(k = 0; k < N; k++)
    {
        s[k] = A * exp(-t[k]*t[k]/(sigma*sigma));
        S[k] = A * sigma * sqrt(M_PI) * exp(-w[k]*w[k]*sigma*sigma*0.25);   
    }     
	writetxt(t, s, N, "dat/fourier_transform_ex_gauss_time_1.0.txt");
    writetxt(w, S, N, "dat/fourier_transform_ex_gauss_freq_1.0.txt");  


    for(k = 0; k < N; k++)
    {
        s[k] = A * exp(-fabs(t[k]) * sigma);
        S[k] = 2 * A * sigma / (sigma*sigma + w[k]*w[k]);   
    }     
	writetxt(t, s, N, "dat/fourier_transform_ex_exp2_time_1.0.txt");
    writetxt(w, S, N, "dat/fourier_transform_ex_exp2_freq_1.0.txt");  

    for(k = 0; k < N; k++)
    {
        s[k] = t[k] >=0 ? A * exp(-t[k] * sigma) : 0.0;
        S[k] = A  / sqrt(sigma*sigma + w[k]*w[k]);
        PHI[k] = -atan(w[k] / sigma);   
    }     
	writetxt(t, s, N, "dat/fourier_transform_ex_exp1_time_1.0.txt");
    writetxt(w, S, N, "dat/fourier_transform_ex_exp_freq_mag_1.0.txt");      
    writetxt(w, PHI, N, "dat/fourier_transform_ex_exp_freq_phi_1.0.txt");   



    sigma = 2;
    for(k = 0; k < N; k++)
    {
        s[k] = A * exp(-t[k]*t[k]/(sigma*sigma));
        S[k] = A * sigma * sqrt(M_PI) * exp(-w[k]*w[k]*sigma*sigma*0.25);   
    }     
	writetxt(t, s, N, "dat/fourier_transform_ex_gauss_time_2.0.txt");
    writetxt(w, S, N, "dat/fourier_transform_ex_gauss_freq_2.0.txt");    
    

    for(k = 0; k < N; k++)
    {
        s[k] = A * exp(-fabs(t[k]) * sigma);
        S[k] = 2 * A * sigma / (sigma*sigma + w[k]*w[k]);   
    }     
	writetxt(t, s, N, "dat/fourier_transform_ex_exp2_time_2.0.txt");
    writetxt(w, S, N, "dat/fourier_transform_ex_exp2_freq_2.0.txt");  
    
    for(k = 0; k < N; k++)
    {
        s[k] = t[k] >=0 ? A * exp(-t[k] * sigma) : 0.0;
        S[k] = A  / sqrt(sigma*sigma + w[k]*w[k]);
        PHI[k] = -atan(w[k] / sigma);   
    }     
	writetxt(t, s, N, "dat/fourier_transform_ex_exp1_time_2.0.txt");
    writetxt(w, S, N, "dat/fourier_transform_ex_exp_freq_mag_2.0.txt");      
    writetxt(w, PHI, N, "dat/fourier_transform_ex_exp_freq_phi_2.0.txt");  


    // remember to free the resource
	dspl_free(handle);
	return system("gnuplot -p  gnuplot/fourier_transform_ex_gauss.plt");;
}
