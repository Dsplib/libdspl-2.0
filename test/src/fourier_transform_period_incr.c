#include <stdio.h>
#include <string.h>
#include <math.h>
#include "dspl.h"

#define N  10000
#define D  10
#define TAU 1
#define K 51

int main(int argc, char* argv[])
{
	double *t  = NULL,  *s  = NULL;
	double *td = NULL,  *sd = NULL;
	complex_t S[K];
	double w[K], mag[K];
	double T;
	void* handle;
	int k, n;



	handle = dspl_load();
	if(!handle)
	{
		printf("cannot to load libdspl!\n");
		return 0;
	}

	t = (double*)malloc(N*sizeof(double));
	s = (double*)malloc(N*sizeof(double));

	td = (double*)malloc(N/D*sizeof(double));
	sd = (double*)malloc(N/D*sizeof(double));

	linspace(-10, 10, N, DSPL_SYMMETRIC, t);
	decimate(t, N, D,td, &n);

	T = 2.0;
	signal_pimp(t, N, 2.0, TAU, 0.0, T, s);
	decimate(s, N, D, sd, &n);
	writetxt(td, sd, n, "dat/fourier_transform_period_2.0_time.txt");


	fourier_series_dec(t, s, N, T, K, w, S);
	for(n = 0; n < K; n++)
		mag[n] = T*ABS(S[n])/20.0;
	writetxt(w, mag, K, "dat/fourier_transform_period_2.0_freq.txt");


	T = 4.0;
	signal_pimp(t, N, 2.0, TAU, 0.0, T, s);
	decimate(s, N, D, sd, &n);
	writetxt(td, sd, n, "dat/fourier_transform_period_4.0_time.txt");


	fourier_series_dec(t, s, N, T, K, w, S);
	for(n = 0; n < K; n++)
		mag[n] = T*ABS(S[n])/20.0;
	writetxt(w, mag, K, "dat/fourier_transform_period_4.0_freq.txt");


	T = 8.0;
	signal_pimp(t, N, 2.0, TAU, 0.0, T, s);
	decimate(s, N, D, sd, &n);
	writetxt(td, sd, n, "dat/fourier_transform_period_8.0_time.txt");


	fourier_series_dec(t, s, N, T, K, w, S);
	for(n = 0; n < K; n++)
		mag[n] = T*ABS(S[n])/20.0;
	writetxt(w, mag, K, "dat/fourier_transform_period_8.0_freq.txt");


	// remember to free the resource
	dspl_free(handle);

	return 0;
}
