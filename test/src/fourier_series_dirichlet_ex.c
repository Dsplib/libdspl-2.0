#include <stdio.h>
#include <string.h>
#include <math.h>
#include "dspl.h"

#define N  10000

int main(int argc, char* argv[])
{
  double t[N],  s[N];
	void* handle;
	int k;
	handle = dspl_load();
	if(!handle)
	{
		printf("cannot to load libdspl!\n");
		return 0;
	}


	linspace(-M_2PI, M_2PI, N, DSPL_PERIODIC, t);
	for(k = 0; k < N; k++)
		s[k] = sin(1.0/sin(t[k]));
	writetxt(t, s, N, "dat/fourier_series_signal0.txt");

	linspace(-0.2, 0.2, N, DSPL_PERIODIC, t);
	for(k = 0; k < N; k++)
		s[k] = sin(1.0/sin(t[k]));
	writetxt(t, s, N, "dat/fourier_series_signal1.txt");
	// remember to free the resource
	dspl_free(handle);


	return 0;
}
