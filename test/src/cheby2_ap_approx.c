#include <stdio.h>
#include <string.h>
#include <math.h>
#include "dspl.h"

#define N  1000
#define ORD 5
int main(int argc, char* argv[]) 
{	
	double w[N],  h[N];
	double a[ORD+1], b[ORD+1];
	void* handle;
	int k;
	handle = dspl_load();
	if(!handle)
	{
		printf("cannot to load libdspl!\n");
		return 0;
	}
	
	
	linspace(0, 2.5,  N, DSPL_PERIODIC, w);
	
	cheby2_ap(20.0, ORD, b,a);
	freqs_resp(b, a, ORD, w, N, 0, h, NULL, NULL);
	
	for( k =0; k< N; k++)
	{
		h[k] *= 6.0;
		w[k] *= 4.0;
	}
	writetxt(w, h, N, "dat/cheby2_approx.txt");   
	
	
	// remember to free the resource
	dspl_free(handle);
	
	return 0;
}
