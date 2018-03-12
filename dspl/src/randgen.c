

#include <stdlib.h>      
#include <math.h>
#include <time.h>

#include "dspl.h"


#define DSPL_RAND_MOD_X1  2147483647 
#define DSPL_RAND_MOD_X2  2145483479


/**************************************************************************************************
Uniform random numbers generator
***************************************************************************************************/
int randu(double* x, int n)
{
 	int k,m;
 	unsigned int x1[4], x2[4], y;

 	if(!x)
		return  ERROR_PTR;
 	if(n<1)
		return ERROR_SIZE;

 	x1[1] = rand();
 	x2[1] = rand();
 	x1[2] = rand();
 	x2[2] = rand();
 	x1[3] = rand();
 	x2[3] = rand();
 	for(k = 0; k<n; k++)
 	{
 		x1[0] = (63308 * x1[2] - 183326*x1[3]) % DSPL_RAND_MOD_X1;
 		x2[0] = (86098 * x2[1] - 539608*x2[3]) % DSPL_RAND_MOD_X2;
 		y = (x1[0] - x2[0]) %  DSPL_RAND_MOD_X1;
 		for(m = 3; m > 0; m--)
 		{
 			x1[m] = x1[m-1];
 			x2[m] = x2[m-1];
 		}

 		x[k] = (double)y/DSPL_RAND_MOD_X1;
 	}
	
	return RES_OK;
}





/**************************************************************************************************
Gaussian random numbers generator
***************************************************************************************************/
int randn(double* x, int n, double mu, double sigma)
{
	int k, m;
	double x1[128], x2[128];
	int res;
	if(!x)
		return ERROR_PTR;

 	if(n<1)
		return ERROR_SIZE;

	if(sigma < 0.0)
		return ERROR_RAND_SIGMA;	

	k=0;
	while(k < n)
	{
		res = randu(x1, 128);
		if(res != RES_OK)
			goto exit_label;
		
		res = randu(x2, 128);
		if(res != RES_OK)
			goto exit_label;
		m = 0 ;
		while(k<n && m < 128)
		{
			x[k] = sqrt(-2.0*log(x1[m]))*cos(M_2PI*x2[m])*sigma + mu;
			k++;
			m++;
			if(k<n && m < 128)
			{
				x[k] = sqrt(-2.0*log(x1[m]))*sin(M_2PI*x2[m])*sigma + mu;
				k++;
				m++;			
			}			
		}		
	}
	
	res = RES_OK;		
exit_label:
	return res;
}

