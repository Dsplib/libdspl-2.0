/*
* Copyright (c) 2015-2017 Sergey Bakhurin
* Digital Signal Processing Library [http://dsplib.org]
*
* This file is part of DSPL.
*  
* is free software: you can redistribute it and/or modify
* it under the terms of the GNU General Public License as published by
* the Free Software Foundation, either version 3 of the License, or
* (at your option) any later version.
*
* DSPL is distributed in the hope that it will be useful,
* but WITHOUT ANY WARRANTY; without even the implied warranty of
* MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
* GNU General Public License for more details.
*
* You should have received a copy of the GNU General Public License
* along with Foobar.  If not, see <http://www.gnu.org/licenses/>.
*/

#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "dspl.h"




int butter_ap_zp(int ord, double rp, complex_t *z, int* nz, complex_t *p, int* np)
{
  	double alpha;
	double theta; 
    double ep; 
    int r;
	int L;
    int ind = 0, k;

    if(rp < 0 || rp == 0)
		return ERROR_FILTER_RP;
	if(ord < 1)
		return ERROR_FILTER_ORD;
	if(!z || !p || !nz || !np)
		return ERROR_PTR;

    ep = sqrt(pow(10.0, rp*0.1) - 1.0);
	r = ord % 2;
	L = (int)((ord-r)/2);

    alpha = pow(ep, -1.0/(double)ord);
    if(r)
    {
        RE(p[ind]) = -alpha;
        IM(p[ind]) = 0.0;  
        ind++;      
    }
    for(k = 0; k < L; k++)
    {
        theta = M_PI*(double)(2*k + 1)/(double)(2*ord);
        RE(p[ind]) = RE(p[ind+1]) = -alpha *  sin(theta);
        IM(p[ind]) =   alpha *  cos(theta);
        IM(p[ind+1]) =  -alpha *  cos(theta);
        ind+=2;
    }    
    return RES_OK;
}




/**************************************************************************************************
Analog Normalized Butterworth filter  
***************************************************************************************************/
int butter_ap(double Rp, int ord, double* b, double* a)
{
	double ep;
	double p[3] = {0.0, 0.0, 1.0};
	double alpha;
	double teta;
	double *acc = NULL;
	int r;
	int L;
	int n;
	int k;
	int res;
	
	if(Rp < 0 || Rp == 0)
		return ERROR_FILTER_RP;
	if(ord < 1)
		return ERROR_FILTER_ORD;
	if(!b || !a)
		return ERROR_PTR;

	acc = (double*)malloc((ord+1)*sizeof(double));

	memset(acc, 0,  (ord+1)*sizeof(double));
	memset(a,   0,  (ord+1)*sizeof(double));
	memset(b,   0,  (ord+1)*sizeof(double));


	ep = sqrt(pow(10.0, Rp*0.1) - 1.0);
	r = ord % 2;
	L = (int)((ord-r)/2);
	
	b[0] = 1.0/ep;
	
	alpha = pow(ep, -1.0/(double)ord);
	/* first pole according to filter order */
	if(r)
	{
		/* we have one real pole if filter order is odd */
		a[0] = alpha;
		a[1] = 1.0;
		k = 2;
	}
	else
	{
		/* all poles are complex if filter order is even   */
		a[0] = 1.0;
		k = 1;
	}
	
	/* coeff calculation */
	for(n = 0; n < L; n++)
	{
		memcpy(acc, a, k*sizeof(double));
		teta = M_PI*(double)(2*n + 1)/(double)(2*ord);
		p[0] = alpha * alpha;
		p[1] = 2.0 * alpha * sin(teta);
		
		res = conv(p, 3, acc, k, a);
		if(res!=RES_OK)			
			goto exit_label;	
	
		k+=2;
	}
    
exit_label:
	if(acc)
		free(acc);
	return res;
}



/**************************************************************************************************
Analog Normalized Chebyshev type 1 filter  
***************************************************************************************************/
int cheby1_ap(double Rp, int ord, double* b, double* a)
{
	double ep;
	double alpha;
	double beta; 
	double shbeta;
	double chbeta;
	double sigma;
	double omega;
	double gain;
	
	
	double p[3] = {0.0, 0.0, 1.0};
	
	double *acc = NULL;
	int r;
	int L;
	int n;
	int k;
	int res;
	
	if(Rp < 0 || Rp == 0)
		return ERROR_FILTER_RP;
	if(ord < 1)
		return ERROR_FILTER_ORD;
	if(!b || !a)
		return ERROR_PTR;
	
	acc = (double*)malloc((ord+1)*sizeof(double));

	memset(acc, 0,  (ord+1)*sizeof(double));
	memset(a,   0,  (ord+1)*sizeof(double));
	memset(b,   0,  (ord+1)*sizeof(double));

	ep = sqrt(pow(10.0, Rp*0.1) - 1.0);
	r = ord % 2;
	L = (int)((ord-r)/2);
	
	beta   = asinh(1.0/ep)/(double)ord;
	chbeta = cosh(beta);
	shbeta = sinh(beta);
	
	/* first pole according to filter order */
	if(r)
	{ 
		/* we have one real pole if filter order is odd */
		a[0] = -shbeta;
		a[1] = 1.0;
		k = 2;
		gain = shbeta;
	}
	else
	{
		/* all poles are complex if filter order is even   */
		a[0] = 1.0;
		k = 1;
		gain = 1.0/sqrt(1.0 + ep*ep);
	}
	
	
	/* coeff calculation */
	for(n = 0; n < L; n++)
	{
		memcpy(acc, a, k*sizeof(double));
		alpha = M_PI*(double)(2*n + 1)/(double)(2*ord);
		sigma = -sin(alpha)*shbeta;
		omega =  cos(alpha)*chbeta;
		
		p[0] = sigma * sigma + omega * omega;
		gain*=p[0];
		p[1] = -2.0 * sigma;
		
		res = conv(p, 3, acc, k, a);
		if(res!=RES_OK)
			goto exit_label;
		k+=2;
	}
	b[0] = gain;
	
	res = RES_OK;
exit_label:
	if(acc)
		free(acc);
	return res;
}



/**************************************************************************************************
Analog Normalized Chebyshev type II filter  
***************************************************************************************************/
int cheby2_ap(double Rs, int ord, double *b, double *a)
{
	double es, *acc = NULL, *bcc = NULL;
	double alpha, beta, sigma, omega, sh, ch, so2;
	double p[3] = {0.0, 0.0, 1.0};
	double q[3] = {0.0, 0.0, 1.0};
	double gain;
	int r, L, n, kp;
	int res;
	
	if(Rs < 0 || Rs == 0)
		return ERROR_FILTER_RS;
	if(ord < 1)
		return ERROR_FILTER_ORD;
	if(!b || !a)
		return ERROR_PTR;
	
	acc = (double*)malloc((ord+1)*sizeof(double));
	bcc = (double*)malloc((ord+1)*sizeof(double));

	memset(acc, 0,  (ord+1)*sizeof(double));
	memset(bcc, 0,  (ord+1)*sizeof(double));
	memset(a,   0,  (ord+1)*sizeof(double));
	memset(b,   0,  (ord+1)*sizeof(double));

	es = sqrt(pow(10.0, Rs*0.1) - 1.0);
	r = ord % 2;
	L = (int)((ord-r)/2);

	beta = asinh(es)/(double)ord;

	sh = sinh(beta);
	ch = cosh(beta);

	b[0] = 1.0;
	/* first pole according to filter order */
	if(r)
	{
		/* we have one real pole if filter order is odd  */
		a[0] = 1.0/sh;
		a[1] = 1.0;
		kp = 2;
	}
	else
	{
		/* all poles are TCOMPLEX if filter order is even  */
		a[0] = 1.0;
		kp = 1;
	}

	/* coeff calculation */
	for(n = 0; n < L; n++)
	{
		memcpy(acc, a,  kp * sizeof(double));
		memcpy(bcc, b,  kp * sizeof(double));

		alpha = M_PI*(double)(2*n+1)/(double)(2*ord);
		
		sigma = sh*sin(alpha);
		omega = ch*cos(alpha);

		so2 = 1.0 / (sigma*sigma+omega*omega);

		q[0] = 1.0/(cos(alpha)*cos(alpha));		
		p[0] = so2;
		p[1] = 2.0*sigma*so2;

		res = conv(p, 3, acc, kp, a);
		if(res!=RES_OK)
			goto exit_label;
		
		res = conv(q, 3, bcc, kp, b);
		if(res!=RES_OK)
			goto exit_label;
		
		kp+=2;
		
	}
	gain = b[0] / a[0];
	for(n = 0; n < ord+1; n++)
		b[n] /= gain;
	
exit_label:
	if(acc)
		free(acc);
	if(bcc)
		free(bcc);
	return res;

}




