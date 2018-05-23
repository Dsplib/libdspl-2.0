/*
* Copyright (c) 2015-2018 Sergey Bakhurin
* Digital Signal Processing Library [http://dsplib.org]
*
* This file is part of libdspl-2.0.
*  
* is free software: you can redistribute it and/or modify
* it under the terms of the GNU Lesser  General Public License as published by
* the Free Software Foundation, either version 3 of the License, or
* (at your option) any later version.
*
* DSPL is distributed in the hope that it will be useful,
* but WITHOUT ANY WARRANTY; without even the implied warranty of
* MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
* GNU General Public License for more details.
*
* You should have received a copy of the GNU Lesser General Public License
* along with Foobar.  If not, see <http://www.gnu.org/licenses/>.
*/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "dspl.h"





/******************************************************************************
Analog Normalized Butterworth filter
*******************************************************************************/
int DSPL_API butter_ap(double rp, int ord, double* b, double* a)
{
	int res;
	complex_t *z = NULL;
	complex_t *p = NULL;
	int nz, np;

	if(rp < 0.0)
		return ERROR_FILTER_RP;
	if(ord < 1)
		return ERROR_FILTER_ORD;
	if(!a || !b)
		return ERROR_PTR;

	z = (complex_t*) malloc(ord*sizeof(complex_t));
	p = (complex_t*) malloc(ord*sizeof(complex_t));
	
	
	res = butter_ap_zp(ord, rp, z, &nz, p, &np);
	if(res != RES_OK)
		goto exit_label;
	
	res = filter_zp2ab(z, nz, p, np, ord, b, a);
	if(res != RES_OK)
		goto exit_label;

	b[0] = a[0];


exit_label:
	if(z)
		free(z);
	if(p)
		free(p);
	return res;   
}








/******************************************************************************
Analog Normalized Butterworth filter zeros and poles  
*******************************************************************************/
int DSPL_API butter_ap_zp(int ord, double rp, complex_t *z, int* nz, 
			  				complex_t *p, int* np)
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
	*np = ord;
	*nz = 0;  
	return RES_OK;
}






/******************************************************************************
Analog Normalized Chebyshev type 1 filter  
*******************************************************************************/
int DSPL_API cheby1_ap(double rp, int ord, double* b, double* a)
{
	int res;
	complex_t *z = NULL;
	complex_t *p = NULL;
	int nz, np, k;
	complex_t h0 = {1.0, 0.0};
	double tmp;
	
	
	if(rp < 0.0)
		return ERROR_FILTER_RP;
	if(ord < 1)
		return ERROR_FILTER_ORD;
	if(!a || !b)
		return ERROR_PTR;

	z = (complex_t*) malloc(ord*sizeof(complex_t));
	p = (complex_t*) malloc(ord*sizeof(complex_t));
	
	
	res = cheby1_ap_zp(ord, rp, z, &nz, p, &np);
	if(res != RES_OK)
		goto exit_label;
	
	res = filter_zp2ab(z, nz, p, np, ord, b, a);
	if(res != RES_OK)
		goto exit_label;


	if(!(ord % 2))
		RE(h0) = 1.0 / pow(10.0, rp*0.05);
 
	for(k = 0; k < np; k++)
	{
		tmp	= CMRE(h0, p[k]);
		IM(h0) = CMIM(h0, p[k]);
		RE(h0) = tmp;
	}
		
	b[0] = fabs(RE(h0));

exit_label:
	if(z)
		free(z);
	if(p)
		free(p);
	return res;   
}





/******************************************************************************
Analog Normalized Chebyshev type 1 filter zeros and poles  
*******************************************************************************/
int DSPL_API cheby1_ap_zp(int ord, double rp, complex_t *z, int* nz, 
							complex_t *p, int* np)
{
	double theta; 
	double ep;
	double beta; 
	double shbeta;
	double chbeta;
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
	
	
	beta   = asinh(1.0/ep)/(double)ord;
	chbeta = cosh(beta);
	shbeta = sinh(beta);

	if(r)
	{
		RE(p[ind]) = -shbeta;
		IM(p[ind]) =  0.0;  
		ind++;	  
	}
	for(k = 0; k < L; k++)
	{
		theta = M_PI*(double)(2*k + 1)/(double)(2*ord);
		RE(p[ind]) = RE(p[ind+1]) = -shbeta *  sin(theta);
		IM(p[ind]) =	chbeta *  cos(theta);
		IM(p[ind+1]) = -IM(p[ind]);
		ind+=2;
	}  
	*np = ord;
	*nz = 0;  
	return RES_OK;
}




/******************************************************************************
 * Analog Normalized Chebyshev type 2 filter
 *******************************************************************************/
int DSPL_API cheby2_ap(double rs, int ord, double* b, double* a)
{
	int res;
	complex_t *z = NULL;
	complex_t *p = NULL;
	int nz, np;
	double norm;
	
	
	if(rs < 0.0)
		return ERROR_FILTER_RP;
	if(ord < 1)
		return ERROR_FILTER_ORD;
	if(!a || !b)
		return ERROR_PTR;
	
	z = (complex_t*) malloc(ord*sizeof(complex_t));
	p = (complex_t*) malloc(ord*sizeof(complex_t));
	
	
	res = cheby2_ap_zp(ord, rs, z, &nz, p, &np);
	if(res != RES_OK)
		goto exit_label;
		
	res = filter_zp2ab(z, nz, p, np, ord, b, a);
	if(res != RES_OK)
		goto exit_label;

	norm = a[0] / b[0];
	
	for(nz = 0; nz < ord+1; nz++)
		b[nz]*=norm;
	
exit_label:
	if(z)
		free(z);
	if(p)
		free(p);
	return res;
}




/******************************************************************************
Analog Normalized Chebyshev type 2 filter zeros and poles 
*******************************************************************************/
int DSPL_API cheby2_ap_zp(int ord, double rs, complex_t *z, int* nz, 
			  				complex_t *p, int* np)
{
	double es;
	int L, r, k;
	double beta;
	int iz, ip;
	
	double alpha;
	double chb, shb, sa, ca;
	double ssh2, cch2;
	
	if(rs < 0 || rs == 0)
		return ERROR_FILTER_RS;
	if(ord < 1)
		return ERROR_FILTER_ORD;
	if(!z || !p || !nz || !np)
		return ERROR_PTR;

	es = sqrt(pow(10.0, rs*0.1) - 1.0);
	r = ord % 2;
	L = (int)((ord-r)/2);
	
	beta   = asinh(es)/(double)ord;
	
	chb = cosh(beta);
	shb = sinh(beta);
	
	iz = ip = 0;
	
	if(r)
	{
		RE(p[0]) = -1.0 / sinh(beta);
		IM(p[0]) =  0.0;
		ip = 1;
	}
	
	for(k = 0; k < L; k++)
	{
		alpha = M_PI*(double)(2*k + 1)/(double)(2*ord);
		sa = sin(alpha);
		ca = cos(alpha);
		ssh2  = sa*shb;
		ssh2 *= ssh2;
		
		cch2  = ca*chb;
		cch2 *= cch2;
		
		RE(z[iz]) = RE(z[iz+1]) = 0.0;
		IM(z[iz]) = 1.0 / ca;
		IM(z[iz+1]) = -IM(z[iz]);
		iz+=2;
		
		RE(p[ip]) = RE(p[ip+1]) = -sa*shb / (ssh2 + cch2);
		IM(p[ip]) = ca*chb / (ssh2 + cch2);
		IM(p[ip+1]) = -IM(p[ip]);
		ip+=2;		
	} 
	*nz = iz;
	*np = ip;
	
	return RES_OK;
}










/******************************************************************************
 A *nalog Normalized Chebyshev type 2 filter zeros and poles 
 *******************************************************************************/
int DSPL_API ellip_ap_zp(int ord, double rp, double rs, 
			 	complex_t *z, int* nz, complex_t *p, int* np)
{
	double es, ep;
	int L, r, n, res;
	int iz, ip;
	double ke, k, u, t;
	complex_t tc, v0, jv0;
	
	
	if(rp < 0 || rp == 0)
		return ERROR_FILTER_RP;
	if(rs < 0 || rs == 0)
		return ERROR_FILTER_RS;
	if(ord < 1)
		return ERROR_FILTER_ORD;
	if(!z || !p || !nz || !np)
		return ERROR_PTR;
	
	es = sqrt(pow(10.0, rs*0.1) - 1.0);
	ep = sqrt(pow(10.0, rp*0.1) - 1.0);
	ke = ep / es;
	
	r = ord % 2;
	L = (int)((ord-r)/2);
	
	res = ellip_modulareq(rp, rs, ord, &k);
	if(res != RES_OK)
		return res;
	// v0
	RE(tc) = 0.0;
	IM(tc) = 1.0 / ep;
	
	ellip_asn_cmplx(&tc, 1, ke, &v0);
	
	t = RE(v0);
	RE(v0) = IM(v0) / (double)ord;
	IM(v0) = -t / (double)ord;
	
	RE(jv0) = -IM(v0);
	IM(jv0) =  RE(v0);
	
	
	iz = ip = 0;
	
	if(r)
	{
		res = ellip_sn_cmplx(&jv0, 1, k, &tc);
		if(res != RES_OK)
			return res;
		RE(p[0]) = -IM(tc);
		IM(p[0]) =  RE(tc);
		ip = 1;
	}
	
	for(n = 0; n < L; n++)
	{
		u = (double)(2 * n + 1)/(double)ord;
		
		res = ellip_cd(& u, 1, k, &t);
		if(res != RES_OK)
			return res;
		
		RE(z[iz]) = RE(z[iz+1]) = 0.0;
		IM(z[iz])   =  1.0/(k*t);
		IM(z[iz+1]) = -1.0/(k*t);
		iz+=2;
		
		RE(tc) = u - RE(jv0);
		IM(tc) =   - IM(jv0);
		
		res = ellip_cd_cmplx(&tc, 1, k, p+ip+1);
		if(res != RES_OK)
			return res;
		
		RE(p[ip]) = -IM(p[ip+1]);
		IM(p[ip]) =  RE(p[ip+1]);
		
		RE(p[ip+1]) =   RE(p[ip]);
		IM(p[ip+1]) =  -IM(p[ip]);
		
		ip+=2;
		
		
	} 
	*nz = iz;
	*np = ip;
	
	return RES_OK;
}




/******************************************************************************
Zeros and poles to filter coefficients recalc
*******************************************************************************/
int DSPL_API filter_zp2ab(complex_t *z, int nz, complex_t *p, int np, 
						int ord, double* b, double* a)
{
	complex_t *acc = NULL;
	int res;

	if(!z || !p || !b || !a)
		return ERROR_PTR;
	if(nz < 0 || np < 0)
		return ERROR_SIZE;
	if(nz > ord || np > ord)
		return ERROR_POLY_ORD;

	acc = (complex_t*) malloc((ord+1) * sizeof(complex_t));
	res = poly_z2a_cmplx(z, nz, ord, acc);
	if(res != RES_OK)
		goto exit_label; 
 
	res = cmplx2re(acc, ord+1, b, NULL);
	if(res != RES_OK)
		goto exit_label;

	res = poly_z2a_cmplx(p, np, ord, acc);
	if(res != RES_OK)
		goto exit_label;

	res = cmplx2re(acc, ord+1, a, NULL);
	if(res != RES_OK)
		goto exit_label;


exit_label:
	if(acc)
		free(acc);
	return res;
	
}






