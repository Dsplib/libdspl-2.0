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



/**************************************************************************************************
inverse cd function
***************************************************************************************************/
int ellip_acd(double* w, int n, double k, double* u)
{
    double lnd[ELLIP_ITER], t;    
    int i, m;    

    if(!u || !w)
        return ERROR_PTR;
    if(n<1)
        return ERROR_SIZE;
    if(k < 0.0 || k>= 1.0)
        return ERROR_ELLIP_MODULE; 

    ellip_landen(k,ELLIP_ITER, lnd);
    

    for(m = 0; m < n; m++)
    {
        u[m] = w[m];
        for(i = 1; i < ELLIP_ITER; i++)
        {   
            t = lnd[i-1]*u[m];
            t *= t;
            t = 1.0 + sqrt(1.0 - t);
            u[m]  = 2.0 * u[m] / (t+t*lnd[i]);
        }
        u[m] = 2.0 * acos(u[m]) / M_PI;
    }
    return RES_OK;
}





/**************************************************************************************************
inverse cd function
***************************************************************************************************/
int ellip_acd_cmplx(complex_t* w, int n, double k, complex_t* u)
{
    double lnd[ELLIP_ITER], t;
    complex_t tmp0, tmp1;    
    int i, m;    

    if(!u || !w)
        return ERROR_PTR;
    if(n<1)
        return ERROR_SIZE;
    if(k < 0.0 || k>= 1.0)
        return ERROR_ELLIP_MODULE; 

    ellip_landen(k,ELLIP_ITER, lnd);
    

    for(m = 0; m < n; m++)
    {
        RE(u[m]) = RE(w[m]);
        IM(u[m]) = IM(w[m]);
        for(i = 1; i < ELLIP_ITER; i++)
        {   
            RE(tmp0) = lnd[i-1]*RE(u[m]);
            IM(tmp0) = lnd[i-1]*IM(u[m]);  
            RE(tmp1) = 1.0 - CMRE(tmp0, tmp0);
            IM(tmp1) =     - CMIM(tmp0, tmp0);
            
            sqrt_cmplx(&tmp1, 1, &tmp0);
            RE(tmp0) += 1.0;

            RE(tmp1) = RE(tmp0) * (1.0 + lnd[i]);
            IM(tmp1) = IM(tmp0) * (1.0 + lnd[i]);           

            t = 2.0 / ABSSQR(tmp1);

            RE(tmp0) = t * CMCONJRE(tmp1, u[m]);
            IM(tmp0) = t * CMCONJIM(tmp1, u[m]);

            RE(u[m]) = RE(tmp0);
            IM(u[m]) = IM(tmp0);

        }
        acos_cmplx(&tmp0, 1, u+m);
        t = 2.0 / M_PI;
        RE(u[m]) *= t;
        IM(u[m]) *= t;
    }
    return RES_OK;
}





/**************************************************************************************************
inverse sn function
***************************************************************************************************/
int ellip_asn(double* w, int n, double k, double* u)
{
    double lnd[ELLIP_ITER], t;    
    int i, m;    

    if(!u || !w)
        return ERROR_PTR;
    if(n<1)
        return ERROR_SIZE;
    if(k < 0.0 || k>= 1.0)
        return ERROR_ELLIP_MODULE; 

    ellip_landen(k,ELLIP_ITER, lnd);
    

    for(m = 0; m < n; m++)
    {
        u[m] = w[m];
        for(i = 1; i < ELLIP_ITER; i++)
        {   
            t = lnd[i-1]*u[m];
            t *= t;
            t = 1.0 + sqrt(1.0 - t);
            u[m]  = 2.0 * u[m] / (t+t*lnd[i]);
        }
        u[m] = 2.0 * asin(u[m]) / M_PI;
    }
    return RES_OK;
}





/**************************************************************************************************
inverse sn function
***************************************************************************************************/
int ellip_asn_cmplx(complex_t* w, int n, double k, complex_t* u)
{
    double lnd[ELLIP_ITER], t;
    complex_t tmp0, tmp1;    
    int i, m;    

    if(!u || !w)
        return ERROR_PTR;
    if(n<1)
        return ERROR_SIZE;
    if(k < 0.0 || k>= 1.0)
        return ERROR_ELLIP_MODULE; 

    ellip_landen(k,ELLIP_ITER, lnd);
    

    for(m = 0; m < n; m++)
    {
        RE(u[m]) = RE(w[m]);
        IM(u[m]) = IM(w[m]);
        for(i = 1; i < ELLIP_ITER; i++)
        {   
            RE(tmp0) = lnd[i-1]*RE(u[m]);
            IM(tmp0) = lnd[i-1]*IM(u[m]);  
            RE(tmp1) = 1.0 - CMRE(tmp0, tmp0);
            IM(tmp1) =     - CMIM(tmp0, tmp0);
            
            sqrt_cmplx(&tmp1, 1, &tmp0);
            RE(tmp0) += 1.0;

            RE(tmp1) = RE(tmp0) * (1.0 + lnd[i]);
            IM(tmp1) = IM(tmp0) * (1.0 + lnd[i]);           

            t = 2.0 / ABSSQR(tmp1);

            RE(tmp0) = t * CMCONJRE(tmp1, u[m]);
            IM(tmp0) = t * CMCONJIM(tmp1, u[m]);

            RE(u[m]) = RE(tmp0);
            IM(u[m]) = IM(tmp0);

        }
        asin_cmplx(&tmp0, 1, u+m);
        t = 2.0 / M_PI;
        RE(u[m]) *= t;
        IM(u[m]) *= t;
    }
    return RES_OK;
}






/**************************************************************************************************
Elliptic cd function
***************************************************************************************************/
int ellip_cd(double* u, int n, double k, double* y)
{
    double lnd[ELLIP_ITER];    
    int i, m;    

    if(!u || !y)
        return ERROR_PTR;
    if(n<1)
        return ERROR_SIZE;
    if(k < 0.0 || k>= 1.0)
        return ERROR_ELLIP_MODULE; 

    ellip_landen(k,ELLIP_ITER, lnd);
    

    for(m = 0; m < n; m++)
    {
        y[m] = cos(u[m] * M_PI * 0.5);
        for(i = ELLIP_ITER-1; i>0; i--)
        {
            y[m] = (1.0 + lnd[i]) / (1.0 / y[m] + lnd[i]*y[m]);
        }
    }
    return RES_OK;
}






/**************************************************************************************************
Elliptic cd function
***************************************************************************************************/
int ellip_cd_cmplx(complex_t* u, int n, double k, complex_t* y)
{
    double lnd[ELLIP_ITER], t;    
    int i, m; 
    complex_t tmp;   

    if(!u || !y)
        return ERROR_PTR;
    if(n<1)
        return ERROR_SIZE;
    if(k < 0.0 || k>= 1.0)
        return ERROR_ELLIP_MODULE; 

    ellip_landen(k,ELLIP_ITER, lnd);
    

    for(m = 0; m < n; m++)
    {
        RE(tmp) = RE(u[m]) * M_PI * 0.5;
        IM(tmp) = IM(u[m]) * M_PI * 0.5;
        
        cos_cmplx(&tmp, 1, y+m);

        for(i = ELLIP_ITER-1; i>0; i--)
        {
            t = 1.0 / ABSSQR(y[m]);
    
            RE(tmp) =  RE(y[m]) * t + RE(y[m]) * lnd[i];
            IM(tmp) = -IM(y[m]) * t + IM(y[m]) * lnd[i];    

            t = (1.0 + lnd[i]) / ABSSQR(tmp);

            RE(y[m]) =   RE(tmp) * t;
            IM(y[m]) =  -IM(tmp) * t;
    
        }
    }
    return RES_OK;
}






/**************************************************************************************************
Landen transform
***************************************************************************************************/
int ellip_landen(double k, int n, double* y)
{
    int i;
    y[0] = k;

    if(!y)
        return ERROR_PTR;
    if(n < 1)
        return ERROR_SIZE;
    if(k < 0.0 || k>= 1.0)
        return ERROR_ELLIP_MODULE;

    for(i = 1; i < n; i++)
    {
        y[i] = y[i-1] / (1.0 + sqrt(1.0 - y[i-1] * y[i-1]));
        y[i] *= y[i];
    }

    return RES_OK;
}





/**************************************************************************************************
Elliptic modular equation
***************************************************************************************************/
int ellip_modulareq(double rp, double rs, int ord, double *k)
{
    double ep, es, ke, kp, t, sn;
    int i, L, r;

    if(rp < 0 || rp == 0)
		return ERROR_FILTER_RP;
    if(rs < 0 || rs == 0)
		return ERROR_FILTER_RS;
	if(ord < 1)
		return ERROR_FILTER_ORD;
    if(!k)
        return ERROR_PTR;


    ep = sqrt(pow(10.0, rp*0.1)-1.0);
    es = sqrt(pow(10.0, rs*0.1)-1.0);
    
    ke = ep/es;

    ke = sqrt(1.0 - ke*ke);

    r = ord % 2;
    L = (ord-r)/2;

    kp = 1.0;
    for(i = 0; i < L; i++)
    {
        t = (double)(2*i+1) / (double)ord;
        ellip_sn(&t, 1, ke, &sn);
        sn*=sn;
        kp *= sn*sn;
    }

    kp *= pow(ke, (double)ord);
    *k = sqrt(1.0 - kp*kp);

    return RES_OK;
    
}



/**************************************************************************************************
Elliptic rational function
***************************************************************************************************/
int ellip_rat(double* w, int n, int ord, double k, double* u)
{
    double t, xi, w2, xi2, k2;    
    int i, m, r, L;    

    if(!u || !w)
        return ERROR_PTR;
    if(n<1)
        return ERROR_SIZE;
    if(k < 0.0 || k>= 1.0)
        return ERROR_ELLIP_MODULE; 
    
    r = ord%2;
    L = (ord-r)/2;

    if(r)
        memcpy(u, w, n*sizeof(double));
    else
    {
        for(m = 0; m < n; m++)
        {
            u[m] = 1.0;
        }
    }

    k2 = k*k;
    for(i = 0; i < L; i++)
    {
        t = (double)(2*i+1) / (double)ord; 
        ellip_cd(&t, 1, k, &xi);
        xi2 = xi*xi;
        for(m = 0; m < n; m++)
        {
            w2 = w[m]*w[m];
            u[m] *= (w2 - xi2) / (1.0 - w2 * k2 * xi2);  
            u[m] *= (1.0 - k2*xi2) / (1.0 - xi2);  
        }
    }
    return RES_OK;
}






/**************************************************************************************************
Elliptic sn function
***************************************************************************************************/
int ellip_sn(double* u, int n, double k, double* y)
{
    double lnd[ELLIP_ITER];    
    int i, m;    

    if(!u || !y)
        return ERROR_PTR;
    if(n<1)
        return ERROR_SIZE;
    if(k < 0.0 || k>= 1.0)
        return ERROR_ELLIP_MODULE; 

    ellip_landen(k,ELLIP_ITER, lnd);
    

    for(m = 0; m < n; m++)
    {
        y[m] = sin(u[m] * M_PI * 0.5);
        for(i = ELLIP_ITER-1; i>0; i--)
        {
            y[m] = (1.0 + lnd[i]) / (1.0 / y[m] + lnd[i]*y[m]);
        }
    }
    return RES_OK;
}


/**************************************************************************************************
Elliptic sn function
***************************************************************************************************/
int ellip_sn_cmplx(complex_t* u, int n, double k, complex_t* y)
{
    double lnd[ELLIP_ITER], t;    
    int i, m; 
    complex_t tmp;   

    if(!u || !y)
        return ERROR_PTR;
    if(n<1)
        return ERROR_SIZE;
    if(k < 0.0 || k>= 1.0)
        return ERROR_ELLIP_MODULE; 

    ellip_landen(k,ELLIP_ITER, lnd);
    

    for(m = 0; m < n; m++)
    {
        RE(tmp) = RE(u[m]) * M_PI * 0.5;
        IM(tmp) = IM(u[m]) * M_PI * 0.5;
        
        sin_cmplx(&tmp, 1, y+m);

        for(i = ELLIP_ITER-1; i>0; i--)
        {
            t = 1.0 / ABSSQR(y[m]);
    
            RE(tmp) =  RE(y[m]) * t + RE(y[m]) * lnd[i];
            IM(tmp) = -IM(y[m]) * t + IM(y[m]) * lnd[i];    

            t = (1.0 + lnd[i]) / ABSSQR(tmp);

            RE(y[m]) =   RE(tmp) * t;
            IM(y[m]) =  -IM(tmp) * t;
    
        }
    }
    return RES_OK;
}


