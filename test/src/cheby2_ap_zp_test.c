#include <stdio.h>
#include <stdlib.h>
#include "dspl.h"

#define ORD 7

int main()
{
    void* handle;           // DSPL handle
    handle = dspl_load();   // Load DSPL function

    complex_t z[ORD], p[ORD];    
    int nz, np, k;
    double Rs = 40.0;
    
    int res = cheby2_ap_zp(ORD, Rs, z, &nz, p, &np);
    if(res != RES_OK)
        printf("error code = 0x%8x\n", res);
    
    printf("\nChebyshev type 2 zeros:\n");
    for(k = 0; k < nz; k++)
        printf("z[%2d] = %9.3f  %9.3f j\n", k, RE(z[k]), IM(z[k]));

    printf("\nChebyshev type 2 poles:\n");
    for(k = 0; k < np; k++)
        printf("p[%2d] = %9.3f  %9.3f j\n", k, RE(p[k]), IM(p[k]));

    dspl_free(handle);      // free dspl handle
    return 0;
}


