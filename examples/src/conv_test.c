#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "dspl.h"


int main()
{
    void* handle;           /* DSPL handle        */
    handle = dspl_load();   /* Load DSPL function */
    complex_t ac[3] = {{0.0, 1.0}, {1.0, 1.0}, {2.0, 2.0}};
    complex_t bc[4] = {{3.0, 3.0}, {4.0, 4.0}, {5.0, 5.0}, {6.0, 6.0}};
    complex_t cc[6];

    double ar[3] = {1.0, 2.0, 3.0};
    double br[4] = {3.0, -1.0, 2.0, 4.0};
    double cr[6];

    int n;

    printf("\nconv\n--------------------------------\n");
    conv(ar, 3, br, 4, cr);
    for(n = 0; n < 6; n++)
        printf("cr[%d] = %5.1f\n", n, cr[n]);

    printf("\nconv_cmplx\n--------------------------------\n");
    conv_cmplx(ac, 3, bc, 4, cc);
    for(n = 0; n < 6; n++)
        printf("cc[%d] = %5.1f%+5.1fj\n", n, RE(cc[n]),IM(cc[n]));

    dspl_free(handle);      // free dspl handle
    return 0;
}





