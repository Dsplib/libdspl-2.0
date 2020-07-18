#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "dspl.h"

int main(int argc, char* argv[])
{
    void* handle;           /*  DSPL handle        */


    complex_t x[3] = {{1.0, 2.0}, {3.0, 4.0}, {5.0, 6.0}};
    complex_t y[3];
    int k;

    handle = dspl_load();   /*  Load DSPL function */

    printf("\n\nacos_cmplx\n---------------------------------\n");
    acos_cmplx(x, 3, y);
    for(k = 0; k < 3; k++)
      printf("acos_cmplx(%.1f%+.1fj) = %.3f%+.3fj\n",
              RE(x[k]), IM(x[k]), RE(y[k]), IM(y[k]));

    printf("\n\nasin_cmplx\n---------------------------------\n");
    asin_cmplx(x, 3, y);
    for(k = 0; k < 3; k++)
        printf("asin_cmplx(%.1f%+.1fj) = %.3f%+.3fj\n",
                RE(x[k]), IM(x[k]), RE(y[k]), IM(y[k]));

    printf("\n\ncos_cmplx\n---------------------------------\n");
    cos_cmplx(x, 3, y);
    for(k = 0; k < 3; k++)
        printf("cos_cmplx(%.1f%+.1fj) = %9.3f%+9.3fj\n",
                RE(x[k]), IM(x[k]), RE(y[k]), IM(y[k]));

    printf("\n\nlog_cmplx\n---------------------------------\n");
    log_cmplx(x, 3, y);
    for(k = 0; k < 3; k++)
        printf("log_cmplx(%.1f%+.1fj) = %.3f%+.3fj\n",
                RE(x[k]), IM(x[k]), RE(y[k]), IM(y[k]));

    printf("\n\nsin_cmplx\n---------------------------------\n");
    sin_cmplx(x, 3, y);
    for(k = 0; k < 3; k++)
        printf("sin_cmplx(%.1f%+.1fj) = %9.3f%+9.3fj\n",
                RE(x[k]), IM(x[k]), RE(y[k]), IM(y[k]));

    printf("\n\nsqrt_cmplx\n---------------------------------\n");
    sqrt_cmplx(x, 3, y);
    for(k = 0; k < 3; k++)
        printf("sqrt_cmplx(%.1f%+.1fj) = %.3f%+.3fj\n",
                RE(x[k]), IM(x[k]), RE(y[k]), IM(y[k]));

    /* free dspl handle */
    dspl_free(handle);

    return 0;
}


