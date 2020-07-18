#include <stdio.h>
#include <stdlib.h>
#include "dspl.h"

#define N 14

int main()
{
    void* handle;           /* DSPL handle        */
    handle = dspl_load();   /* Load DSPL function */

    double k[N];
    int i, err;

    err = ellip_landen(0.93, N, k);
    if(err != RES_OK)
    {
        printf("Error code: %8x\n", err);
        return err;
    }

    printf(" i%8sk[i]\n\n", "");
    for(i = 1; i < N; i++)
        printf("%2d  %11.3e\n", i, k[i]);

    dspl_free(handle);      /* free dspl handle */
    return 0;
}


