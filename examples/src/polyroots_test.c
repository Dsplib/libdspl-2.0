#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "dspl.h"

#define N  2

int main(int argc, char* argv[])
{
    void* hdspl;                     /* DSPL handle           */
    double a[N+1] = {2.0, 2.0, 1.0}; /* P(x) = 2 + 2x + x^2   */
    complex_t r[N] = {0};            /* roots                 */
    int err, n, info;
    hdspl = dspl_load();            /* Load DSPL functions    */
    if(!hdspl)
    {
        printf("libdspl loading error!\n");
        return -1;
    }


    /* roots calculation */
    err = polyroots(a, N, r, &info);
    printf("Error code: 0x%.8x\n", err);

    /* print roots */
    for(n = 0; n < N; n++)
        printf("r[%d] = % -8.5f% -8.5f j\n", n, RE(r[n]), IM(r[n]));

    /* free dspl handle  */
    dspl_free(hdspl);
    return 0;
}

