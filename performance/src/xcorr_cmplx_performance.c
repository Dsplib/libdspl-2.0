#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#include "dspl.h"

#define N 16777216
#define P 16777216

int main(int argc, char* argv[])
{
    void* hdspl;            /* DSPL handle        */
    void* hplot;            /* GNUPLOT handles    */
    
    random_t rnd = {0};     /* random structure   */
    
    hdspl = dspl_load();    /* Load DSPL function */
    
    complex_t *x = (complex_t*) malloc (N * sizeof(complex_t));
    complex_t *z = (complex_t*) malloc ((2*P+1) * sizeof(complex_t));
    complex_t *y = (complex_t*) malloc ((2*P+1) * sizeof(complex_t));
    
    int err, i, j, k = 1, size = N;
    clock_t  t;
    double   dt;
    
    err = random_init(&rnd, RAND_TYPE_MRG32K3A, NULL);
    if(err != RES_OK)
        goto exit_label;

    err = randn((double*)x, 2*N, 0.0, 1.0, &rnd);
    if(err != RES_OK)
        goto exit_label;
    
    for(j = 0; j < 22; j++)
    {
        t = clock();
        for(i = 0; i < k; i++)
            xcorr_cmplx(x, size, x, size, DSPL_XCORR_BIASED, P, y, NULL);
        t = clock() - t;
        dt = ((double)t)/CLOCKS_PER_SEC/(double)k;
        printf ("N = %8d  xcorr time:         %12.6f sec\n", size, dt);
        
        k *= 1.7;
        size /= 2;;
    }


exit_label:
    free(x);
    free(y);
    free(z);
    /* free dspl handle */
    dspl_free(hdspl);
    return err;
}

