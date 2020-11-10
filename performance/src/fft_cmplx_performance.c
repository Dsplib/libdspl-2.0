#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#include "dspl.h"

#define NMAX 4194304
#define L 20
#define SIZE_FACTOR 2.3 



int fft_perf_cmplx(int* len, double* dlen, double* mflops)
{
    int k = (int)( 8 * pow(SIZE_FACTOR, L));
    int i, j, err;
    clock_t  t;
    double   dt;
    complex_t *x = NULL;
    complex_t *X = NULL;
    
    x = (complex_t*) malloc (NMAX * sizeof(complex_t));
    if(!x)
    {
        err = ERROR_MALLOC;
        goto exit_label;
    }
    X = (complex_t*) malloc (NMAX * sizeof(complex_t));
    if(!X)
    {
        err = ERROR_MALLOC;
        goto exit_label;
    }
    
    fft_t pfft = {0};       /* fft structure */
    random_t rnd = {0};     /* random structure   */
    
    err = random_init(&rnd, RAND_TYPE_MRG32K3A, NULL);
    if(err != RES_OK)
        goto exit_label;


    err = randn((double*)x, 2*NMAX, 0.0, 1.0, &rnd);
    if(err != RES_OK)
        goto exit_label;

    
    printf("--------------------\n");
    printf("FFT size      MFlops\n");
    printf("--------------------\n");
    for(j = 0; j < L; j++)
    {
        dlen[j] = len[j];
        t = clock();
        /* FFT WORK CYCLE*/
        for(i = 0; i < k; i++)
        {
            fft_cmplx(x, len[j], &pfft, X);
        }
        
        t = clock() - t;
        
        dt = ((double)t)/CLOCKS_PER_SEC/(double)k;
        
        mflops[j] = 5.0 * (double)len[j] * 
                    log2((double)len[j]) / dt / 1E6;
        
        printf ("%8d    %10.1f (k = %d)\n",  len[j], mflops[j], k);
        
        k /= SIZE_FACTOR;
    }

exit_label:
    fft_free(&pfft);
    free(x);
    free(X);
    return err;
}


int main(int argc, char* argv[])
{
    void* hdspl;            /* DSPL handle        */
    void* hplot;            /* GNUPLOT handles    */
    
    
    hdspl = dspl_load();    /* Load DSPL function */
    
    int len_r2[L] = {2, 4, 8, 16, 32, 64, 128, 256, 512, 1024, 2048, 4096,
                    8192, 16384, 32768, 65536, 131072, 262144, 524288, 1048576};
                    
    int len_nr[L] = {6, 9, 12, 15, 18, 24, 36, 80, 100, 108, 210, 504, 1000,
                    1960, 4725, 8000, 10368, 27000, 75600, 165375};

    int err;
    double   mflops[L] = {0};
    double dlen[L];
   
    printf("\n\nDouble precision complex 1D radix-2\n");
    fft_perf_cmplx(len_r2, dlen, mflops);
    writetxt(dlen, mflops, L, "dat/fft_cmplx_dspl_r2.txt");
    
    printf("\n\nDouble precision complex 1D non-powers of two\n");
    fft_perf_cmplx(len_nr, dlen, mflops);
    writetxt(dlen, mflops, L, "dat/fft_cmplx_dspl_nr.txt");
    

exit_label:
    /* free dspl handle */
    dspl_free(hdspl);
    return err;
}

