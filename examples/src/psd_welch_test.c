#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "dspl.h"


#define N     8192
#define NFFT  512

#define FS   0.01



int main(int argc, char* argv[])
{
    void* hdspl;  /* DSPL handle        */
    void* hplot;  /* GNUPLOT handle     */
    hdspl = dspl_load();       // Load DSPL function
    random_t rnd = {0};     /* random structure   */
    
    double *x = NULL;
    double *psd = NULL;
    double *frq = NULL;
    
    int k, err;

    x = (double*)malloc(N*sizeof(double));
    psd = (double*)malloc(N*sizeof(double));
    frq = (double*)malloc(N*sizeof(double));
    
    /* input signal fill as noise -40 dB/Hz power spectrum density */
    random_init(&rnd, RAND_TYPE_MRG32K3A, NULL);
    randn(x, N, 0.0, 0.1, &rnd);
    
    /* x[k] = 0.1*cos(2*pi*0.1*k) + cos(2*pi*0.2*k) + noise */
    for(k = 0; k < N; k++)
        x[k] +=  cos(M_2PI * 0.26 * (double)k) + 0.1*cos(M_2PI*0.2* (double)k); 
               
    /* Twosided PSD with logscale magnitude */
    err = psd_welch(x, N, DSPL_WIN_BLACKMAN , 0, NFFT, NFFT, NULL, FS, 
                    DSPL_FLAG_LOGMAG | DSPL_FLAG_PSD_TWOSIDED, psd, frq);

    printf("error: 0x%.8x\n", err);

    // Save PSD to the "dat/psd_welch.txt" file
    writetxt(frq, psd, NFFT, "dat/psd_welch.txt");
    
    
    /* Twosided PSD with logscale magnitude */
    err = psd_periodogram(x, N, DSPL_WIN_BLACKMAN, 0, NULL, FS,     
                           DSPL_FLAG_LOGMAG | DSPL_FLAG_PSD_TWOSIDED, psd, frq);
    printf("error: 0x%.8x\n", err);
    writetxt(frq, psd, N, "dat/psd_periodogram.txt");
    
    
       /* Twosided PSD with logscale magnitude */
    err = psd_bartlett(x, N, NFFT, NULL, FS,     
                       DSPL_FLAG_LOGMAG | DSPL_FLAG_PSD_TWOSIDED, psd, frq);
    printf("error: 0x%.8x\n", err);

    // Save PSD to the "dat/psd_welch.txt" file
    writetxt(frq, psd, NFFT, "dat/psd_barlett.txt");

    /* plotting by GNUPLOT */
    gnuplot_create(argc, argv, 560, 420, "img/psd_welch.png", &hplot);
    gnuplot_cmd(hplot, "unset key");
    gnuplot_cmd(hplot, "set grid");
    gnuplot_cmd(hplot, "set xlabel 'frequency'");
    gnuplot_cmd(hplot, "set ylabel 'PSD, dB'");

    gnuplot_cmd(hplot, "plot 'dat/psd_periodogram.txt' w l lt 3,\\");
    gnuplot_cmd(hplot, "     'dat/psd_barlett.txt' w l lt 1,\\");
    gnuplot_cmd(hplot, "     'dat/psd_welch.txt' w l lt -1");

    gnuplot_close(hplot);

    dspl_free(hdspl);      // free dspl handle
    if(x)
      free(x);

  return 0;
}


