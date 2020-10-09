#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "dspl.h"


#define N     8192
#define FS    1.0


void psd_plot(int argc, char* argv[],   char* fn_png, 
                        char* fn_data0, char* fn_data1)
{
    char cmd[512] = {0};
    void* hplot;  /* GNUPLOT handle               */
    /* plotting by GNUPLOT */
    gnuplot_create(argc, argv, 680, 480, fn_png, &hplot);
    gnuplot_cmd(hplot, "unset key");
    gnuplot_cmd(hplot, "set grid");
    gnuplot_cmd(hplot, "set xlabel 'frequency, Hz'");
    gnuplot_cmd(hplot, "set ylabel 'Power Spectral Density, [dB/Hz]'");
    gnuplot_cmd(hplot, "set yrange [-60: 40]");
    
    sprintf(cmd, "plot '%s' w l lt -1,\\", fn_data0);
    gnuplot_cmd(hplot, cmd);
    
    memset(cmd, 0, 512);
    sprintf(cmd, "'%s' w l lt 2", fn_data1);
    gnuplot_cmd(hplot, cmd);
   
    gnuplot_close(hplot);
}


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

    x   = (double*)malloc(N*sizeof(double));
    psd = (double*)malloc(N*sizeof(double));
    frq = (double*)malloc(N*sizeof(double));
    
    /* input signal fill as noise -40 dB/Hz power spectrum density */
    random_init(&rnd, RAND_TYPE_MRG32K3A, NULL);
    randn(x, N, 0.0, 0.1, &rnd);
    
    /* x[k] = 0.1*cos(2*pi*0.1*k) + cos(2*pi*0.2*k) + noise */
    for(k = 0; k < N; k++)
        x[k] +=  cos(M_2PI * 0.26 * (double)k) + 0.1*cos(M_2PI*0.2* (double)k); 



    /* Twosided PSD with logscale magnitude 8192 points with rect window */
    err = psd_periodogram(x, N, DSPL_WIN_RECT, 0, NULL, FS,     
                           DSPL_FLAG_LOGMAG | DSPL_FLAG_PSD_TWOSIDED, psd, frq);
    writetxt(frq, psd, N, "dat/psd_periodogram_8192_win_rect.txt");


    /* Twosided PSD with logscale magnitude 8192 points with blackman window */
    err = psd_periodogram(x, N, DSPL_WIN_BLACKMAN, 0, NULL, FS,     
                           DSPL_FLAG_LOGMAG | DSPL_FLAG_PSD_TWOSIDED, psd, frq);
    writetxt(frq, psd, N, "dat/psd_periodogram_8192_win_blackman.txt");

    psd_plot(argc, argv, "img/psd_perodogram_8192.png", 
                         "dat/psd_periodogram_8192_win_rect.txt",
                         "dat/psd_periodogram_8192_win_blackman.txt");


    /* Twosided PSD with logscale magnitude 1024 points with rect window */
    err = psd_periodogram(x, N/8, DSPL_WIN_RECT, 0, NULL, FS,     
                           DSPL_FLAG_LOGMAG | DSPL_FLAG_PSD_TWOSIDED, psd, frq);
    writetxt(frq, psd, N/8, "dat/psd_periodogram_1024_win_rect.txt");
    
    /* Twosided PSD with logscale magnitude 1024 points with blackman window */
    err = psd_periodogram(x, N/8, DSPL_WIN_BLACKMAN, 0, NULL, FS,     
                           DSPL_FLAG_LOGMAG | DSPL_FLAG_PSD_TWOSIDED, psd, frq);
    writetxt(frq, psd, N/8, "dat/psd_periodogram_1024_win_blackman.txt");

    psd_plot(argc, argv, "img/psd_perodogram_1024.png", 
                         "dat/psd_periodogram_1024_win_rect.txt",
                         "dat/psd_periodogram_1024_win_blackman.txt");


    /* Twosided PSD with logscale magnitude 128 points with rect window */
    err = psd_periodogram(x, N/64, DSPL_WIN_RECT, 0, NULL, FS,     
                           DSPL_FLAG_LOGMAG | DSPL_FLAG_PSD_TWOSIDED, psd, frq);
    writetxt(frq, psd, N/64, "dat/psd_periodogram_128_win_rect.txt");
    
    /* Twosided PSD with logscale magnitude 128 points with blackman window */
    err = psd_periodogram(x, N/64, DSPL_WIN_BLACKMAN, 0, NULL, FS,     
                           DSPL_FLAG_LOGMAG | DSPL_FLAG_PSD_TWOSIDED, psd, frq);
    writetxt(frq, psd, N/64, "dat/psd_periodogram_128_win_blackman.txt");

    psd_plot(argc, argv, "img/psd_perodogram_128.png", 
                         "dat/psd_periodogram_128_win_rect.txt",
                         "dat/psd_periodogram_128_win_blackman.txt");



    dspl_free(hdspl);      // free dspl handle
    if(x)
      free(x);
    if(psd)
      free(psd);
    if(frq)
      free(frq);
  return 0;
}


