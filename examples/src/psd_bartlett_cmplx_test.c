#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "dspl.h"


#define N         8192
#define FS        1.0



void psd_plot(int argc, char* argv[], char* fn_png, char* fn_data)
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
    
    sprintf(cmd, "plot '%s' w l lt -1", fn_data);
    
    gnuplot_cmd(hplot, cmd);
    gnuplot_close(hplot);
}


int main(int argc, char* argv[])
{
    void* hdspl;  /* DSPL handle                  */

    hdspl = dspl_load();    /* Load DSPL function */
    random_t rnd = {0};     /* random structure   */
    
    complex_t *x   = NULL;
    double *psd = NULL;
    double *frq = NULL;

    int k, err;

    x   = (complex_t*)malloc(N*sizeof(complex_t));
    psd = (double*)malloc(N*sizeof(double));
    frq = (double*)malloc(N*sizeof(double));
    
    /* input signal fill as noise -20 dB/Hz power spectrum density */
    random_init(&rnd, RAND_TYPE_MRG32K3A, NULL);
    randn_cmplx(x, N, NULL, 0.1, &rnd);
    
    /* x[k] = 0.1 * cos(2 * pi * 0.2 * k) + cos(2 * pi * 0.26 * k) + noise */
    for(k = 0; k < N; k++)
    {
        RE(x[k]) +=     cos(M_2PI * 0.26 * (double)k) + 
                    0.1*cos(M_2PI*0.2* (double)k); 

        IM(x[k]) +=     sin(M_2PI * 0.26 * (double)k) + 
                    0.1*sin(M_2PI*0.2* (double)k); 
    }
    /* Twosided PSD logscale magnitude n = 8192, nfft = 8192. 
       This case is periodogram */
    err = psd_bartlett_cmplx(x, N, N, NULL, FS, 
                       DSPL_FLAG_LOGMAG | DSPL_FLAG_PSD_TWOSIDED, psd, frq);
    writetxt(frq, psd, N, "dat/psd_bartlett_cmplx_8192.txt");
    psd_plot(argc, argv, "img/psd_bartlett_cmplx_8192.png", 
                         "dat/psd_bartlett_cmplx_8192.txt");


    /* Twosided PSD with logscale magnitude n = 8192, nfft = 1024 */
    err = psd_bartlett_cmplx(x, N, N/8, NULL, FS, 
                       DSPL_FLAG_LOGMAG | DSPL_FLAG_PSD_TWOSIDED, psd, frq);
    writetxt(frq, psd, N/8, "dat/psd_bartlett_cmplx_1024.txt");
    psd_plot(argc, argv, "img/psd_bartlett_cmplx_1024.png", 
                         "dat/psd_bartlett_cmplx_1024.txt");
    

    /* Twosided PSD with logscale magnitude n = 8192, nfft = 128 */
    err = psd_bartlett_cmplx(x, N, N/64, NULL, FS, 
                       DSPL_FLAG_LOGMAG | DSPL_FLAG_PSD_TWOSIDED, psd, frq);
    writetxt(frq, psd, N/64, "dat/psd_bartlett_cmplx_128.txt");
    psd_plot(argc, argv, "img/psd_bartlett_cmplx_128.png", 
                         "dat/psd_bartlett_cmplx_128.txt");

    dspl_free(hdspl);      /* free dspl handle */
    if(x)
      free(x);
    if(psd)
      free(psd);
    if(frq)
      free(frq);
  return 0;
}


