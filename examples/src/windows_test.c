#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
 
#include "dspl.h"
 
#define N  64
#define K  1024



void win_norm(double* w, int n)
{
    int k;
    double s = 0.0;
    for(k = 0; k < n; k++)
        s+=w[k];
    for(k = 0; k < n; k++)
        w[k] /= s;
}



void window_plot(double* w, int n, int k, int argc, char* argv[], 
                 double ymin, double ymax,
                 char* png_fn, char* time_fn, char* freq_fn, char* title)
{
    void* hplot;            /* GNUPLOT handles    */
    char str[128] = {0};
    double* t = NULL;
    double* W = NULL;
    double* f = NULL;
    double fs = (double)n;
    fft_t pfft = {0};
    
    t = (double*)malloc(n*sizeof(double));
    W = (double*)malloc(k*sizeof(double));
    f = (double*)malloc(k*sizeof(double));
    
    linspace(0, n, n, DSPL_PERIODIC, t);
    writetxt(t, w, n, time_fn);
    win_norm(w, n);
    fft_mag(w, k, &pfft, fs, DSPL_FLAG_LOGMAG | DSPL_FLAG_FFT_SHIFT, W, f);
    writetxt(f, W, k, freq_fn);


    /* plotting 3d surface by GNUPLOT */
    /* Create window 0 */
    gnuplot_create(argc, argv, 820, 360, png_fn, &hplot);
    
    gnuplot_cmd(hplot, "set grid");
    
    memset(str,0,128);
    sprintf(str, "set title '%s'", title);    
    gnuplot_cmd(hplot, str);
    
    gnuplot_cmd(hplot, "set multiplot layout 1,2 rowsfirst");
    gnuplot_cmd(hplot, "set xlabel 'n'");
    gnuplot_cmd(hplot, "set ylabel 'w(n)'");
    gnuplot_cmd(hplot, "unset key");
    
    gnuplot_cmd(hplot, "set xrange[0:63]");
    
    memset(str,0,128);
    sprintf(str, "set yrange[%f:%f]", ymin, ymax);
    gnuplot_cmd(hplot, str);
    
    memset(str,0,128);
    sprintf(str, "plot '%s' w i lc 1,'%s' w p pt 7 ps 0.5 lc 1", 
            time_fn, time_fn);
    gnuplot_cmd(hplot, str);
    
    gnuplot_cmd(hplot, "set xrange[-20:20]");
    gnuplot_cmd(hplot, "set yrange[-130:5]");
    memset(str,0,128);
    
    gnuplot_cmd(hplot, "set xlabel 'freq [DFT bins]'");
    gnuplot_cmd(hplot, "set ylabel 'W(freq), dB'");
    sprintf(str, "plot '%s' w l lc 2  ", freq_fn);
    gnuplot_cmd(hplot, str);
    gnuplot_close(hplot);
    
    if(t)
        free(t);
    if(W)
        free(W);
    if(f)
        free(f);
}



int main(int argc, char* argv[])
{
    void* hdspl;            /* DSPL handle        */
    
    hdspl = dspl_load();    /* Load DSPL function */
    double w[K] = {0};
    
    /* Rectangular window */
    window(w, N, DSPL_WIN_PERIODIC | DSPL_WIN_RECT, 0.0);
    window_plot(w, N, K, argc, argv, 0.0, 1.1, 
                "img/win_rect.png",
                "dat/win_rect_time.txt", 
                "dat/win_rect_freq.txt", 
                "Rectangular window");

    /* Bartlett window (triangular)*/
    window(w, N, DSPL_WIN_PERIODIC | DSPL_WIN_BARTLETT, 0.0);
    window_plot(w, N, K, argc, argv, 0.0, 1.1, 
                "img/win_bartlett.png", 
                "dat/win_bartlett_time.txt", 
                "dat/win_bartlett_freq.txt",
                "Bartlett window");

    /* Flat top window */
    window(w, N, DSPL_WIN_PERIODIC | DSPL_WIN_FLAT_TOP, 0.0);
    window_plot(w, N, K, argc, argv, -1.0, 5.0, 
                "img/win_flattop.png", 
                "dat/win_flattop_time.txt", 
                "dat/win_flattop_freq.txt",
                "Flat top window");

    /* Bartlett - Hann window*/
    window(w, N, DSPL_WIN_PERIODIC | DSPL_WIN_BARTLETT_HANN, 0.0);
    window_plot(w, N, K, argc, argv, 0.0, 1.1, 
                "img/win_bartletthann.png", 
                "dat/win_bartletthann_time.txt", 
                "dat/win_bartletthann_freq.txt",
                "Bartlett-Hann window");

    /* Blackman  window*/
    window(w, N, DSPL_WIN_PERIODIC | DSPL_WIN_BLACKMAN, 0.0);
    window_plot(w, N, K, argc, argv, 0.0, 1.1, 
                "img/win_blackman.png", 
                "dat/win_blackman_time.txt", 
                "dat/win_blackman_freq.txt",
                "Blackman window");

    /* Blackman - Harris window*/
    window(w, N, DSPL_WIN_PERIODIC | DSPL_WIN_BLACKMAN_HARRIS, 0.0);
    window_plot(w, N, K, argc, argv, 0.0, 1.1, 
                "img/win_blackmanharris.png", 
                "dat/win_blackmanharris_time.txt", 
                "dat/win_blackmanharris_freq.txt",
                "Blackman-Harris window");


    /* Blackman - Nuttall window*/
    window(w, N, DSPL_WIN_PERIODIC | DSPL_WIN_BLACKMAN_NUTTALL, 0.0);
    window_plot(w, N, K, argc, argv, 0.0, 1.1, 
                "img/win_blackmannuttall.png", 
                "dat/win_blackmannuttall_time.txt", 
                "dat/win_blackmannuttall_freq.txt",
                "Blackman-Nuttull window");

    /* Dolph-Chebyshev window (sidelobes level is -50 dB)*/
    window(w, N, DSPL_WIN_PERIODIC | DSPL_WIN_CHEBY, 50.0);
    window_plot(w, N, K, argc, argv, 0.0, 1.1, 
                "img/win_cheby50.png", 
                "dat/win_cheby50_time.txt", 
                "dat/win_cheby50_freq.txt",
                "Dolph-Chebyshev window (Rs = 50dB)");

    /* Dolph-Chebyshev window (sidelobes level is -80 dB)*/
    window(w, N, DSPL_WIN_PERIODIC | DSPL_WIN_CHEBY, 80.0);
    window_plot(w, N, K, argc, argv, 0.0, 1.1, 
                "img/win_cheby80.png", 
                "dat/win_cheby80_time.txt", 
                "dat/win_cheby80_freq.txt",
                "Dolph-Chebyshev window (Rs = 80dB)");

    /* Dolph-Chebyshev window (sidelobes level is -120 dB)*/
    window(w, N, DSPL_WIN_PERIODIC | DSPL_WIN_CHEBY, 120.0);
    window_plot(w, N, K, argc, argv, 0.0, 1.1, 
                "img/win_cheby120.png", 
                "dat/win_cheby120_time.txt", 
                "dat/win_cheby120_freq.txt", 
                "Dolph-Chebyshev window (Rs = 120dB)");

    /* Gaussian window (sigma = 0.5)*/
    window(w, N, DSPL_WIN_PERIODIC | DSPL_WIN_GAUSSIAN, 0.5);
    window_plot(w, N, K, argc, argv, 0.0, 1.1, 
                "img/win_gaussian0p5.png", 
                "dat/win_gaussian0p5_time.txt", 
                "dat/win_gaussian0p5_freq.txt", 
                "Gaussian window (sigma = 0.5)");

    /* Gaussian window (sigma = 0.3)*/
    window(w, N, DSPL_WIN_PERIODIC | DSPL_WIN_GAUSSIAN, 0.3);
    window_plot(w, N, K, argc, argv, 0.0, 1.1, 
                "img/win_gaussian0p3.png", 
                "dat/win_gaussian0p3_time.txt", 
                "dat/win_gaussian0p3_freq.txt", 
                "Gaussian window (sigma = 0.3)");

    /* Hamming window*/
    window(w, N, DSPL_WIN_PERIODIC | DSPL_WIN_HAMMING, 0.0);
    window_plot(w, N, K, argc, argv, 0.0, 1.1, 
                "img/win_hamming.png", 
                "dat/win_hamming_time.txt", 
                "dat/win_hamming_freq.txt", 
                "Hamming window");

    /* Hann window*/
    window(w, N, DSPL_WIN_PERIODIC | DSPL_WIN_HANN, 0.0);
    window_plot(w, N, K, argc, argv, 0.0, 1.1, 
                "img/win_hann.png", 
                "dat/win_hann_time.txt", 
                "dat/win_hann_freq.txt", 
                "Hann window");

    /* Lanczos window */
    window(w, N, DSPL_WIN_PERIODIC | DSPL_WIN_LANCZOS, 0.0);
    window_plot(w, N, K, argc, argv, 0.0, 1.1, 
                "img/win_lanczos.png", 
                "dat/win_lanczos_time.txt", 
                "dat/win_lanczos_freq.txt", 
                "Lanczos window");

    /* Nuttall window */
    window(w, N, DSPL_WIN_PERIODIC | DSPL_WIN_NUTTALL, 0.0);
    window_plot(w, N, K, argc, argv, 0.0, 1.1, 
                "img/win_nuttall.png", 
                "dat/win_nuttall_time.txt", 
                "dat/win_nuttall_freq.txt", 
                "Nuttall window");

    /* Cosine window */
    window(w, N, DSPL_WIN_PERIODIC | DSPL_WIN_COS, 0.0);
    window_plot(w, N, K, argc, argv, 0.0, 1.1, 
                "img/win_cos.png", 
                "dat/win_cos_time.txt", 
                "dat/win_cos_freq.txt", 
                "Cosine window");


    /* Kaiser window  pi * alpha = 4 */
    window(w, N, DSPL_WIN_PERIODIC | DSPL_WIN_KAISER, 4.0);
    window_plot(w, N, K, argc, argv, 0.0, 1.1, 
                "img/win_kaiser4p0.png", 
                "dat/win_kaiser4p0_time.txt", 
                "dat/win_kaiser4p0_freq.txt", 
                "Kaiser window (pi * alpha = 4)");
                
    /* Kaiser window  pi * alpha = 8 */
    window(w, N, DSPL_WIN_PERIODIC | DSPL_WIN_KAISER, 8.0);
    window_plot(w, N, K, argc, argv, 0.0, 1.1, 
                "img/win_kaiser8p0.png", 
                "dat/win_kaiser8p0_time.txt", 
                "dat/win_kaiser8p0_freq.txt", 
                "Kaiser window (pi * alpha = 8)");
                
    /* Kaiser window  pi * alpha = 12 */
    window(w, N, DSPL_WIN_PERIODIC | DSPL_WIN_KAISER, 12.0);
    window_plot(w, N, K, argc, argv, 0.0, 1.1, 
                "img/win_kaiser12p0.png", 
                "dat/win_kaiser12p0_time.txt", 
                "dat/win_kaiser12p0_freq.txt", 
                "Kaiser window (pi * alpha = 12)");


    dspl_free(hdspl);      /* free dspl handle */
    return 0;
}