#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "dspl.h"


/*filter order */
#define ORD 64

#define W0  0.3
#define W1  0.7

/* Frequency response vector size */
#define N  1024


/*******************************************************************************
 * function calculates filter frequency response and save magnitude to
 * the text file.
 * params: b - pointer to the transfer fuction H(z) numerator vector
 *         a - pointer to the transfer fuction H(z) denominator vector
 *         ord - filter order
 *         n - number of magnitude vector size
 *         fn - file name
 ******************************************************************************/
void freq_resp_write2txt(double* b, int ord, int n, char* fn)
{
    double *w = NULL, *mag = NULL;
    int k;

    w   = (double*)malloc(n*sizeof(double));
    mag = (double*)malloc(n*sizeof(double));

    /* Normalized frequency from 0 to pi */
    linspace(0, M_PI, n , DSPL_PERIODIC, w);

    /* Magnitude (dB) calculation */
    filter_freq_resp(b, NULL, ord, w, n, DSPL_FLAG_LOGMAG, mag, NULL, NULL);

    /* Frequency normaliztion from 0 to 1 */
    for(k = 0; k < N; k++)
        w[k] /= M_PI;

    /* Save magnitude to the txt file */
    writetxt(w, mag, n, fn);

    free(w);
    free(mag);
}


/*******************************************************************************
 * Main program
 ******************************************************************************/
int main(int argc, char* argv[])
{
    void* hdspl;           /* DSPL handle         */
    void* hplot;           /* GNUPLOT handle      */
    hdspl = dspl_load();   /* Load DSPL functions */

    /* FIR filter coeff. vectors */
    double h[ORD+1];

    /*------------------------------------------------------------------------*/
    fir_linphase(ORD, W0, W1, DSPL_FILTER_LPF, DSPL_WIN_RECT, 0, h);
    freq_resp_write2txt(h, ORD, N, "dat/fir_lpf_rect.txt");

    fir_linphase(ORD, W0, W1, DSPL_FILTER_LPF, DSPL_WIN_HAMMING, 0, h);
    freq_resp_write2txt(h, ORD, N, "dat/fir_lpf_hamming.txt");

    fir_linphase(ORD, W0, W1, DSPL_FILTER_LPF, DSPL_WIN_BLACKMAN, 0, h);
    freq_resp_write2txt(h, ORD, N, "dat/fir_lpf_blackman.txt");

    fir_linphase(ORD, W0, W1, DSPL_FILTER_LPF, DSPL_WIN_BLACKMAN_HARRIS, 0, h);
    freq_resp_write2txt(h, ORD, N, "dat/fir_lpf_blackman_harris.txt");

    /*------------------------------------------------------------------------*/

    fir_linphase(ORD, W0, W1, DSPL_FILTER_HPF, DSPL_WIN_RECT, 0, h);
    freq_resp_write2txt(h, ORD, N, "dat/fir_hpf_rect.txt");

    fir_linphase(ORD, W0, W1, DSPL_FILTER_HPF, DSPL_WIN_HAMMING, 0, h);
    freq_resp_write2txt(h, ORD, N, "dat/fir_hpf_hamming.txt");

    fir_linphase(ORD, W0, W1, DSPL_FILTER_HPF, DSPL_WIN_BLACKMAN, 0, h);
    freq_resp_write2txt(h, ORD, N, "dat/fir_hpf_blackman.txt");

    fir_linphase(ORD, W0, W1, DSPL_FILTER_HPF, DSPL_WIN_BLACKMAN_HARRIS, 0, h);
    freq_resp_write2txt(h, ORD, N, "dat/fir_hpf_blackman_harris.txt");

    /*------------------------------------------------------------------------*/

    fir_linphase(ORD, W0, W1, DSPL_FILTER_BPASS, DSPL_WIN_RECT, 0, h);
    freq_resp_write2txt(h, ORD, N, "dat/fir_bpass_rect.txt");

    fir_linphase(ORD, W0, W1, DSPL_FILTER_BPASS, DSPL_WIN_HAMMING, 0, h);
    freq_resp_write2txt(h, ORD, N, "dat/fir_bpass_hamming.txt");

    fir_linphase(ORD, W0, W1, DSPL_FILTER_BPASS, DSPL_WIN_BLACKMAN, 0, h);
    freq_resp_write2txt(h, ORD, N, "dat/fir_bpass_blackman.txt");

    fir_linphase(ORD, W0, W1, DSPL_FILTER_BPASS, DSPL_WIN_BLACKMAN_HARRIS, 0, h);
    freq_resp_write2txt(h, ORD, N, "dat/fir_bpass_blackman_harris.txt");

    /*------------------------------------------------------------------------*/

    fir_linphase(ORD, W0, W1, DSPL_FILTER_BSTOP, DSPL_WIN_RECT, 0, h);
    freq_resp_write2txt(h, ORD, N, "dat/fir_bstop_rect.txt");

    fir_linphase(ORD, W0, W1, DSPL_FILTER_BSTOP, DSPL_WIN_HAMMING, 0, h);
    freq_resp_write2txt(h, ORD, N, "dat/fir_bstop_hamming.txt");

    fir_linphase(ORD, W0, W1, DSPL_FILTER_BSTOP, DSPL_WIN_BLACKMAN, 0, h);
    freq_resp_write2txt(h, ORD, N, "dat/fir_bstop_blackman.txt");

    fir_linphase(ORD, W0, W1, DSPL_FILTER_BSTOP, DSPL_WIN_BLACKMAN_HARRIS, 0, h);
    freq_resp_write2txt(h, ORD, N, "dat/fir_bstop_blackman_harris.txt");

    /*------------------------------------------------------------------------*/

    /* plotting by GNUPLOT */
    gnuplot_create(argc, argv, 920, 840, "img/fir_linphase_test.png", &hplot);
    gnuplot_cmd(hplot, "unset key");
    gnuplot_cmd(hplot, "set grid");
    gnuplot_cmd(hplot, "unset xlabel");

    gnuplot_cmd(hplot, "set yrange [-130:5]");
    gnuplot_cmd(hplot, "set xtics 0,1");
    gnuplot_cmd(hplot, "set xtics add ('0.3' 0.3)");
    gnuplot_cmd(hplot, "set xtics add ('0.7' 0.7)");
    gnuplot_cmd(hplot, "set xtics add ('1' 1)");
    gnuplot_cmd(hplot, "set multiplot layout 4,4 rowsfirst");

    gnuplot_cmd(hplot, "set ylabel 'Magnitude, dB'");
    gnuplot_cmd(hplot, "set title 'Rect win'");
    gnuplot_cmd(hplot, "plot 'dat/fir_lpf_rect.txt' with lines");
    gnuplot_cmd(hplot, "unset ylabel");
    gnuplot_cmd(hplot, "set title 'Hamming win'");
    gnuplot_cmd(hplot, "plot 'dat/fir_lpf_hamming.txt' with lines");
    gnuplot_cmd(hplot, "set title 'Blackman win'");
    gnuplot_cmd(hplot, "plot 'dat/fir_lpf_blackman.txt' with lines");
    gnuplot_cmd(hplot, "set title 'Blackman-Harris win'");
    gnuplot_cmd(hplot, "plot 'dat/fir_lpf_blackman_harris.txt' with lines");
    gnuplot_cmd(hplot, "unset title");

    gnuplot_cmd(hplot, "set ylabel 'Magnitude, dB'");
    gnuplot_cmd(hplot, "plot 'dat/fir_hpf_rect.txt' with lines");
    gnuplot_cmd(hplot, "unset ylabel");
    gnuplot_cmd(hplot, "plot 'dat/fir_hpf_hamming.txt' with lines");
    gnuplot_cmd(hplot, "plot 'dat/fir_hpf_blackman.txt' with lines");
    gnuplot_cmd(hplot, "plot 'dat/fir_hpf_blackman_harris.txt' with lines");

    gnuplot_cmd(hplot, "set ylabel 'Magnitude, dB'");
    gnuplot_cmd(hplot, "plot 'dat/fir_bpass_rect.txt' with lines");
    gnuplot_cmd(hplot, "unset ylabel");
    gnuplot_cmd(hplot, "plot 'dat/fir_bpass_hamming.txt' with lines");
    gnuplot_cmd(hplot, "plot 'dat/fir_bpass_blackman.txt' with lines");
    gnuplot_cmd(hplot, "plot 'dat/fir_bpass_blackman_harris.txt' with lines");

    gnuplot_cmd(hplot, "set ylabel 'Magnitude, dB'");
    gnuplot_cmd(hplot, "set xlabel 'normalized frequency'");
    gnuplot_cmd(hplot, "plot 'dat/fir_bstop_rect.txt' with lines");
    gnuplot_cmd(hplot, "unset ylabel");
    gnuplot_cmd(hplot, "plot 'dat/fir_bstop_hamming.txt' with lines");
    gnuplot_cmd(hplot, "plot 'dat/fir_bstop_blackman.txt' with lines");
    gnuplot_cmd(hplot, "plot 'dat/fir_bstop_blackman_harris.txt' with lines");

    gnuplot_cmd(hplot, "unset multiplot");
    gnuplot_close(hplot);

    /* free dspl handle */
    dspl_free(hdspl);

    return 0;
}



