#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "dspl.h"

// Filter order (must be even for bandstop IIR)
#define ORD 14

// Frequency response vector size
#define N   5001


int main(int argc, char* argv[])
{
    void* hdspl;  /* DSPL handle        */
    void* hplot;  /* GNUPLOT handle     */
    hdspl = dspl_load();       // Load DSPL function

    double a[ORD+1], b[ORD+1];  // H(s) coefficients
    double rs = 60.0;  // Bandstop suppression equals 60 dB
    double rp = 1.0;   // Bandpass ripple equals 1 dB

    // Frequency (w), magnitude (mag), phase response (phi)
    // and group delay (tau)
    double w[N], mag[N], phi[N], tau[N];
    int k;

    // Calculate Chebyshev type 2 digital bandstop filter
    int res = iir(rp, rs, ORD, 0.3, 0.7,
                  DSPL_FILTER_ELLIP | DSPL_FILTER_BSTOP, b, a);
    if(res != RES_OK)
        printf("error code = 0x%8x\n", res);

    // Print coefficients
    for(k = 0; k < ORD+1; k++)
        printf("b[%2d] = %9.3f     a[%2d] = %9.3f\n", k, b[k], k, a[k]);

    // Normalized circular frequency from 0 to pi
    linspace(0, M_PI, N , DSPL_PERIODIC, w);

    // Filter magnitude, phase response and group delay
    filter_freq_resp(b, a, ORD, w, N, DSPL_FLAG_LOGMAG|DSPL_FLAG_UNWRAP,
                     mag, phi, tau);

    // Normalized frequency from 0 to 1.
    // w = 1 corresponds to Fs/2
    for(k = 0; k < N; k++)
        w[k] /= M_PI;

    // Save filter frequency response to the txt-files
    // for plotting by GNUPLOT
    writetxt(w, mag, N, "dat/iir_bstop_mag.txt");
    writetxt(w, phi, N, "dat/iir_bstop_phi.txt");
    writetxt(w, tau, N, "dat/iir_bstop_tau.txt");


    /* plotting by GNUPLOT */
    gnuplot_create(argc, argv, 920, 260, "img/iir_bstop.png", &hplot);
    gnuplot_cmd(hplot, "unset key");
    gnuplot_cmd(hplot, "set grid");
    gnuplot_cmd(hplot, "set xlabel 'normalized frequency'");
    gnuplot_cmd(hplot, "set multiplot layout 1,3 rowsfirst");
    gnuplot_cmd(hplot, "set ylabel 'Magnitude, dB'");
    gnuplot_cmd(hplot, "set yrange [-100:5]");
    gnuplot_cmd(hplot, "plot 'dat/iir_bstop_mag.txt' with lines");
    gnuplot_cmd(hplot, "set ylabel 'Phase response, rad'");
    gnuplot_cmd(hplot, "unset yrange");
    gnuplot_cmd(hplot, "plot 'dat/iir_bstop_phi.txt' with lines");
    gnuplot_cmd(hplot, "set ylabel 'Groupdelay, samples'");
    gnuplot_cmd(hplot, "unset yrange");
    gnuplot_cmd(hplot, "plot 'dat/iir_bstop_tau.txt' with lines");
    gnuplot_cmd(hplot, "unset multiplot");
    gnuplot_close(hplot);

    dspl_free(hdspl);      // free dspl handle

    // run GNUPLOT script
    return 0;
}


