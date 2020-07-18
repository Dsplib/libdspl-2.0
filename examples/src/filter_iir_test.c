#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "dspl.h"

#define ORD 6
#define N   2000




int main(int argc, char* argv[])
{
    void* hdspl;  /* DSPL handle        */
    void* hplot;  /* GNUPLOT handle     */

    double b[ORD+1], a[ORD+1];
    double t[N], s[N], n[N], sf[N];
    random_t rnd;
    int k;
    int err;

    /* Load DSPL function  */
    hdspl = dspl_load();

    /* random generator init */
    random_init(&rnd, RAND_TYPE_MT19937, NULL);

    /* fill time vector         */
    linspace(0, N, N, DSPL_PERIODIC, t);

    /* generate noise        */
    randn(n, N, 0, 1.0, &rnd);

    /* input signal s = sin(2*pi*t) + n(t) */
    for(k = 0; k < N; k++)
        s[k] = sin(M_2PI*0.02*t[k]) + n[k];

    /* IIR filter coefficients calculation */
    iir(1.0, 70.0, ORD, 0.06, 0.0, DSPL_FILTER_ELLIP | DSPL_FILTER_LPF, b, a);

    /* input signal filtration */
    filter_iir(b, a, ORD, s, N, sf);

    /* save input signal and filter output to the txt-files */
    writetxt(t,s, N, "dat/s.txt");
    writetxt(t,sf,N, "dat/sf.txt");

    /* plotting by GNUPLOT */
    gnuplot_create(argc, argv, 820, 340, "img/filter_iir_test.png", &hplot);
    gnuplot_cmd(hplot, "unset key");
    gnuplot_cmd(hplot, "set grid");
    gnuplot_cmd(hplot, "set xlabel 'n'");
    gnuplot_cmd(hplot, "set ylabel 's(n)'");
    gnuplot_cmd(hplot, "set yrange [-3:3]");
    gnuplot_cmd(hplot, "set multiplot layout 2,1 rowsfirst");
    gnuplot_cmd(hplot, "plot 'dat/s.txt'  with lines");
    gnuplot_cmd(hplot, "set ylabel 's_f(n)'");
    gnuplot_cmd(hplot, "plot 'dat/sf.txt' with lines");
    gnuplot_cmd(hplot, "unset multiplot");
    gnuplot_close(hplot);

    /* free DSPL handle */
    dspl_free(hdspl);

    return err;
}
