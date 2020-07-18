#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "dspl.h"

#define N      1000
#define ORD    4

int main(int argc, char* argv[])
{
    void* hdspl;    /* DSPL handle */
    void* hplot;    /* GNUPLOT handle */

    double w[N], h[N];
    complex_t hz[N];
    double bs[ORD+1], as[ORD+1];
    double bz[ORD+1], az[ORD+1];
    int err, k;

    hdspl = dspl_load();     /* Load DSPL function */

    /* normalized analog lowpass filter Chebyshev type 1 */
    err = cheby1_ap(1.0, ORD, bs, as);
    if(err != RES_OK)
    {
        printf("cheby1_ap error code %d\n", err);
        return err;
    }

    /* Bilinear transform */
    err = bilinear(bs, as, ORD, bz, az);
    if(err != RES_OK)
    {
        printf("bilinear error code %d\n", err);
        return err;
    }

    /* Print coefficients */
    for(k = 0; k < ORD+1; k++)
        printf("bz[%d] = %7.3f  az[%d] = %7.3f\n", k, bz[k], k, az[k]);

    /* Digital filter magnitude    */
    linspace(0, M_PI, N, DSPL_PERIODIC, w);
    freqz(bz, az, ORD, w, N, hz);
    for(k = 0; k < N; k++)
    {
        h[k] = 10.0 * log10 (ABSSQR(hz[k])); /* Logarithmic scale         */
        w[k] /= M_PI;                        /* frequency from 0 to 1     */
    }

    writetxt(w,h,N,"dat/bilinear.txt");

    /* plotting by GNUPLOT */
    gnuplot_create(argc, argv, 560, 380, "img/bilinear.png", &hplot);
    gnuplot_cmd(hplot, "set grid");
    gnuplot_cmd(hplot, "unset key");
    gnuplot_cmd(hplot, "set xlabel 'normalized frequency'");
    gnuplot_cmd(hplot, "set ylabel 'Magnitude, dB'");
    gnuplot_cmd(hplot, "set yrange [-80:5]");
    gnuplot_cmd(hplot, "plot 'dat/bilinear.txt' with lines");
    gnuplot_close(hplot);

     /* free dspl handle */
    dspl_free(hdspl);

    return err;
}

