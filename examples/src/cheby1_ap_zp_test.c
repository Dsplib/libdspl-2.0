#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "dspl.h"

/* Filter order */
#define ORD    7

int main(int argc, char* argv[])
{
    void* hdspl;  /* DSPL handle        */
    void* hplot;  /* GNUPLOT handle     */

    /* Load DSPL functions */
    hdspl = dspl_load();

    complex_t z[ORD+1]; /* H(s) zeros vector */
    complex_t p[ORD+1]; /* H(s) poles vector */
    double Rp = 0.5;    /* Magnitude ripple from 0 to 1 rad/s */

    int res, k, nz, np;

    /* Zeros and poles vectors calculation */
    res = cheby1_ap_zp(ORD, Rp, z, &nz, p, &np);
    if(res != RES_OK)
        printf("error code = 0x%8x\n", res);

    /* print H(s) zeros values */
    printf("Chebyshev type 1 filter zeros: %d\n", nz);
    for(k = 0; k < nz; k++)
        printf("z[%2d] = %9.3f %+9.3f j\n", k, RE(z[k]), IM(z[k]));

    /* print H(s) poles values */
    printf("Chebyshev type 1 filter poles: %d\n", np);
    for(k = 0; k < np; k++)
        printf("p[%2d] = %9.3f %+9.3f j\n", k, RE(p[k]), IM(p[k]));

    /* Write complex poles to the file */
    writetxt_cmplx(p, ORD, "dat/cheby1_ap_zp.txt");

    /* plotting by GNUPLOT */
    gnuplot_create(argc, argv, 440, 360, "img/cheby1_ap_zp_test.png", &hplot);
    gnuplot_cmd(hplot, "unset key");
    gnuplot_cmd(hplot, "set xzeroaxis");
    gnuplot_cmd(hplot, "set yzeroaxis");
    gnuplot_cmd(hplot, "set xlabel 'sigma'");
    gnuplot_cmd(hplot, "set ylabel 'jw'");
    gnuplot_cmd(hplot, "set size square");
    gnuplot_cmd(hplot, "set xrange [-1.5:1.5]");
    gnuplot_cmd(hplot, "set yrange [-1.5:1.5]");
    gnuplot_cmd(hplot, "plot 'dat/cheby1_ap_zp.txt'  with points pt 2");
    gnuplot_close(hplot);

    /* free dspl handle */
    dspl_free(hdspl);

    return res;
}

