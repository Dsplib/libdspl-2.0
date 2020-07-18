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
    double Rs = 40;     /* Stopband suppression, dB */

    int res, k, nz, np;

    /* Zeros and poles vectors calculation */
    res = cheby2_ap_zp(ORD, Rs, z, &nz, p, &np);
    if(res != RES_OK)
        printf("error code = 0x%8x\n", res);

    /* print H(s) zeros values */
    printf("Chebyshev type 2 filter zeros: %d\n", nz);
    for(k = 0; k < nz; k++)
        printf("z[%2d] = %9.3f %+9.3f j\n", k, RE(z[k]), IM(z[k]));

    /* print H(s) poles values */
    printf("Chebyshev type 2 filter poles: %d\n", np);
    for(k = 0; k < np; k++)
        printf("p[%2d] = %9.3f %+9.3f j\n", k, RE(p[k]), IM(p[k]));

    /* Write complex poles to the file */
    writetxt_cmplx(z, nz, "dat/cheby2_ap_z.txt");
    writetxt_cmplx(p, np, "dat/cheby2_ap_p.txt");

    /* plotting by GNUPLOT */
    gnuplot_create(argc, argv, 440, 360, "img/cheby2_ap_zp_test.png", &hplot);
    gnuplot_cmd(hplot, "unset key");
    gnuplot_cmd(hplot, "set xzeroaxis");
    gnuplot_cmd(hplot, "set yzeroaxis");
    gnuplot_cmd(hplot, "set xlabel 'sigma'");
    gnuplot_cmd(hplot, "set ylabel 'jw'");
    gnuplot_cmd(hplot, "set size square");
    gnuplot_cmd(hplot, "set xrange [-3:3]");
    gnuplot_cmd(hplot, "set yrange [-3:3]");
    gnuplot_cmd(hplot, "plot 'dat/cheby2_ap_p.txt'  with points pt 2, \\");
    gnuplot_cmd(hplot, "     'dat/cheby2_ap_z.txt'  with points pt 6");
    gnuplot_close(hplot);

    /* free dspl handle */
    dspl_free(hdspl);

    return res;
}

