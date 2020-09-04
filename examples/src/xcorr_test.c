#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "dspl.h"

#define N 500
#define P 490

int main(int argc, char* argv[])
{
    void* hdspl;            /* DSPL handle        */
    void* hplot;            /* GNUPLOT handles    */
    random_t rnd = {0};     /* random structure   */
    hdspl = dspl_load();    /* Load DSPL function */

    double x[N];
    double z[2*P+1];
    double t[2*P+1];
    
    int err;
    
    err = random_init(&rnd, RAND_TYPE_MRG32K3A, NULL);
    if(err != RES_OK)
        goto exit_label;

    err = randn(x, N, 0.0, 1.0, &rnd);
    if(err != RES_OK)
        goto exit_label;

    err = xcorr(x, N, x, N, DSPL_XCORR_UNBIASED, P, z, t);

    /* Save to files "dat/randu_mrg32k3a.txt" */
    writetxt(t, z, 2*P+1, "dat/xcorr_unbiased.txt");
    
    err = xcorr(x, N, x, N, DSPL_XCORR_BIASED, P, z, t);

    /* Save to files "dat/randu_mrg32k3a.txt" */
    writetxt(t, z, 2*P+1, "dat/xcorr_biased.txt");

    /**************************************************************************/
    /* plotting by GNUPLOT                                                    */
    /**************************************************************************/
    /* Create window plot */
    err = gnuplot_create(argc, argv, 420, 420, "img/xcorr_test.png", &hplot);
    if(err != RES_OK)
        goto exit_label;

    gnuplot_cmd(hplot, "set xlabel 't'");
    gnuplot_cmd(hplot, "set ylabel 'z'");
    gnuplot_cmd(hplot, "unset key");
    gnuplot_cmd(hplot, "set title 'xcorr'");
    gnuplot_cmd(hplot, "plot 'dat/xcorr_unbiased.txt' with lines, 'dat/xcorr_biased.txt' with lines");


exit_label:
    if(hplot)
        gnuplot_close(hplot);

    /* free dspl handle */
    dspl_free(hdspl);
    return 0;
}

