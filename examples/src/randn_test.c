#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "dspl.h"

#define N  50000

int main(int argc, char* argv[])
{
    void* hdspl;            /* DSPL handle        */
    void* hplot;            /* GNUPLOT handles    */
    random_t rnd = {0};     /* random structure   */
    hdspl = dspl_load();    /* Load DSPL function */

    double *x = NULL;
    double *y = NULL;
    int err;

    x = (double*) malloc(N * sizeof(double));
    y = (double*) malloc(N * sizeof(double));

    if(!x || !y)
    {
        err = ERROR_PTR;
        printf("malloc error!\n");
        goto exit_label;
    }


    err = random_init(&rnd, RAND_TYPE_MRG32K3A, NULL);
    if(err != RES_OK)
        goto exit_label;

    err = randn(x, N, 0.0, 1.0, &rnd);
    if(err != RES_OK)
        goto exit_label;

    err = randn(y, N, 0.0, 1.0, &rnd);
    if(err != RES_OK)
        goto exit_label;

    /* Save to files "dat/randu_mrg32k3a.txt" */
    writetxt(x, y, N, "dat/randn_mrg32k3a.txt");

    /**************************************************************************/
    /* plotting by GNUPLOT                                                    */
    /**************************************************************************/
    /* Create window plot */
    err = gnuplot_create(argc, argv, 420, 420, "img/randn_test.png", &hplot);
    if(err != RES_OK)
        goto exit_label;

    gnuplot_cmd(hplot, "set xlabel 'x'");
    gnuplot_cmd(hplot, "set ylabel 'y'");
    gnuplot_cmd(hplot, "unset key");
    gnuplot_cmd(hplot, "set xrange[-3:3]");
    gnuplot_cmd(hplot, "set yrange[-3:3]");
    gnuplot_cmd(hplot, "set title 'MRG32K3A'");
    gnuplot_cmd(hplot, "plot 'dat/randn_mrg32k3a.txt' with points pointtype 0");


exit_label:
    if(x)
        free(x);
    if(y)
        free(y);
    if(hplot)
        gnuplot_close(hplot);

    /* free dspl handle */
    dspl_free(hdspl);
    return 0;
}

