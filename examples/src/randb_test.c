#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "dspl.h"

#define N  100

int main(int argc, char* argv[])
{
    void* hdspl;            /* DSPL handle        */
    void* hplot;            /* GNUPLOT handles    */
    random_t rnd = {0};     /* random structure   */
    hdspl = dspl_load();    /* Load DSPL function */

    double x[N];
    int err;

    err = random_init(&rnd, RAND_TYPE_MRG32K3A, NULL);
    if(err != RES_OK)
        goto exit_label;

    err = randb(x, N, &rnd);
    if(err != RES_OK)
        goto exit_label;

    writetxt(x, NULL, N, "dat/randb_mrg32k3a.txt");

    err = randb2(x, N, &rnd);
    if(err != RES_OK)
        goto exit_label;

    writetxt(x, NULL, N, "dat/randb2_mrg32k3a.txt");

    /**************************************************************************/
    /* plotting by GNUPLOT                                                    */
    /**************************************************************************/
    /* Create window plot */
    err = gnuplot_create(argc, argv, 820, 320, "img/randb_test.png", &hplot);
    if(err != RES_OK)
        goto exit_label;

    gnuplot_cmd(hplot, "set xlabel 'x'");
    gnuplot_cmd(hplot, "set ylabel 'y'");
    gnuplot_cmd(hplot, "unset key");
    gnuplot_cmd(hplot, "set yrange[-1.2:1.2]");
    gnuplot_cmd(hplot, "set title 'randb, randb2'");
    gnuplot_cmd(hplot, "plot 'dat/randb_mrg32k3a.txt' with lines, \\");
    gnuplot_cmd(hplot, "     'dat/randb2_mrg32k3a.txt' with lines");

exit_label:
    printf("Error code: %x\n", err);
    if(hplot)
        gnuplot_close(hplot);

    /* free dspl handle */
    dspl_free(hdspl);
    return 0;
}

