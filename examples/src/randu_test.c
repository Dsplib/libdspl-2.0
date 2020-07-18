#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "dspl.h"

#define N  10000

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

    /**************************************************************************/
    /* MRG32K3A random numbers generator                                      */
    /**************************************************************************/
    double seed_mrg32k3a = 1234.0;

    err = random_init(&rnd, RAND_TYPE_MRG32K3A, (void*)(&seed_mrg32k3a));
    if(err != RES_OK)
        goto exit_label;

    err = randu(x, N, &rnd);
    if(err != RES_OK)
        goto exit_label;

    err = randu(y, N, &rnd);
    if(err != RES_OK)
        goto exit_label;

    /* Save to files "dat/randu_mrg32k3a.txt" */
    writetxt(x, y, N, "dat/randu_mrg32k3a.txt");



    /**************************************************************************/
    /* MT19937 random numbers generator                                       */
    /**************************************************************************/
    unsigned long long seed_mt19937 = 1234353456;

    err = random_init(&rnd, RAND_TYPE_MT19937, (void*)(&seed_mt19937));
    if(err != RES_OK)
        goto exit_label;

    err = randu(x, N, &rnd);
    if(err != RES_OK)
        goto exit_label;

    err = randu(y, N, &rnd);
    if(err != RES_OK)
        goto exit_label;

    /* Save to files "dat/randu_mrg32k3a.txt" */
    writetxt(x, y, N, "dat/randu_mt19937.txt");


    /***************************************************************************/
    /* Standard C random numbers generator                                     */
    /***************************************************************************/
    err = randu(x, N, NULL);
    if(err != RES_OK)
        goto exit_label;

    err = randu(y, N, NULL);
    if(err != RES_OK)
        goto exit_label;

    /* Save to files "dat/randu_std.txt" */
    writetxt(x, y, N, "dat/randu_std.txt");



    /***************************************************************************/
    /* plotting by GNUPLOT                                                     */
    /***************************************************************************/
    /* Create window plot */
    err = gnuplot_create(argc, argv, 920, 320, "img/randu_test.png", &hplot);
    if(err != RES_OK)
        goto exit_label;

    gnuplot_cmd(hplot, "set multiplot layout 1,3 rowsfirst");
    gnuplot_cmd(hplot, "set xlabel 'x'");
    gnuplot_cmd(hplot, "set ylabel 'y'");
    gnuplot_cmd(hplot, "unset key");
    gnuplot_cmd(hplot, "set title 'MRG32K3A'");
    gnuplot_cmd(hplot, "plot 'dat/randu_mrg32k3a.txt' with points pointtype 0");
    gnuplot_cmd(hplot, "set title 'MT19937'");
    gnuplot_cmd(hplot, "plot 'dat/randu_mt19937.txt'  with points pointtype 0");
    gnuplot_cmd(hplot, "set title 'Standard C'");
    gnuplot_cmd(hplot, "plot 'dat/randu_std.txt'  with points pointtype 0");



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

