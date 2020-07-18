#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "dspl.h"

#define N  1000

int main(int argc, char* argv[])
{
    void* hdspl;            /* DSPL handle        */
    void* hplot;            /* GNUPLOT handles    */
    random_t rnd = {0};     /* random structure   */
    hdspl = dspl_load();    /* Load DSPL function */

    int x[N];
    int err;

    err = random_init(&rnd, RAND_TYPE_MRG32K3A, NULL);
    if(err != RES_OK)
        goto exit_label;

    err = randi(x, N, -4, 3, &rnd);
    if(err != RES_OK)
        goto exit_label;

    /* Save to file */
    err = writetxt_int(x, NULL, N, "dat/randi_mrg32k3a.txt");
    if(err != RES_OK)
        goto exit_label;

    /**************************************************************************/
    /* plotting by GNUPLOT                                                    */
    /**************************************************************************/
    /* Create window plot */
    err = gnuplot_create(argc, argv, 820, 320, "img/randi_test.png", &hplot);
    if(err != RES_OK)
        goto exit_label;

    gnuplot_cmd(hplot, "set xlabel 'n'");
    gnuplot_cmd(hplot, "set ylabel 'x(n)'");
    gnuplot_cmd(hplot, "unset key");
    gnuplot_cmd(hplot, "set yrange[-5:4]");
    gnuplot_cmd(hplot, "set title 'randi example'");
    gnuplot_cmd(hplot, "plot 'dat/randi_mrg32k3a.txt' with impulses");

exit_label:
    printf("Error code = 0x%.8x\n", err);

    if(hplot)
        gnuplot_close(hplot);

    /* free dspl handle */
    dspl_free(hdspl);
    return 0;
}

