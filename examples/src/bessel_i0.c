#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "dspl.h"

#define N 50



/******************************************************************************
 * Main function
 ******************************************************************************/
int main(int argc, char* argv[])
{

    void* hdspl;    /* DSPL handle */
    void* hplot;    /* GNUPLOT handle */
    double x[N], y[N];

    /* Load DSPL function */
    hdspl = dspl_load();

    /* x in [0, 3] */
    linspace(0.0, 3.0, N, DSPL_SYMMETRIC, x);

    /* Bessel I0(x) function */
    bessel_i0(x, N, y);

    /* Write calculated values to the dat/dat0.txt file */
    writetxt(x, y, N, "dat/dat0.txt");


    /* plotting by GNUPLOT */
    gnuplot_create(argc, argv, 560, 380, "img/bessel_i0.png", &hplot);
    gnuplot_cmd(hplot, "set grid");
    gnuplot_cmd(hplot, "set xlabel 'x'");
    gnuplot_cmd(hplot, "set key left top");
    gnuplot_cmd(hplot, "set ylabel 'I_0(x)'");
    gnuplot_cmd(hplot, "set yrange [0:5]");
    gnuplot_cmd(hplot, "plot 'dat/dat0.txt' with lines");
    gnuplot_close(hplot);


    /* free dspl handle */
    dspl_free(hdspl);

    return 0;
}

