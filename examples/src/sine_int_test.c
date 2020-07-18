#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "dspl.h"

#define N 400

int main(int argc, char* argv[])
{
    void* hdspl;           /* DSPL handle        */
    void* hplot;           /* GNUPLOT handle     */
    hdspl = dspl_load();   /* Load DSPL function */

    double x[N], y[N];

    linspace(-6*M_PI, 6*M_PI, N, DSPL_PERIODIC, x);

    sine_int(x, N, y);
    writetxt(x, y, N, "dat/dat0.txt");

    sinc(x, N, 1.0, y);
    writetxt(x, y, N, "dat/dat1.txt");

    /* plotting by GNUPLOT */
    gnuplot_create(argc, argv, 560, 280, "img/sine_int.png", &hplot);
    gnuplot_cmd(hplot, "set grid");
    gnuplot_cmd(hplot, "set xlabel 'x'");
    gnuplot_cmd(hplot, "set lmargin at screen 0.10");
    gnuplot_cmd(hplot, "set key left top");
    gnuplot_cmd(hplot, "set ylabel 'Si(x), sinc(x)'");
    gnuplot_cmd(hplot, "set yrange [-2:2]");
    gnuplot_cmd(hplot, "plot   'dat/dat0.txt' with lines title 'Si(x)', \\");
    gnuplot_cmd(hplot, "       'dat/dat1.txt' with lines title 'sinc(x)'");
    gnuplot_close(hplot);

    dspl_free(hdspl);      // free dspl handle

    return 0;
}