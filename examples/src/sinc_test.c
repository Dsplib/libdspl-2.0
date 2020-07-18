#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "dspl.h"

#define N   1000


int main(int argc, char* argv[])
{
    void* hdspl;           /* DSPL handle        */
    void* hplot;           /* GNUPLOT handle     */
    hdspl = dspl_load();   /* Load DSPL function */

    double x[N], y[N];

    linspace(-10.0, 10.0, N , DSPL_SYMMETRIC, x);
    sinc(x, N, 1.0, y);
    writetxt(x, y, N, "dat/sinc_test_1.0.txt");

    sinc(x, N, M_PI, y);
    writetxt(x, y, N, "dat/sinc_test_pi.txt");

    sinc(x, N, 2.0*M_PI, y);
    writetxt(x, y, N, "dat/sinc_test_2pi.txt");


    /* plotting by GNUPLOT */
    gnuplot_create(argc, argv, 560, 280, "img/sinc_test.png", &hplot);
    gnuplot_cmd(hplot, "set grid");
    gnuplot_cmd(hplot, "set xlabel 'x'");
    gnuplot_cmd(hplot, "set ylabel 'sinc(x,a)'");
    gnuplot_cmd(hplot, "set yrange [-0.25:1.1]");
    gnuplot_cmd(hplot, "plot   'dat/sinc_test_1.0.txt' with lines  title 'a = 1.0'");
    gnuplot_cmd(hplot, "replot 'dat/sinc_test_pi.txt'  with lines  title 'a = pi'");
    gnuplot_cmd(hplot, "replot 'dat/sinc_test_2pi.txt' with lines  title 'a = 2pi'");
    gnuplot_close(hplot);

    dspl_free(hdspl);      // free dspl handle

    return 0;
}

