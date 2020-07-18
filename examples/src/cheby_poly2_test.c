#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "dspl.h"

#define N 250

int main(int argc, char* argv[])
{
    void* hdspl;  /* DSPL handle        */
    void* hplot;  /* GNUPLOT handle     */


    double x[N], y[N];
    int ord;
    char fn[64];
    hdspl = dspl_load();   /* Load DSPL function */
    linspace(-1.0, 1.0, N, DSPL_SYMMETRIC, x);
    for(ord = 1; ord < 5; ord++)
    {
        cheby_poly2(x, N, ord, y);
        sprintf(fn, "dat/cheby_poly2_ord%d.txt", ord);
        writetxt(x,y,N,fn);
    }

    /* plotting by GNUPLOT */
    gnuplot_create(argc, argv, 560, 380, "img/cheby_poly2.png", &hplot);
    gnuplot_cmd(hplot, "set grid");
    gnuplot_cmd(hplot, "set key left top");
    gnuplot_cmd(hplot, "set xlabel 'x'");
    gnuplot_cmd(hplot, "set ylabel 'U_N (x)'");
    gnuplot_cmd(hplot, "set yrange [-3.5:3.5]");
    gnuplot_cmd(hplot, "plot 'dat/cheby_poly2_ord1.txt' with lines, \\");
    gnuplot_cmd(hplot, "     'dat/cheby_poly2_ord2.txt' with lines, \\");
    gnuplot_cmd(hplot, "     'dat/cheby_poly2_ord3.txt' with lines, \\");
    gnuplot_cmd(hplot, "     'dat/cheby_poly2_ord4.txt' with lines");
    gnuplot_close(hplot);

    dspl_free(hdspl);      /* free dspl handle */

    return 0;
}