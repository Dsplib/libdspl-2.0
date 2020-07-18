#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "dspl.h"

#define N  100

int main(int argc, char* argv[])
{
    void* hdspl;            /* DSPL handle        */
    void* hplot[3];         /* GNUPLOT handles    */

    hdspl = dspl_load();    /* Load DSPL function */

    double x[N];
    double s[N];  /* s(x) = sin(x) */
    double c[N];  /* c(x) = cos(x) */
    int n;
    int err;

    /* x vector from -4*pi to 4*pi */
    linspace(-4.0 * M_PI, 4 * M_PI, N , DSPL_SYMMETRIC, x);
    for(n = 0; n < N; n++)
    {
        s[n] = sin(x[n]); /* s(x) = sin(x) */
        c[n] = cos(x[n]); /* c(x) = cos(x) */
    }

    /* Save to files "dat/sine.txt" and "dat/cosine.txt" */
    writetxt(x, s, N, "dat/sine.txt");
    writetxt(x, c, N, "dat/cosine.txt");


    /* plotting by GNUPLOT */
    /* Create window 0 */
    err = gnuplot_create(argc, argv, 560, 280, "img/gnuplot_script_sin.png",
                         hplot);

    printf("GNUPLOT err = %d\n", err);
    gnuplot_cmd(hplot[0], "set grid");
    gnuplot_cmd(hplot[0], "set xlabel 'x'");
    gnuplot_cmd(hplot[0], "set ylabel 'sin(x)'");
    gnuplot_cmd(hplot[0], "plot 'dat/sine.txt' with lines title 'sin(x)'");
    gnuplot_close(hplot[0]);

    /* Create window 1 */
    err = gnuplot_create(argc, argv, 560, 280, "img/gnuplot_script_cos.png",
                         hplot+1);

    printf("GNUPLOT err = %d\n", err);
    gnuplot_cmd(hplot[1], "set grid");
    gnuplot_cmd(hplot[1], "set xlabel 'x'");
    gnuplot_cmd(hplot[1], "set ylabel 'cos(x)'");
    gnuplot_cmd(hplot[1], "plot 'dat/cosine.txt' with lines title 'cos(x)'");
    gnuplot_close(hplot[1]);

    /* Create window 2 */
    err = gnuplot_create(argc, argv, 560, 280, "img/gnuplot_script_sincos.png",
                         hplot+2);

    printf("GNUPLOT err = %d\n", err);
    gnuplot_cmd(hplot[2], "set grid");
    gnuplot_cmd(hplot[2], "set xlabel 'x'");
    gnuplot_cmd(hplot[2], "set ylabel 'sin(x), cos(x)'");
    gnuplot_cmd(hplot[2], "plot   'dat/sine.txt'   with lines title 'sin(x)', \\");
    gnuplot_cmd(hplot[2], "       'dat/cosine.txt' with lines title 'cos(x)");
    gnuplot_close(hplot[2]);

    dspl_free(hdspl);      /* free dspl handle */
    return 0;
}

