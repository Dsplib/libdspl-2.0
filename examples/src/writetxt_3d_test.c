#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
 
#include "dspl.h"
 
#define NX  20
#define NY  30
 
int main(int argc, char* argv[])
{
    void* hdspl;            /* DSPL handle        */
    void* hplot;            /* GNUPLOT handles    */
 
    hdspl = dspl_load();    /* Load DSPL function */
 
    double x[NX];
    double y[NY];
    double z[NX * NY];
    int n, m;
    int err;
 
    /* x vector from -2 to 2 */
    linspace(-2.0, 2.0, NX, DSPL_SYMMETRIC, x);
 
    /* y vector from -2.5 to 2.5 */
    linspace(-2.5, 2.5, NY, DSPL_SYMMETRIC, y);
 
    /* z(x,y) = x * exp(-x^2 - y^2) */
    for(n = 0; n < NX; n++)
    {
        for(m = 0; m < NY; m++)
        {
            z[n + m*NX] = x[n]*exp(-x[n]*x[n] - y[m]*y[m]);
        }
    }
 
    /* Save to files "dat/data3d.txt" */
    err = writetxt_3d(x, NX, y, NY, z, "dat/data3d.txt");
    printf("writetxt_3d error 0x%8x\n", err);
 
    /* plotting 3d surface by GNUPLOT */
    /* Create window 0 */
    err = gnuplot_create(argc, argv, 560, 480, "img/writetxt_3d.png", &hplot);
    printf("GNUPLOT err = %d\n", err);
    gnuplot_cmd(hplot, "set pm3d implicit at s");
    gnuplot_cmd(hplot, "set view 50, 340, 1, 1");
    gnuplot_cmd(hplot, "set xlabel 'x'");
    gnuplot_cmd(hplot, "set ylabel 'y'");
    gnuplot_cmd(hplot, "splot 'dat/data3d.txt' with lines");
    gnuplot_close(hplot);
 
    dspl_free(hdspl);      /* free dspl handle */
    return 0;
}