#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "dspl.h"


#define N   500


int main(int argc, char* argv[])
{
    void* hdspl;           /* DSPL handle        */
    void* hplot;           /* GNUPLOT handle     */
    hdspl = dspl_load();   /* Load DSPL function */

    double u[N];
    double y[N];
    
    /* fill u*K(k) vector from -4 to 4. 
    We will have 2 periods sn(uK(k), k) */
    linspace(-4.0, 4.0, N, DSPL_PERIODIC, u);
    
    /* sn(uK(0), 0) */
    ellip_cd(u, N, 0.0, y);
    writetxt(u,y,N,"dat/ellip_cd_0p00.txt");
    
    /* sn(uK(0.9), 0.9) */
    ellip_cd(u, N, 0.9, y);
    writetxt(u,y,N,"dat/ellip_cd_0p90.txt");
    
    /* sn(uK(0.99), 0.99) */
    ellip_cd(u, N, 0.99, y);
    writetxt(u,y,N,"dat/ellip_cd_0p99.txt");
    
    /* plotting by GNUPLOT                                                    */
    gnuplot_create(argc, argv, 560, 360, "img/ellip_cd.png", &hplot);
    gnuplot_cmd(hplot, "set grid");
    gnuplot_cmd(hplot, "set xlabel 'u'");
    gnuplot_cmd(hplot, "set ylabel 'y = cd(u, k)'");
    gnuplot_cmd(hplot, "set yrange [-1.2:1.2]");
    gnuplot_cmd(hplot, "plot 'dat/ellip_cd_0p00.txt' w l,\\");
    gnuplot_cmd(hplot, "     'dat/ellip_cd_0p90.txt' w l,\\");
    gnuplot_cmd(hplot, "     'dat/ellip_cd_0p99.txt' w l");
    gnuplot_close(hplot);

    /* free dspl handle */
    dspl_free(hdspl);

    return RES_OK;
}


