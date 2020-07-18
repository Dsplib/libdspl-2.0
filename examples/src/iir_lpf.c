#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "dspl.h"

// Порядок фильтра
#define ORD 5

// размер векторов частотной характеристики фильтра
#define N   1000


int main(int argc, char* argv[])
{
    void* hdspl;  /* DSPL handle        */
    void* hplot;  /* GNUPLOT handle     */
    hdspl = dspl_load();       // Load DSPL function

    double a[ORD+1], b[ORD+1];  // коэффициенты H(s)
    double rs = 60.0;  // неравномерность в полосе пропускания 3дБ
    double rp = 1.0;
    // Частота (w), АЧХ (mag), ФЧХ (phi) и ГВЗ (tau)
    double w[N], mag[N], phi[N], tau[N];
    int k;

    // рассчитываем цифровой ФНЧ Чебышева 2 рода
    int res = iir(rp, rs, ORD, 0.3, 0.0,
                  DSPL_FILTER_CHEBY2 | DSPL_FILTER_LPF, b, a);
    if(res != RES_OK)
        printf("error code = 0x%8x\n", res);

    // печать коэффициентов фильтра
    for(k = 0; k < ORD+1; k++)
        printf("b[%2d] = %9.3f     a[%2d] = %9.3f\n", k, b[k], k, a[k]);

    // вектор нормированной частоты от 0 до pi
    linspace(0, M_PI, N , DSPL_PERIODIC, w);

    //частотные характеристика фильтра
    filter_freq_resp(b, a, ORD, w, N, DSPL_FLAG_LOGMAG|DSPL_FLAG_UNWRAP,
                    mag, phi, tau);

    for(k = 0; k < N; k++)
        w[k] /= M_PI;

    // Сохранить характеристики для построения графиков
    writetxt(w, mag, N, "dat/iir_lpf_mag.txt");
    writetxt(w, phi, N, "dat/iir_lpf_phi.txt");
    writetxt(w, tau, N, "dat/iir_lpf_tau.txt");

    /* plotting by GNUPLOT */
    gnuplot_create(argc, argv, 920, 260, "img/iir_lpf.png", &hplot);
    gnuplot_cmd(hplot, "unset key");
    gnuplot_cmd(hplot, "set grid");
    gnuplot_cmd(hplot, "set xlabel 'normalized frequency'");
    gnuplot_cmd(hplot, "set multiplot layout 1,3 rowsfirst");
    gnuplot_cmd(hplot, "set ylabel 'Magnitude, dB'");
    gnuplot_cmd(hplot, "set yrange [-100:5]");
    gnuplot_cmd(hplot, "plot 'dat/iir_lpf_mag.txt' with lines");
    gnuplot_cmd(hplot, "set ylabel 'Phase response, rad'");
    gnuplot_cmd(hplot, "unset yrange");
    gnuplot_cmd(hplot, "plot 'dat/iir_lpf_phi.txt' with lines");
    gnuplot_cmd(hplot, "set ylabel 'Groupdelay, samples'");
    gnuplot_cmd(hplot, "unset yrange");
    gnuplot_cmd(hplot, "plot 'dat/iir_lpf_tau.txt' with lines");
    gnuplot_cmd(hplot, "unset multiplot");
    gnuplot_close(hplot);

    dspl_free(hdspl);      // free dspl handle


  return 0;
}


