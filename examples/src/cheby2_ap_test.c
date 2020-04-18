#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "dspl.h"

/* Порядок фильтра */
#define ORD 4

/* размер векторов частотной характеристики фильтра */
#define N   1000

int main(int argc, char* argv[])
{
  void* hdspl;           /* DSPL handle                                     */
  void* hplot;           /* GNUPLOT handle                                  */
  double a[ORD+1];       /* Коэффициенты числителя   H(s)                   */
  double b[ORD+1];       /* Коэффициенты знаменателя H(s)                   */
  double Rs = 60.0;      /* Подавление АЧХ в полосе заграждения             */
  double w[N];           /* Циклическая частота рад/c                       */
  double mag[N];         /* АЧХ фильтра                                     */
  double phi[N];         /* ФЧХ фильтра                                     */
  double tau[N];         /* ГВЗ фильтра                                     */
  int k;

  /* Load DSPL function  */
  hdspl = dspl_load();
  if(!hdspl)
  {
    printf("libdspl loading error!\n");
    return -1;
  }
  
  /* рассчитываем нормированный ФНЧ Чебышева 2 рода                         */
  int res = cheby2_ap(Rs, ORD, b, a);
  if(res != RES_OK)
    printf("error code = 0x%8x\n", res);

  /* печать коэффициентов фильтра                                           */
  for(k = 0; k < ORD+1; k++)
    printf("b[%2d] = %9.3f     a[%2d] = %9.3f\n", k, b[k], k, a[k]);

  /* вектор частоты в логарифмической шакале от 0.01 до 100 рад/c           */
  logspace(-2.0, 2.0, N , DSPL_SYMMETRIC, w);

  /* частотные характеристика фильтра                                       */
  filter_freq_resp(b, a, ORD, w, N, 
                   DSPL_FLAG_LOGMAG|DSPL_FLAG_UNWRAP | DSPL_FLAG_ANALOG, 
                   mag, phi, tau);

  /* Сохранить характеристики для построения графиков                       */
  writetxt(w, mag, N, "dat/cheby2_ap_test_mag.txt");
  writetxt(w, phi, N, "dat/cheby2_ap_test_phi.txt");
  writetxt(w, tau, N, "dat/cheby2_ap_test_tau.txt");

  /* plotting by GNUPLOT                                                    */
  gnuplot_create(argc, argv, 920, 260, "img/cheby2_ap_test.png", &hplot);
  gnuplot_cmd(hplot, "set logscale x");
  gnuplot_cmd(hplot, "unset key");
  gnuplot_cmd(hplot, "set grid");
  gnuplot_cmd(hplot, "set xlabel 'frequency, rad/s'");
  gnuplot_cmd(hplot, "set multiplot layout 1,3 rowsfirst");
  gnuplot_cmd(hplot, "set ylabel 'Magnitude, dB'");
  gnuplot_cmd(hplot, "set yrange [-100:5]");
  gnuplot_cmd(hplot, "plot 'dat/cheby2_ap_test_mag.txt' with lines");
  gnuplot_cmd(hplot, "set ylabel 'Phase response, rad'");
  gnuplot_cmd(hplot, "unset yrange");
  gnuplot_cmd(hplot, "plot 'dat/cheby2_ap_test_phi.txt' with lines");
  gnuplot_cmd(hplot, "set ylabel 'Groupdelay, sec'");
  gnuplot_cmd(hplot, "unset yrange");
  gnuplot_cmd(hplot, "plot 'dat/cheby2_ap_test_tau.txt' with lines");
  gnuplot_cmd(hplot, "unset multiplot");
  gnuplot_close(hplot); 

  /* free dspl handle */
  dspl_free(hdspl);

  return res;
}


