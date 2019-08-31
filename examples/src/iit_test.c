#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "dspl.h"

// Порядок фильтра
#define LPF_ORD 6
#define HPF_ORD 6
#define BPF_ORD 12
#define BSF_ORD 12
#define MAX_ORD BPF_ORD
#define RS 60.0
#define RP 1.0

// размер векторов частотной характеристики фильтра
#define N   2048

void freq_resp_write2txt(double* b, double* a, int ord, int n, char* fn)
{
  double w[N], mag[N];
  int k;
  
  // вектор нормированной частоты от 0 до pi
  linspace(0, M_PI, N , DSPL_PERIODIC, w);
  
  //частотные характеристика фильтра
  filter_freq_resp(b, a, ord, w, n, DSPL_FLAG_LOGMAG, mag, NULL, NULL);
  
  for(k = 0; k < N; k++)
    w[k] /= M_PI;
  
  // Сохранить характеристики для построения графиков
  writetxt(w, mag, n, fn);
}


int main()
{
  void* handle;           // DSPL handle
  handle = dspl_load();   // Load DSPL function
  
  double a[MAX_ORD+1], b[MAX_ORD+1];  // коэффициенты H(s)
  int err;
  
  // рассчитываем цифровой ФНЧ  Баттерворта
  iir(RP, RS, LPF_ORD, 0.3, 0.0, DSPL_FILTER_BUTTER | DSPL_FILTER_LPF, b, a);
  freq_resp_write2txt(b, a, LPF_ORD, N, "dat/iir_butter_lpf.txt");
  
  // рассчитываем цифровой ФВЧ  Баттерворта
  iir(RP, RS, HPF_ORD, 0.3, 0.0, DSPL_FILTER_BUTTER | DSPL_FILTER_HPF, b, a);
  freq_resp_write2txt(b, a, HPF_ORD, N, "dat/iir_butter_hpf.txt");
  
  // рассчитываем цифровой полосовой фильтр Баттерворта
  iir(RP, RS, BPF_ORD, 0.3, 0.7, DSPL_FILTER_BUTTER | DSPL_FILTER_BPASS, b, a);
  freq_resp_write2txt(b, a, BPF_ORD, N, "dat/iir_butter_bpf.txt");
  
  // рассчитываем цифровой режекторный фильтр Баттерворта
  iir(RP, RS, BSF_ORD, 0.3, 0.7, DSPL_FILTER_BUTTER | DSPL_FILTER_BSTOP, b, a);
  freq_resp_write2txt(b, a, BSF_ORD, N, "dat/iir_butter_bsf.txt");
  
  
  
  
  // рассчитываем цифровой ФНЧ  Чебышева 1 рода
  iir(RP, RS, LPF_ORD, 0.3, 0.0, DSPL_FILTER_CHEBY1 | DSPL_FILTER_LPF, b, a);
  freq_resp_write2txt(b, a, LPF_ORD, N, "dat/iir_cheby1_lpf.txt");
  
  // рассчитываем цифровой ФВЧ  Чебышева 1 рода
  iir(RP, RS, HPF_ORD, 0.3, 0.0, DSPL_FILTER_CHEBY1 | DSPL_FILTER_HPF, b, a);
  freq_resp_write2txt(b, a, HPF_ORD, N, "dat/iir_cheby1_hpf.txt");
  
  // рассчитываем цифровой полосовой фильтр Чебышева 1 рода
  iir(RP, RS, BPF_ORD, 0.3, 0.7, DSPL_FILTER_CHEBY1 | DSPL_FILTER_BPASS, b, a);
  freq_resp_write2txt(b, a, BPF_ORD, N, "dat/iir_cheby1_bpf.txt");
  
  // рассчитываем цифровой режекторный фильтр Чебышева 1 рода
  iir(RP, RS, BSF_ORD, 0.3, 0.7, DSPL_FILTER_CHEBY1 | DSPL_FILTER_BSTOP, b, a);
  freq_resp_write2txt(b, a, BSF_ORD, N, "dat/iir_cheby1_bsf.txt");
  
  
  
  
  // рассчитываем цифровой ФНЧ  Чебышева 2 рода
  iir(RP, RS, LPF_ORD, 0.3, 0.0, DSPL_FILTER_CHEBY2 | DSPL_FILTER_LPF, b, a);
  freq_resp_write2txt(b, a, LPF_ORD, N, "dat/iir_cheby2_lpf.txt");
  
  // рассчитываем цифровой ФВЧ  Чебышева 2 рода
  iir(RP, RS, HPF_ORD, 0.3, 0.0, DSPL_FILTER_CHEBY2 | DSPL_FILTER_HPF, b, a);
  freq_resp_write2txt(b, a, HPF_ORD, N, "dat/iir_cheby2_hpf.txt");
  
  // рассчитываем цифровой полосовой фильтр Чебышева 2 рода
  iir(RP, RS, BPF_ORD, 0.3, 0.7, DSPL_FILTER_CHEBY2 | DSPL_FILTER_BPASS, b, a);
  freq_resp_write2txt(b, a, BPF_ORD, N, "dat/iir_cheby2_bpf.txt");
  
  // рассчитываем цифровой режекторный фильтр Чебышева 2 рода
  iir(RP, RS, BSF_ORD, 0.3, 0.7, DSPL_FILTER_CHEBY2 | DSPL_FILTER_BSTOP, b, a);
  freq_resp_write2txt(b, a, BSF_ORD, N, "dat/iir_cheby2_bsf.txt");
  
  
  
  
  // рассчитываем цифровой эллиптический ФНЧ
  iir(RP, RS, LPF_ORD, 0.3, 0.0, DSPL_FILTER_ELLIP | DSPL_FILTER_LPF, b, a);
  freq_resp_write2txt(b, a, LPF_ORD, N, "dat/iir_ellip_lpf.txt");
  
  // рассчитываем цифровой эллиптический ФВЧ 
  iir(RP, RS, HPF_ORD, 0.3, 0.0, DSPL_FILTER_ELLIP | DSPL_FILTER_HPF, b, a);
  freq_resp_write2txt(b, a, HPF_ORD, N, "dat/iir_ellip_hpf.txt");
  
  // рассчитываем цифровой полосовой эллиптический фильтр
  iir(RP, RS, BPF_ORD, 0.3, 0.7, DSPL_FILTER_ELLIP | DSPL_FILTER_BPASS, b, a);
  freq_resp_write2txt(b, a, BPF_ORD, N, "dat/iir_ellip_bpf.txt");
  
  // рассчитываем цифровой режекторный эллиптический фильтр
  iir(RP, RS, BSF_ORD, 0.3, 0.7, DSPL_FILTER_ELLIP | DSPL_FILTER_BSTOP, b, a);
  freq_resp_write2txt(b, a, BSF_ORD, N, "dat/iir_ellip_bsf.txt");
  
  
  
  
  dspl_free(handle);      // free dspl handle
  
  // выполнить скрипт GNUPLOT для построения графиков 
  // по рассчитанным данным
  return system("gnuplot -p gnuplot/iir_test.plt");;
}


