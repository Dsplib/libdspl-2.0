#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "dspl.h"

#define N    1000
#define ORD  4

int main()
{
  void* handle;           // DSPL handle
  handle = dspl_load();   // Load DSPL function
  
  double w[N], h[N];
  complex_t hz[N];
  double bs[ORD+1], as[ORD+1];
  double bz[ORD+1], az[ORD+1];
  int err, k;
  
  // расчет аналогового фильтра Чебышева первого рода
  err = cheby1_ap(1.0, ORD, bs, as);
  if(err != RES_OK)
  {
    printf("cheby1_ap error code %d\n", err);
    return err;
  }
  
  // Билинейное преобразование
  err = bilinear(bs, as, ORD, bz, az);
  if(err != RES_OK)
  {
    printf("bilinear error code %d\n", err);
    return err;
  }
  
  // Печать коэффициентов
  for(k = 0; k < ORD+1; k++)
    printf("bz[%d] = %7.3f    az[%d] = %7.3f\n", k, bz[k], k, az[k]);
  
  
  // Расчет АЧХ полученного цифрового фильтра
  linspace(0, M_PI, N, DSPL_PERIODIC, w);
  freqz(bz, az, ORD, w, N, hz);
  for(k = 0; k < N; k++)
  {
    h[k] = 10.0 * log10 (ABSSQR(hz[k])); // АЧХ в дБ
    w[k] /= M_PI;                        // нормировка частоты от 0 до 1
  }
  
  writetxt(w,h,N,"dat/bilinear.txt");
  
  dspl_free(handle);      // free dspl handle
  
  // выполнить скрипт GNUPLOT для построения графиков
  // по рассчитанным данным
  err =  system("gnuplot -p  gnuplot/bilinear_test.plt");
  printf("err = %d\n", err);
  
  return err;
}

