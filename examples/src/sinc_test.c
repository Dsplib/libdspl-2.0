#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "dspl.h"

#define N   1000


int main()
{
  void* handle;           // DSPL handle
  handle = dspl_load();   // Load DSPL function
  
  double x[N], y[N];
  
  linspace(-10.0, 10.0, N , DSPL_SYMMETRIC, x);
  sinc(x, N, 1.0, y);
  writetxt(x, y, N, "dat/sinc_test_1.0.txt");
  
  sinc(x, N, M_PI, y);
  writetxt(x, y, N, "dat/sinc_test_pi.txt");
  
  sinc(x, N, 2.0*M_PI, y);
  writetxt(x, y, N, "dat/sinc_test_2pi.txt");
  
  dspl_free(handle);      // free dspl handle
  
  // выполнить скрипт GNUPLOT для построения графиков
  // по рассчитанным данным
  int err =  system("gnuplot -p  gnuplot/sinc_test.plt");
  printf("err = %d\n", err);
  
  
  return err;
}

