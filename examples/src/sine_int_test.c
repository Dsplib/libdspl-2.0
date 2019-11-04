#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "dspl.h"

#define N 400

int main(int argc, char* argv[])
{
  void* handle;           // DSPL handle
  handle = dspl_load();   // Load DSPL function

  double x[N], y[N];

  linspace(-6*M_PI, 6*M_PI, N, DSPL_PERIODIC, x);

  sine_int(x, N, y);
  writetxt(x, y, N, "dat/dat0.txt");
  
  sinc(x, N, 1.0, y);  
  writetxt(x, y, N, "dat/dat1.txt");
  
  gnuplot_script(argc, argv, "gnuplot/sine_int.plt");
  
  dspl_free(handle);      // free dspl handle
  
  return 0;
}