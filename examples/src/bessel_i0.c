#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "dspl.h"

#define N 50

int main(int argc, char* argv[])
{
  void* handle;           // DSPL handle
  handle = dspl_load();   // Load DSPL function

  double x[N], y[N];

  linspace(0.0, 3.0, N, DSPL_SYMMETRIC, x);

  bessel_i0(x, N, y);

  writetxt(x, y, N, "dat/dat0.txt");

  /* run GNUPLOT script */
  
  gnuplot_script(argc, argv, "gnuplot/bessel_i0.plt");

  dspl_free(handle);      // free dspl handle

  return 0;
}