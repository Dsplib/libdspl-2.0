#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "dspl.h"

#define N 250

int main(int argc, char* argv[])
{
  void* handle;           // DSPL handle
  handle = dspl_load();   // Load DSPL function

  double x[N], y[N];
  int ord;
  char fn[64];
  
  linspace(-1.0, 1.0, N, DSPL_SYMMETRIC, x);
  for(ord = 1; ord < 5; ord++)
  {
    cheby_poly2(x, N, ord, y);
    sprintf(fn, "dat/cheby_poly2_ord%d.txt", ord);
    writetxt(x,y,N,fn);
  }
  
  /* run GNUPLOT script */
  gnuplot_script(argc, argv, "gnuplot/cheby_poly2.plt");
  
  dspl_free(handle);      // free dspl handle

  return 0;
}