#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "dspl.h"

#define N 250

int main()
{
  void* handle;           // DSPL handle
  handle = dspl_load();   // Load DSPL function

  double x[N], y[N];
  int ord;
  char fn[64];
  
  linspace(-1.0, 1.0, N, DSPL_SYMMETRIC, x);
  for(ord = 1; ord < 5; ord++)
  {
    cheby_poly1(x, N, ord, y);
    sprintf(fn, "dat/cheby_poly1_ord%d.txt", ord);
    writetxt(x,y,N,fn);
  }
  
  dspl_free(handle);      // free dspl handle

  return system("gnuplot  gnuplot/cheby_poly1.plt");
}