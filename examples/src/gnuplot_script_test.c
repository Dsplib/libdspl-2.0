#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "dspl.h"

#define N  100

int main(int argc, char* argv[])
{
  void* handle;           /* DSPL handle        */
  handle = dspl_load();   /* Load DSPL function */
  
  double x[N];
  double s[N];  /* s(x) = sin(x) */
  double c[N];  /* c(x) = cos(x) */
  int n;
  int err;
  
  
  /* x vector from -4*pi to 4*pi */
  linspace(-4.0 * M_PI, 4 * M_PI, N , DSPL_SYMMETRIC, x);
  for(n = 0; n < N; n++)
  {
    s[n] = sin(x[n]); /* s(x) = sin(x) */
    c[n] = cos(x[n]); /* c(x) = cos(x) */
  }    

  /* Save to files "dat/sine.txt" and "dat/cosine.txt" */
  writetxt(x, s, N, "dat/sine.txt");
  writetxt(x, c, N, "dat/cosine.txt");
  
  /* GNUPLOT script gnuplot/gnuplot_script.plt */
  err = gnuplot_script(argc, argv, "gnuplot/gnuplot_script.plt");
  
  /* Print output */
  printf("GNUPLOT err = %d\n", err);
  
  dspl_free(handle);      /* free dspl handle */
  return 0;
}

