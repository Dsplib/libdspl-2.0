#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "dspl.h"

#define ORD 6
#define N   2000




int main(int argc, char* argv[])
{
  void* handle;   /* DSPL handle         */

  double b[ORD+1], a[ORD+1];
  double t[N], s[N], n[N], sf[N];
  random_t rnd;
  int k;
  int err;

  /* Load DSPL function  */
  handle = dspl_load();

  /* random generator init */
  random_init(&rnd, RAND_TYPE_MT19937, NULL);

  /* fill time vector         */
  linspace(0, N, N, DSPL_PERIODIC, t);

  /* generate noise        */
  randn(n, N, 0, 1.0, &rnd);

  /* input signal s = sin(2*pi*t) + n(t) */
  for(k = 0; k < N; k++)
    s[k] = sin(M_2PI*0.02*t[k]) + n[k];

  /* IIR filter coefficients calculation */
  iir(1.0, 70.0, ORD, 0.06, 0.0, DSPL_FILTER_ELLIP | DSPL_FILTER_LPF, b, a);

  /* input signal filtration */
  filter_iir(b, a, ORD, s, N, sf);

  /* save input signal and filter output to the txt-files */
  writetxt(t,s, N, "dat/s.txt");
  writetxt(t,sf,N, "dat/sf.txt");

  /* run GNUPLOT script */
  err = gnuplot_script(argc, argv, "gnuplot/filter_iir.plt");

  /* free DSPL handle */
  dspl_free(handle);


  return err;
}
