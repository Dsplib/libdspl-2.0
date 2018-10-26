#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#include "dspl.h"
#define N 8388608
int main()
{
  void* handle;           // DSPL handle
  handle = dspl_load();   // Load DSPL function
  int n, k, m;

  complex_t *x = NULL;
  complex_t *y = NULL;         // DFT
  fft_t pfft;             // FFT object

  clock_t t0, t1;
  double dt;

  x = (complex_t*)malloc(N*sizeof(complex_t));
  y = (complex_t*)malloc(N*sizeof(complex_t));

  //clear fft object
  memset(&pfft, 0, sizeof(fft_t));

  for(n = 0; n < N; n++)
  {
    RE(x[n]) = (double)n;
    IM(x[n]) = 0.0;
  }

  n = N;
  k = 4;
  while(n > 4)
  {
    printf("n= %16d      ", n);
    fftn_create(&pfft, n);

    t0 = clock();
    for(m = 0; m < k; m++)
      fftn_krn(x, y, &pfft, n, 0);
    t1 = clock();
    dt =  (1000.0 * (double) (t1 - t0)) / CLOCKS_PER_SEC / (double)m;

    printf("%12.6f ms    ", dt);

    fft_free(&pfft);        // clear FFT object


    fft_create(&pfft, n);
    t0 = clock();
    for(m = 0; m < k; m++)
      fft_cmplx(x, n, &pfft, y);
    t1 = clock();
    dt =  (1000.0 * (double) (t1 - t0)) / CLOCKS_PER_SEC / (double)m;

    printf("%12.6f ms \n", dt);

    fft_free(&pfft);        // clear FFT object


    n/=2;
    k = (int)(k*1.5);
  }
  dspl_free(handle);      // free dspl handle
  free(x);
  free(y);
  return 0;
}

