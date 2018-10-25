#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "dspl.h"

#define N 16
int main()
{
  void* handle;           // DSPL handle
  handle = dspl_load();   // Load DSPL function
  int n;

  complex_t x[N];
  complex_t y[N];         // DFT
  fft_t pfft;             // FFT object

  //clear fft object
  memset(&pfft, 0, sizeof(fft_t));

  // Create FFT object
  fftn_create(&pfft, N);

  for(n = 0; n < N; n++)
  {
    RE(x[n]) = (double)n;
    IM(x[n]) = 0.0;
  }

  fftn_krn(x, y, &pfft, N, 0);


  fft_free(&pfft);        // clear FFT object

  dspl_free(handle);      // free dspl handle
  return 0;
}

