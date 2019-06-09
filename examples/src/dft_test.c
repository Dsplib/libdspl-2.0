#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "dspl.h"

#define N 16
int main()
{
  
  void* handle;           // DSPL handle                                   
  handle = dspl_load();   // Load DSPL function

  double    x[N];         // real input signal
  complex_t y[N];         // DFT
  int k;

  for(k = 0; k < N; k++)
    x[k] = (double)k;

  dft(x, N, y);

  for(k = 0; k < N; k++)
    printf("y[%2d] = %9.3f%9.3f\n", k, RE(y[k]), IM(y[k]));

  dspl_free(handle);  // remember to free the resource
  return 0;
}


