#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "dspl.h"

int main()
{
  void* handle;           // DSPL handle
  handle = dspl_load();   // Load DSPL function
  
  double a[5] = {0.0, 1.0, 2.0, 3.0, 4.0};
  double b[5] = {5.0, 6.0, 7.0, 8.0, 9.0};
  double c[10], d[5];
  int err, k, n;
  
  // Concatenate arrays a and b. Result keeps to the array c
  err = concat((void*)a, 5*sizeof(double), 
               (void*)b, 5*sizeof(double), (void*)c);
  printf("\n\nconcatenation result: %d\n\narray c = ", err);
  for(k = 0; k < 10; k++)
    printf("%6.1f", c[k]);
  
  // Decimate array c 2 times. Result keeps to the array d
  err = decimate(c, 10, 2, d, &n);
  printf("\n\ndecimation result:    %d\n\narray d = ", err);
  for (k = 0; k < n; k++)
    printf("%6.1f", d[k]);
  
  dspl_free(handle);      // free dspl handle
  
  return err;
}


