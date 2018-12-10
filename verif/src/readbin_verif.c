#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "dspl.h"

#define N 10
int main()
{
  void* handle;           // DSPL handle
  handle = dspl_load();   // Load DSPL function
  
  
  double    *xr = NULL, *yr = NULL;
  complex_t *xc = NULL, *yc = NULL;
  double err;
  int vr;
  xr = (double*)    malloc(N * sizeof(double));
  xc = (complex_t*) malloc(N * sizeof(complex_t));
  
  randn(xr,   N, 0, 1.0);
  randn((double*)xc, 2*N, 0, 1.0);
  
  writebin(xr, N, DAT_DOUBLE,  "dat/in_real.dat");
  writebin(xc, N, DAT_COMPLEX, "dat/in_cmplx.dat");
  
  readbin("dat/in_real.dat",  (void**)&yr, NULL, NULL);
  readbin("dat/in_cmplx.dat", (void**)&yc, NULL, NULL);  
  
  vr = verif(xr, yr, N, 1E-12, &err);
  printf("readbin real verification error:          %12.4e\n", err);
  
  vr = verif_cmplx(xc, yc, N, 1E-12, &err);
  printf("readbin cmplx verification error:         %12.4e\n", err);
  
  dspl_free(handle);      // free dspl handle

  if(xr) free(xr);
  if(xc) free(xc);
  if(yr) free(yr);
  if(yc) free(yc);
  return 0;
}


