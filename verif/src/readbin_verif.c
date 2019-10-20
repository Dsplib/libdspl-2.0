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
  
  random_t rnd;
  random_init(&rnd, RAND_TYPE_MT19937);
  
  randn(xr, N, 0, 1.0, &rnd);
  randn((double*)xc, 2*N, 0, 1.0,&rnd);
  
  
  writebin(xr, N, DAT_DOUBLE,  "dat/in.dat");
  system("octave octave/readbin_verif.m");
  readbin("dat/out.dat",  (void**)&yr, NULL, NULL);
  
  vr = verif(xr, yr, N, 1E-12, &err);
  if(vr == DSPL_VERIF_SUCCESS)
    printf("readbin real verification OK:          %12.4e\n", err);
  else
    printf("readbin real verification ERROR:       %12.4e\n", err);
  
  
  writebin(xc, N, DAT_COMPLEX,  "dat/in.dat");
  system("octave octave/readbin_verif.m");
  readbin("dat/out.dat",  (void**)&yc, NULL, NULL); 
  
  
  vr = verif_cmplx(xc, yc, N, 1E-12, &err);
  if(vr == DSPL_VERIF_SUCCESS)
    printf("readbin cmplx verification OK:         %12.4e\n", err);
  else
    printf("readbin cmplx verification ERROR:      %12.4e\n", err);
  
  dspl_free(handle);      // free dspl handle

  if(xr) free(xr);
  if(xc) free(xc);
  if(yr) free(yr);
  if(yc) free(yc);
  return 0;
}


