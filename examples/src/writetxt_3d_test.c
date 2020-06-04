#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "dspl.h"

#define NX  20
#define NY  30

int main(int argc, char* argv[])
{
  void* hdspl;            /* DSPL handle        */
  void* hplot;         /* GNUPLOT handles    */ 
  
  hdspl = dspl_load();    /* Load DSPL function */
  
  double x[NX];
  double y[NY];
  double z[NX * NY]; 
  int n, m;
  int err;

  /* x vector from -4 to 4 */
  linspace(-2.0, 2.0, NX, DSPL_SYMMETRIC, x);
  
  /* y vector from -3 to 3 */
  linspace(-2.0, 2.0, NY, DSPL_SYMMETRIC, y);
  
  
  for(n = 0; n < NX; n++)
  {
    for(m = 0; m < NY; m++)
    {
      z[n + m*NX] = x[n]*exp(-x[n]*x[n] - y[m]*y[m]);
    }
  }    

  /* Save to files "dat/sine.txt" and "dat/cosine.txt" */
  err = writetxt_3d(x, NX, y, NY, z, "dat/data3d.txt");
  printf("writetxt_3d error 0x%8x\n", err);
    
  /* plotting by GNUPLOT */
  /* Create window 0 */
  err = gnuplot_create(argc, argv, 560, 480, "img/writetxt_3d.png", &hplot);  
  printf("GNUPLOT err = %d\n", err);
  gnuplot_cmd(hplot, "splot 'dat/data3d.txt' with lines");
  gnuplot_close(hplot); 
  
  dspl_free(hdspl);      /* free dspl handle */
  return 0;
}

