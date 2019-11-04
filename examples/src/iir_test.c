#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "dspl.h"

/* Low-pass filter order */
#define LPF_ORD 6

/* High-pass filter order */
#define HPF_ORD 6

/* band-pass filter order */
#define BPF_ORD 12

/* Band-stop filter order */
#define BSF_ORD 12

/* Maximum filter order */
#define MAX_ORD BPF_ORD

/* Stopband suppression (dB) */
#define RS 60.0

/* Pass-band maximum distortion (dB) */
#define RP 2.0

/* Frequency response vector size */
#define N  1024


/******************************************************************************* 
 * function calculates filter frequency response and save magnitude to
 * the text file.
 * params: b - pointer to the transfer fuction H(z) numerator vector
 *         a - pointer to the transfer fuction H(z) denominator vector
 *         ord - filter order
 *         n - number of magnitude vector size
 *         fn - file name
 ******************************************************************************/
void freq_resp_write2txt(double* b, double* a, int ord, int n, char* fn)
{
  double *w = NULL, *mag = NULL;
  int k;
  
  w   = (double*)malloc(n*sizeof(double));
  mag = (double*)malloc(n*sizeof(double));
  
  /* Normalized frequency from 0 to pi */
  linspace(0, M_PI, n , DSPL_PERIODIC, w);
  
  /* Magnitude (dB) calculation */
  filter_freq_resp(b, a, ord, w, n, DSPL_FLAG_LOGMAG, mag, NULL, NULL);
  
  /* Frequency normaliztion from 0 to 1 */
  for(k = 0; k < N; k++)
    w[k] /= M_PI;
  
  /* Save magnitude to the txt file */
  writetxt(w, mag, n, fn);
  
  free(w);
  free(mag);
}


/*******************************************************************************
 * Main program
 ******************************************************************************/
int main(int argc, char* argv[])
{
  void* handle;           /* DSPL handle        */
  handle = dspl_load();   /* Load DSPL function */
  
  /* Transfer function H(z) coeff. vectors */
  double a[MAX_ORD+1], b[MAX_ORD+1];  

  /*--------------------------------------------------------------------------*/
  
  /* LPF Batterworth */
  iir(RP, RS, LPF_ORD, 0.3, 0.0, DSPL_FILTER_BUTTER | DSPL_FILTER_LPF, b, a);
  freq_resp_write2txt(b, a, LPF_ORD, N, "dat/iir_butter_lpf.txt");
  
  /* HPF Batterworth */
  iir(RP, RS, HPF_ORD, 0.3, 0.0, DSPL_FILTER_BUTTER | DSPL_FILTER_HPF, b, a);
  freq_resp_write2txt(b, a, HPF_ORD, N, "dat/iir_butter_hpf.txt");
  
  /* Band-pass Batterworth */
  iir(RP, RS, BPF_ORD, 0.3, 0.7, DSPL_FILTER_BUTTER | DSPL_FILTER_BPASS, b, a);
  freq_resp_write2txt(b, a, BPF_ORD, N, "dat/iir_butter_bpf.txt");
  
  /* Band-stop Batterworth */
  iir(RP, RS, BSF_ORD, 0.3, 0.7, DSPL_FILTER_BUTTER | DSPL_FILTER_BSTOP, b, a);
  freq_resp_write2txt(b, a, BSF_ORD, N, "dat/iir_butter_bsf.txt");
   
  /*--------------------------------------------------------------------------*/
  
  /* LPF Chebyshev type 1 */
  iir(RP, RS, LPF_ORD, 0.3, 0.0, DSPL_FILTER_CHEBY1 | DSPL_FILTER_LPF, b, a);
  freq_resp_write2txt(b, a, LPF_ORD, N, "dat/iir_cheby1_lpf.txt");
  
  /* HPF Chebyshev type 1 */
  iir(RP, RS, HPF_ORD, 0.3, 0.0, DSPL_FILTER_CHEBY1 | DSPL_FILTER_HPF, b, a);
  freq_resp_write2txt(b, a, HPF_ORD, N, "dat/iir_cheby1_hpf.txt");
  
  /* Bnad-pass Chebyshev type 1 */
  iir(RP, RS, BPF_ORD, 0.3, 0.7, DSPL_FILTER_CHEBY1 | DSPL_FILTER_BPASS, b, a);
  freq_resp_write2txt(b, a, BPF_ORD, N, "dat/iir_cheby1_bpf.txt");
  
  /* Bnad-stop Chebyshev type 1 */
  iir(RP, RS, BSF_ORD, 0.3, 0.7, DSPL_FILTER_CHEBY1 | DSPL_FILTER_BSTOP, b, a);
  freq_resp_write2txt(b, a, BSF_ORD, N, "dat/iir_cheby1_bsf.txt");
  
  /*--------------------------------------------------------------------------*/
  
  /* LPF Chebyshev type 2 */
  iir(RP, RS, LPF_ORD, 0.3, 0.0, DSPL_FILTER_CHEBY2 | DSPL_FILTER_LPF, b, a);
  freq_resp_write2txt(b, a, LPF_ORD, N, "dat/iir_cheby2_lpf.txt");
  
  /* HPF Chebyshev type 2 */
  iir(RP, RS, HPF_ORD, 0.3, 0.0, DSPL_FILTER_CHEBY2 | DSPL_FILTER_HPF, b, a);
  freq_resp_write2txt(b, a, HPF_ORD, N, "dat/iir_cheby2_hpf.txt");
  
  /* Band-pass Chebyshev type 2 */
  iir(RP, RS, BPF_ORD, 0.3, 0.7, DSPL_FILTER_CHEBY2 | DSPL_FILTER_BPASS, b, a);
  freq_resp_write2txt(b, a, BPF_ORD, N, "dat/iir_cheby2_bpf.txt");
  
  /* Band-stop Chebyshev type 2 */
  iir(RP, RS, BSF_ORD, 0.3, 0.7, DSPL_FILTER_CHEBY2 | DSPL_FILTER_BSTOP, b, a);
  freq_resp_write2txt(b, a, BSF_ORD, N, "dat/iir_cheby2_bsf.txt");
  
  /*--------------------------------------------------------------------------*/  
  
  /* LPF Elliptic */
  iir(RP, RS, LPF_ORD, 0.3, 0.0, DSPL_FILTER_ELLIP | DSPL_FILTER_LPF, b, a);
  freq_resp_write2txt(b, a, LPF_ORD, N, "dat/iir_ellip_lpf.txt");
  
  /* HPF Elliptic */
  iir(RP, RS, HPF_ORD, 0.3, 0.0, DSPL_FILTER_ELLIP | DSPL_FILTER_HPF, b, a);
  freq_resp_write2txt(b, a, HPF_ORD, N, "dat/iir_ellip_hpf.txt");
  
  /* Band-pass Elliptic */
  iir(RP, RS, BPF_ORD, 0.3, 0.7, DSPL_FILTER_ELLIP | DSPL_FILTER_BPASS, b, a);
  freq_resp_write2txt(b, a, BPF_ORD, N, "dat/iir_ellip_bpf.txt");
  
  /* Band-stop Elliptic */
  iir(RP, RS, BSF_ORD, 0.3, 0.7, DSPL_FILTER_ELLIP | DSPL_FILTER_BSTOP, b, a);
  freq_resp_write2txt(b, a, BSF_ORD, N, "dat/iir_ellip_bsf.txt");
  
  /*--------------------------------------------------------------------------*/
  
  gnuplot_script(argc, argv, "gnuplot/iir_test.plt");
  
  /* free dspl handle */
  dspl_free(handle);
  
  return 0;
}



