#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "dspl.h"

// Filter order (must be even for bandstop IIR)
#define ORD 14

// Frequency response vector size
#define N   1001


int main()
{
  void* handle;               // DSPL handle
  handle = dspl_load();       // Load DSPL function
  
  double a[ORD+1], b[ORD+1];  // H(s) coefficients
  double rs = 60.0;  // Bandstop suppression equals 60 dB
  double rp = 1.0;   // Bandpass ripple equals 1 dB
  
  // Frequency (w), magnitude (mag), phase response (phi) 
  // and group delay (tau)
  double w[N], mag[N], phi[N], tau[N];
  int k;

  // Calculate Chebyshev type 2 digital bandstop filter 
  int res = iir(rp, rs, ORD, 0.3, 0.7, 
                DSPL_FILTER_CHEBY2 | DSPL_FILTER_BSTOP, b, a);
  if(res != RES_OK)
    printf("error code = 0x%8x\n", res);
  
  // Print coefficients
  for(k = 0; k < ORD+1; k++)
    printf("b[%2d] = %9.3f     a[%2d] = %9.3f\n", k, b[k], k, a[k]);
  
  // Normalized circular frequency from 0 to pi
  linspace(0, M_PI, N , DSPL_PERIODIC, w);
  
  // Filter magnitude, phase response and group delay 
  filter_freq_resp(b, a, ORD, w, N, DSPL_FLAG_LOGMAG|DSPL_FLAG_UNWRAP,  
                   mag, phi, tau);
  
  // Normalized frequency from 0 to 1.
  // w = 1 corresponds to Fs/2
  for(k = 0; k < N; k++)
    w[k] /= M_PI;
  
  // Save filter frequency response to the txt-files
  // for plotting by GNUPLOT
  writetxt(w, mag, N, "dat/iir_bstop_mag.txt");
  writetxt(w, phi, N, "dat/iir_bstop_phi.txt");
  writetxt(w, tau, N, "dat/iir_bstop_tau.txt");
  
  dspl_free(handle);      // free dspl handle
  
  // run GNUPLOT script
  return system("gnuplot -p gnuplot/iir_bstop.plt");;
}


