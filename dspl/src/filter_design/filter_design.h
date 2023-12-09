#ifndef FILTER_DESIGN_H
#define FILTER_DESIGN_H


int fir_linphase_lpf(int ord, double wp, int wintype, 
                     double winparam, double* h);
                     

int iir_ap(double rp, double rs, int ord, int type, double* b, double* a);

#endif

