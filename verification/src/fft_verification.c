#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#include "dspl.h"


/* 2^15 */
#define FFT_R2_SIZE   32768
/* 2 * 3 * 5 * 7 * 11 * 13 */  
#define FFT_CT_SIZE   30030

int main(int argc, char* argv[])
{
    void* hdspl;            /* DSPL handle        */
    fft_t pfft = {0};

    int verr, nx, type;
    double derr;
    complex_t *yout = NULL;
    complex_t *xc   = NULL;
    double    *xd   = NULL;
    
    hdspl = dspl_load();    /* Load DSPL function */
    
    verif_data_gen(FFT_R2_SIZE, DAT_DOUBLE,  "dat/x_fft_double_radix2.dat");
    verif_data_gen(FFT_R2_SIZE, DAT_COMPLEX, "dat/x_fft_complex_radix2.dat");
    verif_data_gen(FFT_CT_SIZE, DAT_DOUBLE,  "dat/x_fft_double_common.dat");
    verif_data_gen(FFT_CT_SIZE, DAT_COMPLEX, "dat/x_fft_complex_common.dat");

    yout = (complex_t*)malloc(FFT_R2_SIZE * sizeof(complex_t  ));

    system("octave octave/fft_verification.m");
    
    
    /*------------------------------------------------------------------------*/
    readbin("dat/x_fft_double_radix2.dat", (void**)(&xd), &nx, &type);
    fft(xd, FFT_R2_SIZE, &pfft, yout);
    verif_str_cmplx(yout, FFT_R2_SIZE, "fft radix-2 for double dat", 
                                        "dat/y_fft_double_radix2.dat", 
                                        "verification.log");

    /*------------------------------------------------------------------------*/
    readbin("dat/x_fft_complex_radix2.dat", (void**)(&xc), &nx, &type);
    fft_cmplx(xc, FFT_R2_SIZE, &pfft, yout);
    verif_str_cmplx(yout, FFT_R2_SIZE, "fft radix-2 for complex dat", 
                                         "dat/y_fft_complex_radix2.dat", 
                                         "verification.log");

    /*------------------------------------------------------------------------*/
    readbin("dat/x_fft_double_common.dat", (void**)(&xd), &nx, &type);
    fft(xd, FFT_CT_SIZE, &pfft, yout);
    verif_str_cmplx(yout, FFT_CT_SIZE, "fft composite for double dat", 
                                        "dat/y_fft_double_common.dat", 
                                        "verification.log");

    /*------------------------------------------------------------------------*/
    readbin("dat/x_fft_complex_common.dat", (void**)(&xc), &nx, &type);
    fft_cmplx(xc, FFT_CT_SIZE, &pfft, yout);
    verif_str_cmplx(yout, FFT_CT_SIZE, "fft composite for complex dat", 
                                         "dat/y_fft_complex_common.dat", 
                                         "verification.log");
  
    /* free dspl handle */
    dspl_free(hdspl);
    
    if(yout)
        free(yout);
    if(xc)
        free(xc);
    if(xd)
        free(xd);

    return 0;
}

