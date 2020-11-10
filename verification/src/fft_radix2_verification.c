#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#include "dspl.h"

#define FFT_SIZE   65536


int main(int argc, char* argv[])
{
    void* hdspl;            /* DSPL handle        */
    fft_t pfft = {0};

    int verr, nx, type;
    double derr;
    complex_t *yout = NULL;
    complex_t *xc   = NULL;

    hdspl = dspl_load();    /* Load DSPL function */
    

    verif_data_gen(4,     DAT_COMPLEX, "dat/x_fft_4.dat");
    verif_data_gen(8,     DAT_COMPLEX, "dat/x_fft_8.dat");
    verif_data_gen(16,    DAT_COMPLEX, "dat/x_fft_16.dat");
    verif_data_gen(32,    DAT_COMPLEX, "dat/x_fft_32.dat");
    verif_data_gen(64,    DAT_COMPLEX, "dat/x_fft_64.dat");
    verif_data_gen(128,   DAT_COMPLEX, "dat/x_fft_128.dat");
    verif_data_gen(256,   DAT_COMPLEX, "dat/x_fft_256.dat");
    verif_data_gen(512,   DAT_COMPLEX, "dat/x_fft_512.dat");
    verif_data_gen(1024,  DAT_COMPLEX, "dat/x_fft_1024.dat");
    verif_data_gen(2048,  DAT_COMPLEX, "dat/x_fft_2048.dat");
    verif_data_gen(4096,  DAT_COMPLEX, "dat/x_fft_4096.dat");
    verif_data_gen(8192,  DAT_COMPLEX, "dat/x_fft_8192.dat");
    verif_data_gen(16384, DAT_COMPLEX, "dat/x_fft_16384.dat");
    verif_data_gen(32768, DAT_COMPLEX, "dat/x_fft_32768.dat");
    verif_data_gen(65536, DAT_COMPLEX, "dat/x_fft_65536.dat");
    
    yout = (complex_t*)malloc(FFT_SIZE * sizeof(complex_t  ));

    system("octave octave/fft_radix2_verification.m");
    

    /*------------------------------------------------------------------------*/
    readbin("dat/x_fft_4.dat", (void**)(&xc), &nx, &type);
    fft_cmplx(xc, 4, &pfft, yout);
    verif_str_cmplx(yout, 4, "fft 4    for complex dat", 
                             "dat/y_fft_4.dat", 
                             "verification.log");
    /*------------------------------------------------------------------------*/
    readbin("dat/x_fft_8.dat", (void**)(&xc), &nx, &type);
    fft_cmplx(xc, 8, &pfft, yout);
    verif_str_cmplx(yout, 8, "fft 8    for complex dat", 
                             "dat/y_fft_8.dat", 
                             "verification.log");
    /*------------------------------------------------------------------------*/
    readbin("dat/x_fft_16.dat", (void**)(&xc), &nx, &type);
    fft_cmplx(xc, 16, &pfft, yout);
    verif_str_cmplx(yout, 16, "fft 16   for complex dat", 
                             "dat/y_fft_16.dat", 
                             "verification.log");
                             
    /*------------------------------------------------------------------------*/
    readbin("dat/x_fft_32.dat", (void**)(&xc), &nx, &type);
    fft_cmplx(xc, 32, &pfft, yout);
    verif_str_cmplx(yout, 32, "fft 32   for complex dat", 
                             "dat/y_fft_32.dat", 
                             "verification.log");

    /*------------------------------------------------------------------------*/
    readbin("dat/x_fft_64.dat", (void**)(&xc), &nx, &type);
    fft_cmplx(xc, 64, &pfft, yout);
    verif_str_cmplx(yout, 64, "fft 64   for complex dat", 
                             "dat/y_fft_64.dat", 
                             "verification.log");
    /*------------------------------------------------------------------------*/
    readbin("dat/x_fft_128.dat", (void**)(&xc), &nx, &type);
    fft_cmplx(xc, 128, &pfft, yout);
    verif_str_cmplx(yout, 128, "fft 128  for complex dat", 
                             "dat/y_fft_128.dat", 
                             "verification.log");
    /*------------------------------------------------------------------------*/
    readbin("dat/x_fft_256.dat", (void**)(&xc), &nx, &type);
    fft_cmplx(xc, 256, &pfft, yout);
    verif_str_cmplx(yout, 256, "fft 256   for complex dat", 
                             "dat/y_fft_256.dat", 
                             "verification.log");
                             
    /*------------------------------------------------------------------------*/
    readbin("dat/x_fft_512.dat", (void**)(&xc), &nx, &type);
    fft_cmplx(xc, 512, &pfft, yout);
    verif_str_cmplx(yout, 512, "fft 512   for complex dat", 
                             "dat/y_fft_512.dat", 
                             "verification.log");

    /*------------------------------------------------------------------------*/
    readbin("dat/x_fft_1024.dat", (void**)(&xc), &nx, &type);
    fft_cmplx(xc, 1024, &pfft, yout);
    verif_str_cmplx(yout, 1024, "fft 1024  for complex dat", 
                             "dat/y_fft_1024.dat", 
                             "verification.log");
                             
    /*------------------------------------------------------------------------*/
    readbin("dat/x_fft_2048.dat", (void**)(&xc), &nx, &type);
    fft_cmplx(xc, 2048, &pfft, yout);
    verif_str_cmplx(yout, 2048, "fft 2048  for complex dat", 
                             "dat/y_fft_2048.dat", 
                             "verification.log");
                             
    /*------------------------------------------------------------------------*/
    readbin("dat/x_fft_4096.dat", (void**)(&xc), &nx, &type);
    fft_cmplx(xc, 4096, &pfft, yout);
    verif_str_cmplx(yout, 4096, "fft 4096  for complex dat", 
                             "dat/y_fft_4096.dat", 
                             "verification.log");
  
    /*------------------------------------------------------------------------*/
    readbin("dat/x_fft_8192.dat", (void**)(&xc), &nx, &type);
    fft_cmplx(xc, 8192, &pfft, yout);
    verif_str_cmplx(yout, 8192, "fft 8192  for complex dat", 
                             "dat/y_fft_8192.dat", 
                             "verification.log");

    /*------------------------------------------------------------------------*/
    readbin("dat/x_fft_16384.dat", (void**)(&xc), &nx, &type);
    fft_cmplx(xc, 16384, &pfft, yout);
    verif_str_cmplx(yout, 16384, "fft 16384  for complex dat", 
                             "dat/y_fft_16384.dat", 
                             "verification.log");

    /*------------------------------------------------------------------------*/
    readbin("dat/x_fft_32768.dat", (void**)(&xc), &nx, &type);
    fft_cmplx(xc, 32768, &pfft, yout);
    verif_str_cmplx(yout, 32768, "fft 32768  for complex dat", 
                             "dat/y_fft_32768.dat", 
                             "verification.log");
                             
    /*------------------------------------------------------------------------*/
    readbin("dat/x_fft_65536.dat", (void**)(&xc), &nx, &type);
    fft_cmplx(xc, 65536, &pfft, yout);
    verif_str_cmplx(yout, 65536, "fft 65536  for complex dat", 
                             "dat/y_fft_65536.dat", 
                             "verification.log");


    /* free dspl handle */
    dspl_free(hdspl);
    
    if(yout)
        free(yout);
    if(xc)
        free(xc);

    return 0;
}

