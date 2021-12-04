#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#include "dspl.h"


#define SIZE   500


int main(int argc, char* argv[])
{
    void* hdspl;            /* DSPL handle        */

    int err, verr, nx, type;
    double derr;
    complex_t  yc;
    double     yd;
    complex_t *xc = NULL;
    double    *xd = NULL;


    hdspl = dspl_load();    /* Load DSPL function */
    
    verif_data_gen(SIZE, DAT_DOUBLE, "dat/real.dat");
    verif_data_gen(SIZE, DAT_COMPLEX, "dat/complex.dat");

    system("octave octave/statistic_verification.m");
    
    
    /*------------------------------------------------------------------------*/
    readbin("dat/real.dat", (void**)(&xd), &nx, &type);
    mean(xd, SIZE, &yd);
    verif_str(&yd, 1, "mean for double data:", 
                      "dat/mean_real.dat", 
                      "verification.log");
    stat_std(xd, SIZE, &yd);
    verif_str(&yd, 1, "stat_std for double data:", "dat/std_real.dat", 
              "verification.log");


    /*------------------------------------------------------------------------*/
    readbin("dat/complex.dat", (void**)(&xc), &nx, &type);
    mean_cmplx(xc, SIZE, &yc);
    verif_str_cmplx(&yc, 1, "mean for complex data:", 
                      "dat/mean_cmplx.dat", 
                      "verification.log");
                      

    /*------------------------------------------------------------------------*/
    stat_std_cmplx(xc, SIZE, &yd);
    verif_str(&yd, 1, "stat_std for complex data:", 
                      "dat/std_cmplx.dat", 
                      "verification.log");


    /* free dspl handle */
    dspl_free(hdspl);
    if(xc)
        free(xc);
    if(xd)
        free(xd);

    return 0;
}

