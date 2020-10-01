#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#include "dspl.h"

#define SIZE 10000

int main(int argc, char* argv[])
{
    void* hdspl;            /* DSPL handle        */
    complex_t *xc = NULL;
    double    *xd = NULL;

    int nx, type, err, verr;
    double derr;
   
    
    hdspl = dspl_load();    /* Load DSPL function */
    
    /* generation real random data for verification */ 
    verif_data_gen(SIZE, DAT_DOUBLE, "dat/real.dat");
    
    /* generation complex random data for verification */ 
    verif_data_gen(SIZE, DAT_COMPLEX, "dat/complex.dat");

    /* RUN verification in octave */
    system("octave octave/writebin_readbin_verification.m");


    /*------------------------------------------------------------------------*/
    /* Read real input data from the file */ 
    readbin("dat/real.dat", (void**)(&xd), &nx, &type);
    
    /* Processing */
    
    /* verification libdspl output and octave output */
    verif_str(xd, SIZE, "writebin and readbin for real data:", 
                        "dat/yreal.dat", 
                        "verification.log");


    /*------------------------------------------------------------------------*/
    /* Read complex input data from the file */ 
    readbin("dat/complex.dat", (void**)(&xc), &nx, &type);
    
    /* Processing */
    
    /* verification libdspl output and octave output */
    verif_str_cmplx(xc, SIZE, "writebin and readbin for complex data:", 
                               "dat/ycomplex.dat", 
                               "verification.log");


    /* free dspl handle */
    dspl_free(hdspl);


    if(xc)
        free(xc);
    if(xd)
        free(xd);
    
    return 0;
}

