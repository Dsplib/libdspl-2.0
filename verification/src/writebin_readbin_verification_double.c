#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#include "dspl.h"

int main(int argc, char* argv[])
{
    void* hdspl;            /* DSPL handle        */
    double    *pxd = NULL;
    double    *pyd = NULL;
    
    int nx, ny, tx, err, verr;
    double derr;
    
    random_t rnd = {0};     /* random structure   */
    
    hdspl = dspl_load();    /* Load DSPL function */
    
    printf("writebin and readbin for double data:..........");
    
    err = random_init(&rnd, RAND_TYPE_MRG32K3A, NULL);
    if(err != RES_OK)
        goto exit_error_code;

    err = randi(&nx, 1, 1, 1000000, &rnd);
    if(err != RES_OK)
        goto exit_error_code;
    
    pxd = (double*) malloc(nx * sizeof(double));
    
    err  = randn(pxd, nx, 1.0, 10.0, &rnd);
    if(err != RES_OK)
        goto exit_error_code;
    
    err = writebin((void*)pxd, nx, DAT_DOUBLE, "dat/x.dat");
    if(err != RES_OK)
        goto exit_error_code;

    err = system("octave octave/writebin_readbin_verification.m");

    err = readbin("dat/y.dat", (void**)(&pyd), &ny, &tx);
    if(err != RES_OK)
        goto exit_error_code;
      
    if(tx!=DAT_DOUBLE)
    {
        err = ERROR_DAT_TYPE;
        goto exit_error_code;
    }
    verr = verif(pxd, pyd, ny, 1E-12, &derr);
    if(verr != DSPL_VERIF_SUCCESS)
        goto exit_error_verif;

     printf("ok (err = %12.4E)", derr);
     goto exit_label;
     
exit_error_code:
     printf("FAILED (with code = 0x%8x)", err);
     goto exit_label;
exit_error_verif:
     printf("FAILED (err = %12.4E)", derr);
    
exit_label:
    if(pxd)
        free(pxd);
    if(pyd)
        free(pyd);
    /* free dspl handle */
    dspl_free(hdspl);
    return 0;
}

