#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#include "dspl.h"

int main(int argc, char* argv[])
{
    void* hdspl;            /* DSPL handle        */
    complex_t    *pxd = NULL;
    complex_t    *pyd = NULL;
    char str[512] = {0};
    
    int nx, ny, tx, err, verr;
    double derr;
    
    random_t rnd = {0};     /* random structure   */
    
    hdspl = dspl_load();    /* Load DSPL function */
    
    sprintf(str, "writebin and readbin for complex data:");
    while(strlen(str) < 48)
      str[strlen(str)] = 46;
    
    
    err = random_init(&rnd, RAND_TYPE_MRG32K3A, NULL);
    if(err != RES_OK)
        goto exit_error_code;

    err = randi(&nx, 1, 1, 1000000, &rnd);
    if(err != RES_OK)
        goto exit_error_code;
    
    pxd = (complex_t*) malloc(nx * sizeof(complex_t));
    
    
    /****************** verifiсation function start  **************************/
    err  = randn((double*)pxd, 2*nx, 1.0, 10.0, &rnd);
    if(err != RES_OK)
        goto exit_error_code;
      
    /************ Write input verification data to the file  ******************/
    err = writebin((void*)pxd, nx, DAT_COMPLEX, "dat/x.dat");
    if(err != RES_OK)
        goto exit_error_code;
    
    /************ RUN external verificator (octave or python)  ****************/
    err = system("octave octave/writebin_readbin_verification.m");

    /***************** Read external verificator output  **********************/
    err = readbin("dat/y.dat", (void**)(&pyd), &ny, &tx);
    if(err != RES_OK)
        goto exit_error_code;
    if(tx!=DAT_COMPLEX)
    {
        err = ERROR_DAT_TYPE;
        goto exit_error_code;
    }
    
    /**************************** Verification  *******************************/
    verr = verif_cmplx(pxd, pyd, ny, 1E-12, &derr);
    if(verr != DSPL_VERIF_SUCCESS)
        goto exit_error_verif;

     sprintf(str, "%s ok (err = %12.4E)", str,  derr);
     goto exit_label;

exit_error_code:
     sprintf(str, "%s FAILED (with code = 0x%8x)", str, err);
     goto exit_label;

exit_error_verif:
     sprintf(str, "%s FAILED (err = %12.4E)", str, derr);

exit_label:
    /************************ write str to log file  **************************/
    
    
    if(pxd)
        free(pxd);
    if(pyd)
        free(pyd);
      
    /* free dspl handle */
    dspl_free(hdspl);
    
    printf("%s\n", str);
    return 0;
}

