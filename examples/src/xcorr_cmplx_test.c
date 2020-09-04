#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "dspl.h"

#define N 200
#define P 20

int main(int argc, char* argv[])
{
    void* hdspl;            /* DSPL handle        */
    void* hplot;            /* GNUPLOT handles    */
    
    random_t rnd = {0};     /* random structure   */
    fft_t    fft = {0};
    
    hdspl = dspl_load();    /* Load DSPL function */

    complex_t x[N];
    complex_t z[2*P+1];
    double t[2*P+1];
    
    int err;
    
    err = random_init(&rnd, RAND_TYPE_MRG32K3A, NULL);
    printf("err = %.8x\n", err);
    if(err != RES_OK)
        goto exit_label;

    err = randn((double*)x, 2*N, 0.0, 1.0, &rnd);
    printf("err = %.8x\n", err);
    if(err != RES_OK)
        goto exit_label;

    err = xcorr_cmplx(x, N, x, N, DSPL_XCORR_BIASED, P, z, t);
    printf("err = %.8x\n", err);
    writetxt_cmplx_re(t, z, 2*P+1, "dat/xcorr_b_re.txt");
    writetxt_cmplx_im(t, z, 2*P+1, "dat/xcorr_b_im.txt");

    err = xcorr_cmplx(x, N, x, N, DSPL_XCORR_UNBIASED, P, z, t);
    printf("err = %.8x\n", err);
    writetxt_cmplx_re(t, z, 2*P+1, "dat/xcorr_u_re.txt");
    writetxt_cmplx_im(t, z, 2*P+1, "dat/xcorr_u_im.txt");
 
    /**************************************************************************/
    /* plotting by GNUPLOT                                                    */
    /**************************************************************************/
    /* Create window plot */
    err = gnuplot_create(argc, argv, 420, 420, "img/xcorr_re.png", &hplot);
    if(err != RES_OK)
        goto exit_label;

    gnuplot_cmd(hplot, "set xlabel 't'");
    gnuplot_cmd(hplot, "set ylabel 'z'");
    gnuplot_cmd(hplot, "unset key");
    gnuplot_cmd(hplot, "set title 'xcorr'");
    gnuplot_cmd(hplot, "plot 'dat/xcorr_u_re.txt' with lines,\\");
    gnuplot_cmd(hplot, "     'dat/xcorr_b_re.txt' with lines");
    gnuplot_close(hplot);

    err = gnuplot_create(argc, argv, 420, 420, "img/xcorr_im.png", &hplot);
    if(err != RES_OK)
        goto exit_label;

    gnuplot_cmd(hplot, "set xlabel 't'");
    gnuplot_cmd(hplot, "set ylabel 'z'");
    gnuplot_cmd(hplot, "unset key");
    gnuplot_cmd(hplot, "set title 'xcorr'");
    gnuplot_cmd(hplot, "plot 'dat/xcorr_u_im.txt' with lines,\\");
    gnuplot_cmd(hplot, "     'dat/xcorr_b_im.txt' with lines");
    gnuplot_close(hplot);

exit_label:
    /* free dspl handle */
    dspl_free(hdspl);
    return err;
}

