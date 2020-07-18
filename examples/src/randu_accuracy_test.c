#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "dspl.h"

#define N  10000
#define L  1.0E-3


int randu_accuracy(double* x, int n, double acc, random_t* prnd)
{
    double z[512];
    int cnt = 0;
    int i;
    int err = RES_OK;
    if(!x)
        return ERROR_PTR;
    if(n < 1)
        return ERROR_SIZE;
    if(acc < 0.0)
        return ERROR_NEGATIVE;

    while(cnt < n)
    {
        randu(z, 512, prnd);
        for(i = 0; i < 512; i++)
        {
            if((z[i] < acc) && (cnt < n))
            {
                x[cnt] = z[i];
                cnt++;
            }
        }
    }
    return RES_OK;
}



int main(int argc, char* argv[])
{
    void* hdspl;            /* DSPL handle        */
    void* hplot;            /* GNUPLOT handles    */
    random_t rnd = {0};     /* random structure   */
    hdspl = dspl_load();    /* Load DSPL function */
    int err;


    double *x = NULL;
    double *y = NULL;

    x = (double*) malloc(N * sizeof(double));
    y = (double*) malloc(N * sizeof(double));

    if(!x || !y)
    {
        err = ERROR_PTR;
        printf("malloc error!\n");
        goto exit_label;
    }


    /* MRG32K3A random numbers generator */
    double seed_mrg32k3a = 1234.0;

    random_init(&rnd, RAND_TYPE_MRG32K3A, (void*)(&seed_mrg32k3a));
    randu_accuracy(x, N, L, &rnd); /* accumulate in x numbers from 0 to L */
    randu_accuracy(y, N, L, &rnd); /* accumulate in y numbers from 0 to L */
    writetxt(x, y, N, "dat/randu_acc_mrg32k3a.txt");


    /* MT19937 random numbers generator  */
    unsigned long long seed_mt19937 = 1234353456;
    random_init(&rnd, RAND_TYPE_MT19937, (void*)(&seed_mt19937));
    randu_accuracy(x, N, L, &rnd); /* accumulate in x numbers from 0 to L */
    randu_accuracy(y, N, L, &rnd); /* accumulate in y numbers from 0 to L */
    writetxt(x, y, N, "dat/randu_acc_mt19937.txt");


    /* Standard C random numbers generator */
    randu_accuracy(x, N, L, NULL); /* accumulate in x numbers from 0 to L */
    randu_accuracy(y, N, L, NULL); /* accumulate in y numbers from 0 to L */
    writetxt(x, y, N, "dat/randu_acc_std.txt");


    /* plotting by GNUPLOT   */
    /* Create window plot    */
    err = gnuplot_create(argc, argv, 920, 320, "img/randu_acc_test.png", &hplot);
    if(err != RES_OK)
        goto exit_label;

    gnuplot_cmd(hplot, "set multiplot layout 1,3 rowsfirst");
    gnuplot_cmd(hplot, "set xlabel 'x'");
    gnuplot_cmd(hplot, "set ylabel 'y'");
    gnuplot_cmd(hplot, "set xtics ('0' 0,'5E-4' 0.0005,' 1E-3' 0.001)");
    gnuplot_cmd(hplot, "set ytics ('0' 0,'5E-4' 0.0005,' 1E-3' 0.001)");
    gnuplot_cmd(hplot, "unset key");
    gnuplot_cmd(hplot, "set title 'MRG32K3A'");
    gnuplot_cmd(hplot, "plot 'dat/randu_acc_mrg32k3a.txt' with points pointtype 0");
    gnuplot_cmd(hplot, "set title 'MT19937'");
    gnuplot_cmd(hplot, "plot 'dat/randu_acc_mt19937.txt'  with points pointtype 0");
    gnuplot_cmd(hplot, "set title 'Standard C'");
    gnuplot_cmd(hplot, "plot 'dat/randu_acc_std.txt'  with points pointtype 0");

exit_label:
    if(x)
        free(x);
    if(y)
        free(y);

    if(hplot)
        gnuplot_close(hplot);

    /* free dspl handle */
    dspl_free(hdspl);
    return 0;
}

