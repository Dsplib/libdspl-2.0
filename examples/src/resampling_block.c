#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "dspl.h"

/* Input signal length */
#define N   256

/* Block size */
#define L   16

/* P/Q */
#define P   160
#define Q   147


int main(int argc, char* argv[])
{
    void* hdspl;  /* DSPL handle        */
    void* hplot;  /* GNUPLOT handle     */

    double s[N], ts[N];   /* input signal and time array               */
    double dt;            /* Discretization time on the output Q/P/Fs  */
    double *ty = NULL;    /* Time array for output signal              */
    double *y = NULL;     /* Block output resampler array              */
    double *z = NULL;     /* Full signal resampling                    */

    int k, ny;
    int err;

    /* Load DSPL function  */
    hdspl = dspl_load();

    /* Fill input signal array */
    for(k = 0; k < N; k++)
    {
        ts[k] = (double)k;
        s[k] = sin(M_2PI*0.2176870748*ts[k]);
    }
    /* Discretization time on the output Q/P/Fs  */
    dt = (double)Q/(double)P;
    
    /* Resampler output samples number */
    ny = (int)((double)(N-1)/dt)+1;
    
    y  = (double*) malloc(ny * sizeof(double));
    ty = (double*) malloc(ny * sizeof(double));
    
     /* Fill resampler output time array */
    for(k = 0; k < ny; k++)
        ty[k] = (double)k*dt;
    
    /* Vraiblaes for time recalc. frd - current block fractional delay */
    double ts0,ts1, ty0, ty1, frd;
    /* input and output array position */
    int pos, ypos; 

    ts0 = ty0 = ts1 = ty1 = 0.0;
    frd = 0.0;
    pos = ypos =  0;
    
    /* current block output tmp array */
    double* tmp = NULL;
    int ntmp;
    
    while(pos+L < N)
    {
        /* Current block resampling from pos (L - Block size) */
        farrow_spline(s+pos, L, P, Q, frd, &tmp, &ntmp);
        
        /*  Copy to the output array. 2 samples overlapping */ 
        if(ypos)
            memcpy(y+ypos+2, tmp+2, (ntmp-2)*sizeof(double));
        else
            memcpy(y, tmp, (ntmp)*sizeof(double));

        /* time recalc to the end of block */
        ty1 = ty0+(double)ntmp*dt;
        ts1 = ts0+(double)L;
        printf("pos = %d, ts1 = %.3f, ty1 = %.3f\n", pos,ts1,ty1);
        
        /* Output array position for the next */
        ypos += ntmp;
        /* 4 samples input array overlapping */
        ts0  = ts1 - 4.0;
        
        /* Resampler output position search */
        ty0 = ty1;
        while(ty0>ts0)
        {
            ty0-=dt;
            ypos--;
        }
        
        /* frac delay recalc */
        frd = ts0 - ty0; 
        printf("ts0 = %.3f, ty0 = %.3f frd = %.3f\n", ts0,ty0, frd);
        
        /* 4 samples input array overlapping */
        pos += L-4;
    }
    /* residual block resampling */
    farrow_spline(s+pos, N-pos, P, Q, frd, &tmp, &ntmp);
    memcpy(y+ypos+2, tmp+2, (ntmp-2)*sizeof(double));

    /* full signal resampling */
    farrow_spline(s, N, P, Q, 0, &z, &k);

    /* 
    save input signal and filter output to the txt-files.
    y.txt and z.txt are identical 
    */
    writetxt(ts, s, N,  "dat/s.txt");
    writetxt(ty, y, ny, "dat/y.txt");
    writetxt(ty, z, ny, "dat/z.txt");

    /* plotting by GNUPLOT*/
    gnuplot_create(argc, argv, 820, 540, "img/resampling_block.png", &hplot);
    gnuplot_cmd(hplot, "unset key");
    gnuplot_cmd(hplot, "set grid");
    gnuplot_cmd(hplot, "set xlabel 'n'");
    gnuplot_cmd(hplot, "set ylabel 's(n)'");
    gnuplot_cmd(hplot, "set yrange [-1.5:1.5]");
    gnuplot_cmd(hplot, "set multiplot layout 2,1 rowsfirst");
    gnuplot_cmd(hplot, "plot 'dat/s.txt'  with lines");

    gnuplot_cmd(hplot, "plot 'dat/y.txt' with lines,\\");
    gnuplot_cmd(hplot, "     'dat/z.txt' with lines");
    gnuplot_cmd(hplot, "unset multiplot");
    gnuplot_close(hplot);
    
    /* free DSPL handle */
    dspl_free(hdspl);

    if(ty)
        free(ty);
    if(y)
        free(y);
    if(z)
        free(z);
    if(tmp)
        free(tmp);
    return err;
}
