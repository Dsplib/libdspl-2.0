#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "dspl.h"

void roots_sort(complex_t* r, int n)
{
    int i, k;
    i = 0;
    complex_t tmp, diff; 
    double min;
    while(i < n)
    {
        RE(diff) = RE(r[i]) - RE(r[i+1]);
        IM(diff) = IM(r[i]) + IM(r[i+1]);
        min = ABS(diff);
        for(k = i+2; k < n; k++)
        { 
            RE(diff) = RE(r[i]) - RE(r[k]);
            IM(diff) = IM(r[i]) + IM(r[k]);
            if(ABS(diff) < min)
            {
                min = ABS(diff);
                RE(tmp) = RE(r[i+1]);
                IM(tmp) = IM(r[i+1]);
                RE(r[i+1]) = RE(r[k]);
                IM(r[i+1]) = IM(r[k]);
                RE(r[k]) = RE(tmp);
                IM(r[k]) = IM(tmp);
            }
        }
        i+=2;
    }
}


int ab2sos(double* b, double* a, int ord, double** sosb, double** sosa)
{
    double *af = NULL;
    double *bf = NULL;
    complex_t *rb = NULL;
    complex_t *ra = NULL;
    double *bg = NULL;
    double gain;
    int err, info, k, n;
    
    complex_t tmpa[2];
    complex_t tmpb[2];
    complex_t tmpc[3];
    
    af = (double*)malloc((ord+1)*sizeof(double));
    bf = (double*)malloc((ord+1)*sizeof(double));
    
    rb = (complex_t*)malloc((ord)*sizeof(complex_t));
    ra = (complex_t*)malloc((ord)*sizeof(complex_t));
    
    memcpy(af, a, (ord+1)*sizeof(double));
    memcpy(bf, b, (ord+1)*sizeof(double));
    
    err = flipip(af, ord+1);
    if(err!=RES_OK)
        goto error_proc;
    
    err = flipip(bf, ord+1);
    if(err!=RES_OK)
        goto error_proc;
    
    
    
    // roots calculation
    err = polyroots(bf, ord, rb, &info);
    if(err!=RES_OK)
        goto error_proc;
    err = polyroots(af, ord, ra, &info);
    if(err!=RES_OK)
        goto error_proc;

    // roots conjugate sorting
    roots_sort(rb, ord);
    roots_sort(ra, ord);
     
    // SOS generation
    for (n = 0; n < ord/2; n++)
    { 
        RE(tmpa[0]) = RE(tmpb[0]) = 1.0;
        IM(tmpa[0]) = IM(tmpb[0]) = 0.0;
        RE(tmpa[1]) = -RE(ra[n*2]);
        RE(tmpb[1]) = -RE(ra[n*2+1]);
        IM(tmpa[1]) = -IM(ra[n*2]);
        IM(tmpb[1]) = -IM(ra[n*2+1]);
        
        err = conv_cmplx(tmpa, 2, tmpb, 2, tmpc);
        if(err!=RES_OK)
            goto error_proc;
    
        for(k = 0; k<3; k++)
            sosa[n][k] = RE(tmpc[k]);
        
        RE(tmpa[0]) = RE(tmpb[0]) = 1.0;
        IM(tmpa[0]) = IM(tmpb[0]) = 0.0;
        RE(tmpa[1]) = -RE(rb[n*2]);
        RE(tmpb[1]) = -RE(rb[n*2+1]);
        IM(tmpa[1]) = -IM(rb[n*2]);
        IM(tmpb[1]) = -IM(rb[n*2+1]);
        
        
        err = conv_cmplx(tmpa, 2, tmpb, 2, tmpc);
        if(err!=RES_OK)
            goto error_proc;
        
        for(k = 0; k<3; k++)
            sosb[n][k] = RE(tmpc[k]);
    }
    
    // SOS gain balancing
    bg = (double*)malloc((ord+1)*sizeof(double));
    memset(bg, 0, (ord+1)*sizeof(double));
    bg[0] = 1.0;
    k = 1;
    for(n = 0; n < ord/2; n++)
    {
        conv(bg, k, sosb[n], 3, bg);
        k+=2;
    }
    
    gain = pow(b[0]/bg[0], 2.0/ord);
    for(n = 0; n < ord/2; n++)
        for(k = 0; k < 3; k++)
            sosb[n][k] *= gain; 
        
error_proc:        
    if(af)
        free(af);
    if(bf)
        free(bf);
    if(rb)
        free(rb);
    if(ra)
        free(ra);
    if(bg)
        free(bg);
    return err;
    
}





#define N      1000

#define ORD   6
/* Stopband suppression (dB) */
#define RS 60.0
 
/* Pass-band maximum distortion (dB) */
#define RP 2.0


int main(int argc, char* argv[])
{
    void* hdspl;    /* DSPL handle */
    void* hplot;    /* GNUPLOT handle */

    double w[N], h[N], g[N];
    
    double b[ORD+1];
    double a[ORD+1];
    
    double B[ORD+1];
    double A[ORD+1];

    double **sosa = NULL;
    double **sosb = NULL;
        
    int err, k, n;

    hdspl = dspl_load();     
    
    
    // Генерирую фильтр порядка ORD (должно быть четным числом)
    iir(RP, RS, ORD, 0.1, 0.3, DSPL_FILTER_BUTTER | DSPL_FILTER_BSTOP, b, a);
    
    // Фильтр порядка ORD раскладывается в ORD/2 секций SOS, в каждой
    // по 3 к-та a и b
    // выделяю память под двумерный массив [ORD/2 x 3] коэффициентов a и b
    sosa = (double**)malloc((ORD/2)*sizeof(double*));
    sosb = (double**)malloc((ORD/2)*sizeof(double*));
    for (n = 0; n < ORD/2; n++)
    {
        sosa[n] = (double*)malloc(3*sizeof(double));
        sosb[n] = (double*)malloc(3*sizeof(double));
    }

    // Эта функция возьмет коэффициенты b и a и разложит их в SOS массивы
    err = ab2sos(b, a, ORD, sosb, sosa);
    
    // Все SOS выведем на печать
    for (n = 0; n < ORD/2; n++) 
    {
        printf("SOS %d\n ------------------------ \n", n);
        for(k = 0; k < 3; k++)
        {
            printf("b[%d] = %12.6f  a[%d] = %12.6f\n", k, sosb[n][k], k, sosa[n][k]);
        }
    }
    
    
    // А вот тут проверка, что перемножив  H(z) всех SOS 
    // получим к-ты исходного фильтра
    memset(B,0,(ORD+1)*sizeof(double));
    memset(A,0,(ORD+1)*sizeof(double));
    B[0] = 1.0;
    A[0] = 1.0;
    k = 1;
    for(n = 0; n < ORD/2; n++)
    {
        conv(B, k, sosb[n], 3, B);
        conv(A, k, sosa[n], 3, A);
        k+=2;
    }
    
    printf(" ------------------------ \n");
    printf("Test coefficients\n");
    for (n = 0; n < ORD+1; n++)
        printf("B[%d] = %12.6f  b[%d] = %12.6f   A[%d] = %12.6f  a[%d] = %12.6f\n", n, B[n], n, b[n],n, A[n], n, a[n]);
    
    
    // память чистим
    for (n = 0; n < ORD/2; n++)
    {
        if(sosa[n])
            free(sosa[n]);
        if(sosb[n])
            free(sosb[n]);
    }
    if(sosa)
        free(sosa);
    if(sosb)
        free(sosb);
    
    /* free dspl handle */
    dspl_free(hdspl);

    return err;
}

