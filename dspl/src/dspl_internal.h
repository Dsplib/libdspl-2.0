#ifndef DSPL_INTERNAL_H
#define DSPL_INTERNAL_H


#define DSPL_FARROW_LAGRANGE_COEFF  0.16666666666666666666666666666667

#define DSPL_RAND_MOD_X1            2147483647
#define DSPL_RAND_MOD_X2            2145483479


void transpose(double* a, int n, int m, double* b);
void transpose_cmplx(complex_t* a, int n, int m, complex_t* b);
void transpose_hermite(complex_t* a, int n, int m, complex_t* b);

int fft_bit_reverse(complex_t* x, complex_t* y, int n, int p2);
int fft_dit(fft_t *pfft, int n, complex_t* y);
void fft_dit_krn(complex_t *x0, complex_t *x1, complex_t *w, int n,
                 complex_t *y0, complex_t *y1);

int fft_p2(int n);
void dft2 (complex_t *x,  complex_t* y);

#endif

