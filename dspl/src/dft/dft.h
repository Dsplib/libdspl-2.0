#ifndef DFT_H
#define DFT_H


/* sqrt(2^31) */
#define FFT_COMPOSITE_MAX           46340


/* FFT kernel */
int fft_krn(complex_t* t0, complex_t* t1, fft_t* p, int n, int addr);

/* DFT 2 points */
void dft2 (complex_t *x,  complex_t* y);

/* DFT 3 points */
#define DFT3_W         0.866025403784439
void dft3 (complex_t *x,  complex_t* y);

/* DFT 4 points */
void dft4 (complex_t *x,  complex_t* y);


/* DFT 5 points */
#define DFT5_W1       -1.250000000000000
#define DFT5_W2        0.559016994374947
#define DFT5_W3        1.538841768587630
#define DFT5_W4        0.587785252292473
#define DFT5_W5        0.363271264002680
void dft5 (complex_t *x,  complex_t* y);


/* DFT 7 points */
#define DFT7_W1      -1.166666666666666518636930
#define DFT7_W2       0.790156468525400224045541
#define DFT7_W3       0.055854267289647742400494
#define DFT7_W4       0.734302201235752405317214
#define DFT7_W5      -0.440958551844098378680314
#define DFT7_W6      -0.340872930623931369442658
#define DFT7_W7       0.533969360337725240661655
#define DFT7_W8      -0.874842290961656665615465
void dft7 (complex_t *x,  complex_t* y);


/* DFT 8 points */
#define DFT8_W         0.707106781186548
void dft8 (complex_t *x,  complex_t* y);
void transpose2x4(complex_t *x,  complex_t* y);
void transpose4x2(complex_t *x,  complex_t* y);

/* DFT 16 points */
#define DFT16_W1       0.923879532511287
#define DFT16_W2       0.382683432365090
#define DFT16_W3       0.707106781186548
void dft16 (complex_t *x,  complex_t* y);
void transpose4x4(complex_t *x,  complex_t* y);

/* DFT 32 points */
void dft32(complex_t *x, complex_t* y, complex_t* w);
void transpose8x4(complex_t *x, complex_t* y);
void transpose4x8(complex_t *x, complex_t* y);

/* DFT 64 points */
void dft64(complex_t *x, complex_t* y, complex_t* w);
void transpose8x8(complex_t *x, complex_t* y);

/* DFT 128 points */
void dft128(complex_t *x, complex_t* y, complex_t* w);

/* DFT 256 points */
void dft256(complex_t *x, complex_t* y, complex_t* w);
void transpose16x16(complex_t* x, complex_t* y);

/* DFT 512 points */
void dft512(complex_t *x, complex_t* y, complex_t* w, complex_t* w32);

/* DFT 1024 points */
void dft1024(complex_t *x, complex_t* y, complex_t* w, complex_t* w32);


/* DFT 2048 points */
void dft2048(complex_t *x, complex_t* y, complex_t* w, 
             complex_t* w32, complex_t* w64);

/* DFT 4096 points */
void dft4096(complex_t *x, complex_t* y, complex_t* w,  complex_t* w256);


#endif
