#include <stdio.h>
#include <stdlib.h>
#ifndef BLAS_H
#define BLAS_H

#define FORTRAN_FUNC(FUNC) FUNC##_

int    FORTRAN_FUNC(xerbla)(const char*, int*info, int);

float  FORTRAN_FUNC(sdot)  (int*, float*, int*, float*, int*);
float  FORTRAN_FUNC(sdsdot)(int*, float*,        float*, int*, float*, int*);

double FORTRAN_FUNC(dsdot) (int*, float*, int*, float*, int*);
double FORTRAN_FUNC(ddot)  (int*, double*, int*, double*, int*);
double FORTRAN_FUNC(qdot)  (int*, double*, int*, double*, int*);

int    FORTRAN_FUNC(cdotuw)  (int*, float*, int*, float*, int*, float*);
int    FORTRAN_FUNC(cdotcw)  (int*, float*, int*, float*, int*, float*);
int    FORTRAN_FUNC(zdotuw)  (int*, double*, int*, double*, int*, double*);
int    FORTRAN_FUNC(zdotcw)  (int*, double*, int*, double*, int*, double*);

int    FORTRAN_FUNC(saxpy) (const int*, const float*, const float*, const int*, float*, const int*);
int    FORTRAN_FUNC(daxpy) (const int*, const double*, const double*, const int*, double*, const int*);
int    FORTRAN_FUNC(qaxpy) (const int*, const double*, const double*, const int*, double*, const int*);
int    FORTRAN_FUNC(caxpy) (const int*, const float*, const float*, const int*, float*, const int*);
int    FORTRAN_FUNC(zaxpy) (const int*, const double*, const double*, const int*, double*, const int*);
int    FORTRAN_FUNC(xaxpy) (const int*, const double*, const double*, const int*, double*, const int*);
int    FORTRAN_FUNC(caxpyc)(const int*, const float*, const float*, const int*, float*, const int*);
int    FORTRAN_FUNC(zaxpyc)(const int*, const double*, const double*, const int*, double*, const int*);
int    FORTRAN_FUNC(xaxpyc)(const int*, const double*, const double*, const int*, double*, const int*);

int    FORTRAN_FUNC(scopy) (int*, float*, int*, float*, int*);
int    FORTRAN_FUNC(dcopy) (int*, double*, int*, double*, int*);
int    FORTRAN_FUNC(qcopy) (int*, double*, int*, double*, int*);
int    FORTRAN_FUNC(ccopy) (int*, float*, int*, float*, int*);
int    FORTRAN_FUNC(zcopy) (int*, double*, int*, double*, int*);
int    FORTRAN_FUNC(xcopy) (int*, double*, int*, double*, int*);

int    FORTRAN_FUNC(sswap) (int*, float*, int*, float*, int*);
int    FORTRAN_FUNC(dswap) (int*, double*, int*, double*, int*);
int    FORTRAN_FUNC(qswap) (int*, double*, int*, double*, int*);
int    FORTRAN_FUNC(cswap) (int*, float*, int*, float*, int*);
int    FORTRAN_FUNC(zswap) (int*, double*, int*, double*, int*);
int    FORTRAN_FUNC(xswap) (int*, double*, int*, double*, int*);

float  FORTRAN_FUNC(sasum) (int*, float*, int*);
float  FORTRAN_FUNC(scasum)(int*, float*, int*);
double FORTRAN_FUNC(dasum) (int*, double*, int*);
double FORTRAN_FUNC(qasum) (int*, double*, int*);
double FORTRAN_FUNC(dzasum)(int*, double*, int*);
double FORTRAN_FUNC(qxasum)(int*, double*, int*);

int    FORTRAN_FUNC(isamax)(int*, float*, int*);
int    FORTRAN_FUNC(idamax)(int*, double*, int*);
int    FORTRAN_FUNC(iqamax)(int*, double*, int*);
int    FORTRAN_FUNC(icamax)(int*, float*, int*);
int    FORTRAN_FUNC(izamax)(int*, double*, int*);
int    FORTRAN_FUNC(ixamax)(int*, double*, int*);

int    FORTRAN_FUNC(ismax) (int*, float*, int*);
int    FORTRAN_FUNC(idmax) (int*, double*, int*);
int    FORTRAN_FUNC(iqmax) (int*, double*, int*);
int    FORTRAN_FUNC(icmax) (int*, float*, int*);
int    FORTRAN_FUNC(izmax) (int*, double*, int*);
int    FORTRAN_FUNC(ixmax) (int*, double*, int*);

int    FORTRAN_FUNC(isamin)(int*, float*, int*);
int    FORTRAN_FUNC(idamin)(int*, double*, int*);
int    FORTRAN_FUNC(iqamin)(int*, double*, int*);
int    FORTRAN_FUNC(icamin)(int*, float*, int*);
int    FORTRAN_FUNC(izamin)(int*, double*, int*);
int    FORTRAN_FUNC(ixamin)(int*, double*, int*);

int    FORTRAN_FUNC(ismin)(int*, float*, int*);
int    FORTRAN_FUNC(idmin)(int*, double*, int*);
int    FORTRAN_FUNC(iqmin)(int*, double*, int*);
int    FORTRAN_FUNC(icmin)(int*, float*, int*);
int    FORTRAN_FUNC(izmin)(int*, double*, int*);
int    FORTRAN_FUNC(ixmin)(int*, double*, int*);

float  FORTRAN_FUNC(samax) (int*, float*, int*);
double FORTRAN_FUNC(damax) (int*, double*, int*);
double FORTRAN_FUNC(qamax) (int*, double*, int*);
float  FORTRAN_FUNC(scamax)(int*, float*, int*);
double FORTRAN_FUNC(dzamax)(int*, double*, int*);
double FORTRAN_FUNC(qxamax)(int*, double*, int*);

float  FORTRAN_FUNC(samin) (int*, float*, int*);
double FORTRAN_FUNC(damin) (int*, double*, int*);
double FORTRAN_FUNC(qamin) (int*, double*, int*);
float  FORTRAN_FUNC(scamin)(int*, float*, int*);
double FORTRAN_FUNC(dzamin)(int*, double*, int*);
double FORTRAN_FUNC(qxamin)(int*, double*, int*);

float  FORTRAN_FUNC(smax)  (int*, float*, int*);
double FORTRAN_FUNC(dmax)  (int*, double*, int*);
double FORTRAN_FUNC(qmax)  (int*, double*, int*);
float  FORTRAN_FUNC(scmax) (int*, float*, int*);
double FORTRAN_FUNC(dzmax) (int*, double*, int*);
double FORTRAN_FUNC(qxmax) (int*, double*, int*);

float  FORTRAN_FUNC(smin)  (int*, float*, int*);
double FORTRAN_FUNC(dmin)  (int*, double*, int*);
double FORTRAN_FUNC(qmin)  (int*, double*, int*);
float  FORTRAN_FUNC(scmin) (int*, float*, int*);
double FORTRAN_FUNC(dzmin) (int*, double*, int*);
double FORTRAN_FUNC(qxmin) (int*, double*, int*);

int    FORTRAN_FUNC(sscal) (int*,  float*, float*, int*);
int    FORTRAN_FUNC(dscal) (int*,  double*, double*, int*);
int    FORTRAN_FUNC(qscal) (int*,  double*, double*, int*);
int    FORTRAN_FUNC(cscal) (int*,  float*, float*, int*);
int    FORTRAN_FUNC(zscal) (int*,  double*, double*, int*);
int    FORTRAN_FUNC(xscal) (int*,  double*, double*, int*);
int    FORTRAN_FUNC(csscal)(int*,  float*, float*, int*);
int    FORTRAN_FUNC(zdscal)(int*,  double*, double*, int*);
int    FORTRAN_FUNC(xqscal)(int*,  double*, double*, int*);

float  FORTRAN_FUNC(snrm2) (int*, float*, int*);
float  FORTRAN_FUNC(scnrm2)(int*, float*, int*);

double FORTRAN_FUNC(dnrm2) (int*, double*, int*);
double FORTRAN_FUNC(qnrm2) (int*, double*, int*);
double FORTRAN_FUNC(dznrm2)(int*, double*, int*);
double FORTRAN_FUNC(qxnrm2)(int*, double*, int*);

int    FORTRAN_FUNC(srot)  (int*, float*, int*, float*, int*, float*, float*);
int    FORTRAN_FUNC(drot)  (int*, double*, int*, double*, int*, double*, double*);
int    FORTRAN_FUNC(qrot)  (int*, double*, int*, double*, int*, double*, double*);
int    FORTRAN_FUNC(csrot) (int*, float*, int*, float*, int*, float*, float*);
int    FORTRAN_FUNC(zdrot) (int*, double*, int*, double*, int*, double*, double*);
int    FORTRAN_FUNC(xqrot) (int*, double*, int*, double*, int*, double*, double*);

int    FORTRAN_FUNC(srotg) (float*, float*, float*, float*);
int    FORTRAN_FUNC(drotg) (double*, double*, double*, double*);
int    FORTRAN_FUNC(qrotg) (double*, double*, double*, double*);
int    FORTRAN_FUNC(crotg) (float*, float*, float*, float*);
int    FORTRAN_FUNC(zrotg) (double*, double*, double*, double*);
int    FORTRAN_FUNC(xrotg) (double*, double*, double*, double*);

int    FORTRAN_FUNC(srotmg)(float*, float*, float*, float*, float*);
int    FORTRAN_FUNC(drotmg)(double*, double*, double*, double*, double*);

int    FORTRAN_FUNC(srotm) (int*, float*, int*, float*, int*, float*);
int    FORTRAN_FUNC(drotm) (int*, double*, int*, double*, int*, double*);
int    FORTRAN_FUNC(qrotm) (int*, double*, int*, double*, int*, double*);

/* Level 2 routines*/

int FORTRAN_FUNC(sger)(int*,    int*, float*,  float*, int*,  float*,  int*, float*,  int*);
int FORTRAN_FUNC(dger)(int*,    int*, double*, double*, int*, double*, int*, double*, int*);
int FORTRAN_FUNC(qger)(int*,    int*, double*, double*, int*, double*, int*, double*, int*);
int FORTRAN_FUNC(cgeru)(int*,    int*, float*,  float*, int*, float*,  int*, float*,  int*);
int FORTRAN_FUNC(cgerc)(int*,    int*, float*,  float*, int*, float*,  int*, float*,  int*);
int FORTRAN_FUNC(zgeru)(int*,    int*, double*, double*, int*, double*, int*, double*, int*);
int FORTRAN_FUNC(zgerc)(int*,    int*, double*, double*, int*, double*, int*, double*, int*);
int FORTRAN_FUNC(xgeru)(int*,    int*, double*, double*, int*, double*, int*, double*, int*);
int FORTRAN_FUNC(xgerc)(int*,    int*, double*, double*, int*, double*, int*, double*, int*);

int FORTRAN_FUNC(sgemv)(const char*, const int*, const int*, const float*, const float*, const int*, const float*, const int*, const float*, float*, const int*);
int FORTRAN_FUNC(dgemv)(const char*, const int*, const int*, const double*, const double*, const int*, const double*, const int*, const double*, double*, const int*);
int FORTRAN_FUNC(qgemv)(const char*, const int*, const int*, const double*, const double*, const int*, const double*, const int*, const double*, double*, const int*);
int FORTRAN_FUNC(cgemv)(const char*, const int*, const int*, const float*, const float*, const int*, const float*, const int*, const float*, float*, const int*);
int FORTRAN_FUNC(zgemv)(const char*, const int*, const int*, const double*, const double*, const int*, const double*, const int*, const double*, double*, const int*);
int FORTRAN_FUNC(xgemv)(const char*, const int*, const int*, const double*, const double*, const int*, const double*, const int*, const double*, double*, const int*);

int FORTRAN_FUNC(strsv) (const char*, const char*, const char*, const int*, const float*, const int*, float*, const int*);
int FORTRAN_FUNC(dtrsv) (const char*, const char*, const char*, const int*, const double*, const int*, double*, const int*);
int FORTRAN_FUNC(qtrsv) (const char*, const char*, const char*, const int*, const double*, const int*, double*, const int*);
int FORTRAN_FUNC(ctrsv) (const char*, const char*, const char*, const int*, const float*, const int*, float*, const int*);
int FORTRAN_FUNC(ztrsv) (const char*, const char*, const char*, const int*, const double*, const int*, double*, const int*);
int FORTRAN_FUNC(xtrsv) (const char*, const char*, const char*, const int*, const double*, const int*, double*, const int*);

int FORTRAN_FUNC(stpsv) (char*, char*, char*, int*, float*, float*, int*);
int FORTRAN_FUNC(dtpsv) (char*, char*, char*, int*, double*, double*, int*);
int FORTRAN_FUNC(qtpsv) (char*, char*, char*, int*, double*, double*, int*);
int FORTRAN_FUNC(ctpsv) (char*, char*, char*, int*, float*, float*, int*);
int FORTRAN_FUNC(ztpsv) (char*, char*, char*, int*, double*, double*, int*);
int FORTRAN_FUNC(xtpsv) (char*, char*, char*, int*, double*, double*, int*);

int FORTRAN_FUNC(strmv) (const char*, const char*, const char*, const int*, const float*, const int*, float*, const int*);
int FORTRAN_FUNC(dtrmv) (const char*, const char*, const char*, const int*, const double*, const int*, double*, const int*);
int FORTRAN_FUNC(qtrmv) (const char*, const char*, const char*, const int*, const double*, const int*, double*, const int*);
int FORTRAN_FUNC(ctrmv) (const char*, const char*, const char*, const int*, const float*, const int*, float*, const int*);
int FORTRAN_FUNC(ztrmv) (const char*, const char*, const char*, const int*, const double*, const int*, double*, const int*);
int FORTRAN_FUNC(xtrmv) (const char*, const char*, const char*, const int*, const double*, const int*, double*, const int*);

int FORTRAN_FUNC(stpmv) (char*, char*, char*, int*, float*, float*, int*);
int FORTRAN_FUNC(dtpmv) (char*, char*, char*, int*, double*, double*, int*);
int FORTRAN_FUNC(qtpmv) (char*, char*, char*, int*, double*, double*, int*);
int FORTRAN_FUNC(ctpmv) (char*, char*, char*, int*, float*, float*, int*);
int FORTRAN_FUNC(ztpmv) (char*, char*, char*, int*, double*, double*, int*);
int FORTRAN_FUNC(xtpmv) (char*, char*, char*, int*, double*, double*, int*);

int FORTRAN_FUNC(stbmv) (char*, char*, char*, int*, int*, float*, int*, float*, int*);
int FORTRAN_FUNC(dtbmv) (char*, char*, char*, int*, int*, double*, int*, double*, int*);
int FORTRAN_FUNC(qtbmv) (char*, char*, char*, int*, int*, double*, int*, double*, int*);
int FORTRAN_FUNC(ctbmv) (char*, char*, char*, int*, int*, float*, int*, float*, int*);
int FORTRAN_FUNC(ztbmv) (char*, char*, char*, int*, int*, double*, int*, double*, int*);
int FORTRAN_FUNC(xtbmv) (char*, char*, char*, int*, int*, double*, int*, double*, int*);

int FORTRAN_FUNC(stbsv) (char*, char*, char*, int*, int*, float*, int*, float*, int*);
int FORTRAN_FUNC(dtbsv) (char*, char*, char*, int*, int*, double*, int*, double*, int*);
int FORTRAN_FUNC(qtbsv) (char*, char*, char*, int*, int*, double*, int*, double*, int*);
int FORTRAN_FUNC(ctbsv) (char*, char*, char*, int*, int*, float*, int*, float*, int*);
int FORTRAN_FUNC(ztbsv) (char*, char*, char*, int*, int*, double*, int*, double*, int*);
int FORTRAN_FUNC(xtbsv) (char*, char*, char*, int*, int*, double*, int*, double*, int*);

int FORTRAN_FUNC(ssymv) (const char*, const int*, const float*, const float*, const int*, const float*, const int*, const float*, float*, const int*);
int FORTRAN_FUNC(dsymv) (const char*, const int*, const double*, const double*, const int*, const double*, const int*, const double*, double*, const int*);
int FORTRAN_FUNC(qsymv) (const char*, const int*, const double*, const double*, const int*, const double*, const int*, const double*, double*, const int*);

int FORTRAN_FUNC(sspmv) (char*, int*, float*, float*, float*, int*, float*, float*, int*);
int FORTRAN_FUNC(dspmv) (char*, int*, double*, double*, double*, int*, double*, double*, int*);
int FORTRAN_FUNC(qspmv) (char*, int*, double*, double*, double*, int*, double*, double*, int*);

int FORTRAN_FUNC(ssyr) (const char*, const int*, const float *, const float*, const int*, float*, const int*);
int FORTRAN_FUNC(dsyr) (const char*, const int*, const double*, const double*, const int*, double*, const int*);
int FORTRAN_FUNC(qsyr) (const char*, const int*, const double*, const double*, const int*, double*, const int*);

int FORTRAN_FUNC(ssyr2) (const char*, const int*, const float *, const float*, const int*, const float*, const int*, float*, const int*);
int FORTRAN_FUNC(dsyr2) (const char*, const int*, const double*, const double*, const int*, const double*, const int*, double*, const int*);
int FORTRAN_FUNC(qsyr2) (const char*, const int*, const double*, const double*, const int*, const double*, const int*, double*, const int*);
int FORTRAN_FUNC(csyr2) (const char*, const int*, const float *, const float*, const int*, const float*, const int*, float*, const int*);
int FORTRAN_FUNC(zsyr2) (const char*, const int*, const double*, const double*, const int*, const double*, const int*, double*, const int*);
int FORTRAN_FUNC(xsyr2) (const char*, const int*, const double*, const double*, const int*, const double*, const int*, double*, const int*);

int FORTRAN_FUNC(sspr) (char*, int*, float *, float*, int*, float*);
int FORTRAN_FUNC(dspr) (char*, int*, double*, double*, int*, double*);
int FORTRAN_FUNC(qspr) (char*, int*, double*, double*, int*, double*);

int FORTRAN_FUNC(sspr2) (char*, int*, float *, float*, int*, float*, int*, float*);
int FORTRAN_FUNC(dspr2) (char*, int*, double*, double*, int*, double*, int*, double*);
int FORTRAN_FUNC(qspr2) (char*, int*, double*, double*, int*, double*, int*, double*);
int FORTRAN_FUNC(cspr2) (char*, int*, float *, float*, int*, float*, int*, float*);
int FORTRAN_FUNC(zspr2) (char*, int*, double*, double*, int*, double*, int*, double*);
int FORTRAN_FUNC(xspr2) (char*, int*, double*, double*, int*, double*, int*, double*);

int FORTRAN_FUNC(cher) (char*, int*, float *, float*, int*, float*, int*);
int FORTRAN_FUNC(zher) (char*, int*, double*, double*, int*, double*, int*);
int FORTRAN_FUNC(xher) (char*, int*, double*, double*, int*, double*, int*);

int FORTRAN_FUNC(chpr) (char*, int*, float *, float*, int*, float*);
int FORTRAN_FUNC(zhpr) (char*, int*, double*, double*, int*, double*);
int FORTRAN_FUNC(xhpr) (char*, int*, double*, double*, int*, double*);

int FORTRAN_FUNC(cher2) (char*, int*, float *, float*, int*, float*, int*, float*, int*);
int FORTRAN_FUNC(zher2) (char*, int*, double*, double*, int*, double*, int*, double*, int*);
int FORTRAN_FUNC(xher2) (char*, int*, double*, double*, int*, double*, int*, double*, int*);

int FORTRAN_FUNC(chpr2) (char*, int*, float *, float*, int*, float*, int*, float*);
int FORTRAN_FUNC(zhpr2) (char*, int*, double*, double*, int*, double*, int*, double*);
int FORTRAN_FUNC(xhpr2) (char*, int*, double*, double*, int*, double*, int*, double*);

int FORTRAN_FUNC(chemv) (const char*, const int*, const float*, const float*, const int*, const float*, const int*, const float*, float*, const int*);
int FORTRAN_FUNC(zhemv) (const char*, const int*, const double*, const double*, const int*, const double*, const int*, const double*, double*, const int*);
int FORTRAN_FUNC(xhemv) (const char*, const int*, const double*, const double*, const int*, const double*, const int*, const double*, double*, const int*);

int FORTRAN_FUNC(chpmv) (char*, int*, float*, float*, float*, int*, float*, float*, int*);
int FORTRAN_FUNC(zhpmv) (char*, int*, double*, double*, double*, int*, double*, double*, int*);
int FORTRAN_FUNC(xhpmv) (char*, int*, double*, double*, double*, int*, double*, double*, int*);

int FORTRAN_FUNC(snorm)(char*, int*, int*, float*, int*);
int FORTRAN_FUNC(dnorm)(char*, int*, int*, double*, int*);
int FORTRAN_FUNC(cnorm)(char*, int*, int*, float*, int*);
int FORTRAN_FUNC(znorm)(char*, int*, int*, double*, int*);

int FORTRAN_FUNC(sgbmv)(char*, int*, int*, int*, int*, float*, float*, int*, float*, int*, float*, float*, int*);
int FORTRAN_FUNC(dgbmv)(char*, int*, int*, int*, int*, double*, double*, int*, double*, int*, double*, double*, int*);
int FORTRAN_FUNC(qgbmv)(char*, int*, int*, int*, int*, double*, double*, int*, double*, int*, double*, double*, int*);
int FORTRAN_FUNC(cgbmv)(char*, int*, int*, int*, int*, float*, float*, int*, float*, int*, float*, float*, int*);
int FORTRAN_FUNC(zgbmv)(char*, int*, int*, int*, int*, double*, double*, int*, double*, int*, double*, double*, int*);
int FORTRAN_FUNC(xgbmv)(char*, int*, int*, int*, int*, double*, double*, int*, double*, int*, double*, double*, int*);

int FORTRAN_FUNC(ssbmv)(char*, int*, int*, float*, float*, int*, float*, int*, float*, float*, int*);    
int FORTRAN_FUNC(dsbmv)(char*, int*, int*, double*, double*, int*, double*, int*, double*, double*, int*);    
int FORTRAN_FUNC(qsbmv)(char*, int*, int*, double*, double*, int*, double*, int*, double*, double*, int*);    
int FORTRAN_FUNC(csbmv)(char*, int*, int*, float*, float*, int*, float*, int*, float*, float*, int*);    
int FORTRAN_FUNC(zsbmv)(char*, int*, int*, double*, double*, int*, double*, int*, double*, double*, int*);    
int FORTRAN_FUNC(xsbmv)(char*, int*, int*, double*, double*, int*, double*, int*, double*, double*, int*);
                                                                     
int FORTRAN_FUNC(chbmv)(char*, int*, int*, float*, float*, int*, float*, int*, float*, float*, int*);
int FORTRAN_FUNC(zhbmv)(char*, int*, int*, double*, double*, int*, double*, int*, double*, double*, int*);
int FORTRAN_FUNC(xhbmv)(char*, int*, int*, double*, double*, int*, double*, int*, double*, double*, int*);

/* Level 3 routines*/

int FORTRAN_FUNC(sgemm)(const char*, const char*, const int*, const int*, const int*, const float*, const float*, const int*, const float*, const int*, const float*, float*, const int*);
int FORTRAN_FUNC(dgemm)(const char*, const char*, const int*, const int*, const int*, const double*, const double*, const int*, const double*, const int*, const double*, double*, const int*);
int FORTRAN_FUNC(qgemm)(const char*, const char*, const int*, const int*, const int*, const double*, const double*, const int*, const double*, const int*, const double*, double*, const int*);
int FORTRAN_FUNC(cgemm)(const char*, const char*, const int*, const int*, const int*, const float*, const float*, const int*, const float*, const int*, const float*, float*, const int*);
int FORTRAN_FUNC(zgemm)(const char*, const char*, const int*, const int*, const int*, const double*, const double*, const int*, const double*, const int*, const double*, double*, const int*);
int FORTRAN_FUNC(xgemm)(const char*, const char*, const int*, const int*, const int*, const double*, const double*, const int*, const double*, const int*, const double*, double*, const int*);

int FORTRAN_FUNC(cgemm3m)(char*, char*, int*, int*, int*, float*,  float*, int*, float*, int*, float*, float*, int*);
int FORTRAN_FUNC(zgemm3m)(char*, char*, int*, int*, int*, double*, double*, int*, double*, int*, double*, double*, int*);
int FORTRAN_FUNC(xgemm3m)(char*, char*, int*, int*, int*, double*, double*, int*, double*, int*, double*, double*, int*);

int FORTRAN_FUNC(sge2mm)(char*, char*, char*, int*, int*, float*, float*, int*, float*, int*, float*, float*, int*);
int FORTRAN_FUNC(dge2mm)(char*, char*, char*, int*, int*, double*, double*, int*, double*, int*, double*, double*, int*);
int FORTRAN_FUNC(cge2mm)(char*, char*, char*, int*, int*, float*, float*, int*, float*, int*, float*, float*, int*);
int FORTRAN_FUNC(zge2mm)(char*, char*, char*, int*, int*, double*, double*, int*, double*, int*, double*, double*, int*);

int FORTRAN_FUNC(strsm)(const char*, const char*, const char*, const char*, const int*, const int*, const float*,  const float*,  const int*, float*,  const int*);
int FORTRAN_FUNC(dtrsm)(const char*, const char*, const char*, const char*, const int*, const int*, const double*, const double*, const int*, double*, const int*);
int FORTRAN_FUNC(qtrsm)(const char*, const char*, const char*, const char*, const int*, const int*, const double*, const double*, const int*, double*, const int*);
int FORTRAN_FUNC(ctrsm)(const char*, const char*, const char*, const char*, const int*, const int*, const float*,  const float*,  const int*, float*,  const int*);
int FORTRAN_FUNC(ztrsm)(const char*, const char*, const char*, const char*, const int*, const int*, const double*, const double*, const int*, double*, const int*);
int FORTRAN_FUNC(xtrsm)(const char*, const char*, const char*, const char*, const int*, const int*, const double*, const double*, const int*, double*, const int*);

int FORTRAN_FUNC(strmm)(const char*, const char*, const char*, const char*, const int*, const int*, const float*,  const float*,  const int*, float*,  const int*);
int FORTRAN_FUNC(dtrmm)(const char*, const char*, const char*, const char*, const int*, const int*, const double*, const double*, const int*, double*, const int*);
int FORTRAN_FUNC(qtrmm)(const char*, const char*, const char*, const char*, const int*, const int*, const double*, const double*, const int*, double*, const int*);
int FORTRAN_FUNC(ctrmm)(const char*, const char*, const char*, const char*, const int*, const int*, const float*,  const float*,  const int*, float*,  const int*);
int FORTRAN_FUNC(ztrmm)(const char*, const char*, const char*, const char*, const int*, const int*, const double*, const double*, const int*, double*, const int*);
int FORTRAN_FUNC(xtrmm)(const char*, const char*, const char*, const char*, const int*, const int*, const double*, const double*, const int*, double*, const int*);

int FORTRAN_FUNC(ssymm)(const char*, const char*, const int*, const int*, const float*, const float*, const int*, const float*, const int*, const float*, float*, const int*);
int FORTRAN_FUNC(dsymm)(const char*, const char*, const int*, const int*, const double*, const double*, const int*, const double*, const int*, const double*, double*, const int*);
int FORTRAN_FUNC(qsymm)(const char*, const char*, const int*, const int*, const double*, const double*, const int*, const double*, const int*, const double*, double*, const int*);
int FORTRAN_FUNC(csymm)(const char*, const char*, const int*, const int*, const float*, const float*, const int*, const float*, const int*, const float*, float*, const int*);
int FORTRAN_FUNC(zsymm)(const char*, const char*, const int*, const int*, const double*, const double*, const int*, const double*, const int*, const double*, double*, const int*);
int FORTRAN_FUNC(xsymm)(const char*, const char*, const int*, const int*, const double*, const double*, const int*, const double*, const int*, const double*, double*, const int*);

int FORTRAN_FUNC(csymm3m)(char*, char*, int*, int*, float*, float*, int*, float*, int*, float*, float*, int*);
int FORTRAN_FUNC(zsymm3m)(char*, char*, int*, int*, double*, double*, int*, double*, int*, double*, double*, int*);
int FORTRAN_FUNC(xsymm3m)(char*, char*, int*, int*, double*, double*, int*, double*, int*, double*, double*, int*);

int FORTRAN_FUNC(ssyrk)(const char*, const char*, const int*, const int*, const float*, const float*, const int*, const float*, float*, const int*);
int FORTRAN_FUNC(dsyrk)(const char*, const char*, const int*, const int*, const double*, const double*, const int*, const double*, double*, const int*);
int FORTRAN_FUNC(qsyrk)(const char*, const char*, const int*, const int*, const double*, const double*, const int*, const double*, double*, const int*);
int FORTRAN_FUNC(csyrk)(const char*, const char*, const int*, const int*, const float*, const float*, const int*, const float*, float*, const int*);
int FORTRAN_FUNC(zsyrk)(const char*, const char*, const int*, const int*, const double*, const double*, const int*, const double*, double*, const int*);
int FORTRAN_FUNC(xsyrk)(const char*, const char*, const int*, const int*, const double*, const double*, const int*, const double*, double*, const int*);

int FORTRAN_FUNC(ssyr2k)(const char*, const char*, const int*, const int*, const float*, const float*, const int*, const float*, const int*, const float*, float*, const int*);
int FORTRAN_FUNC(dsyr2k)(const char*, const char*, const int*, const int*, const double*, const double*, const int*, const double*, const int*, const double*, double*, const int*);
int FORTRAN_FUNC(qsyr2k)(const char*, const char*, const int*, const int*, const double*, const double*, const int*, const double*, const int*, const double*, double*, const int*);
int FORTRAN_FUNC(csyr2k)(const char*, const char*, const int*, const int*, const float*, const float*, const int*, const float*, const int*, const float*, float*, const int*);
int FORTRAN_FUNC(zsyr2k)(const char*, const char*, const int*, const int*, const double*, const double*, const int*, const double*, const int*, const double*, double*, const int*);
int FORTRAN_FUNC(xsyr2k)(const char*, const char*, const int*, const int*, const double*, const double*, const int*, const double*, const int*, const double*, double*, const int*);

int FORTRAN_FUNC(chemm)(const char*, const char*, const int*, const int*, const float*, const float*, const int*, const float*, const int*, const float*, float*, const int*);
int FORTRAN_FUNC(zhemm)(const char*, const char*, const int*, const int*, const double*, const double*, const int*, const double*, const int*, const double*, double*, const int*);
int FORTRAN_FUNC(xhemm)(const char*, const char*, const int*, const int*, const double*, const double*, const int*, const double*, const int*, const double*, double*, const int*);

int FORTRAN_FUNC(chemm3m)(char*, char*, int*, int*, float*, float*, int*, float*, int*, float*, float*, int*);
int FORTRAN_FUNC(zhemm3m)(char*, char*, int*, int*, double*, double*, int*, double*, int*, double*, double*, int*);
int FORTRAN_FUNC(xhemm3m)(char*, char*, int*, int*, double*, double*, int*, double*, int*, double*, double*, int*);

int FORTRAN_FUNC(cherk)(const char*, const char*, const int*, const int*, const float*, const float*, const int*, const float*, float*, const int*);
int FORTRAN_FUNC(zherk)(const char*, const char*, const int*, const int*, const double*, const double*, const int*, const double*, double*, const int*);
int FORTRAN_FUNC(xherk)(const char*, const char*, const int*, const int*, const double*, const double*, const int*, const double*, double*, const int*);

int FORTRAN_FUNC(cher2k)(const char*, const char*, const int*, const int*, const float*, const float*, const int*, const float*, const int*, const float*, float*, const int*);
int FORTRAN_FUNC(zher2k)(const char*, const char*, const int*, const int*, const double*, const double*, const int*, const double*, const int*, const double*, double*, const int*);
int FORTRAN_FUNC(xher2k)(const char*, const char*, const int*, const int*, const double*, const double*, const int*, const double*, const int*, const double*, double*, const int*);
int FORTRAN_FUNC(cher2m)(const char*, const char*, const char*, const int*, const int*, const float*, const float*, const int*, const float*, const int*, const float*, float*, const int*);
int FORTRAN_FUNC(zher2m)(const char*, const char*, const char*, const int*, const int*, const double*, const double*, const int*, const double*, const int*, const double*, double*, const int*);
int FORTRAN_FUNC(xher2m)(const char*, const char*, const char*, const int*, const int*, const double*, const double*, const int*, const double*, const int*, const double*, double*, const int*);


void FORTRAN_FUNC(zgees)(const char*, const char*, void*, int*, complex_t*, int*, int*, complex_t*, complex_t*, int*, complex_t*, int*, double*, int*, int*);

void FORTRAN_FUNC(dgesdd)(const char*, int*, int*, double*, int*, double*, 
                          double*, int*, double*, int*, double*, 
                          int*, int*, int*);


#endif

