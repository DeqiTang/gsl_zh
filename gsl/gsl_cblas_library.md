# GSL CBLAS Library

The prototypes for the low-level CBLAS functions are declared in the file `gsl_cblas.h`. For the definition of the functions consult the documentation available from Netlib ([see BLAS References and Further Reading](https://www.gnu.org/software/gsl/doc/html/blas.html#sec-blas-references)).

## Level 1

- float `cblas_sdsdot`(const int *N*, const float *alpha*, const float * *x*, const int *incx*, const float * *y*, const int *incy*)

  

- double `cblas_dsdot`(const int *N*, const float * *x*, const int *incx*, const float * *y*, const int *incy*)

  

- float `cblas_sdot`(const int *N*, const float * *x*, const int *incx*, const float * *y*, const int *incy*)

  

- double `cblas_ddot`(const int *N*, const double * *x*, const int *incx*, const double * *y*, const int *incy*)

  

- void `cblas_cdotu_sub`(const int *N*, const void * *x*, const int *incx*, const void * *y*, const int *incy*, void * *dotu*)

  

- void `cblas_cdotc_sub`(const int *N*, const void * *x*, const int *incx*, const void * *y*, const int *incy*, void * *dotc*)

  

- void `cblas_zdotu_sub`(const int *N*, const void * *x*, const int *incx*, const void * *y*, const int *incy*, void * *dotu*)

  

- void `cblas_zdotc_sub`(const int *N*, const void * *x*, const int *incx*, const void * *y*, const int *incy*, void * *dotc*)

  

- float `cblas_snrm2`(const int *N*, const float * *x*, const int *incx*)

  

- float `cblas_sasum`(const int *N*, const float * *x*, const int *incx*)

  

- double `cblas_dnrm2`(const int *N*, const double * *x*, const int *incx*)

  

- double `cblas_dasum`(const int *N*, const double * *x*, const int *incx*)

  

- float `cblas_scnrm2`(const int *N*, const void * *x*, const int *incx*)

  

- float `cblas_scasum`(const int *N*, const void * *x*, const int *incx*)

  

- double `cblas_dznrm2`(const int *N*, const void * *x*, const int *incx*)

  

- double `cblas_dzasum`(const int *N*, const void * *x*, const int *incx*)

  

- CBLAS_INDEX `cblas_isamax`(const int *N*, const float * *x*, const int *incx*)

  

- CBLAS_INDEX `cblas_idamax`(const int *N*, const double * *x*, const int *incx*)

  

- CBLAS_INDEX `cblas_icamax`(const int *N*, const void * *x*, const int *incx*)

  

- CBLAS_INDEX `cblas_izamax`(const int *N*, const void * *x*, const int *incx*)

  

- void `cblas_sswap`(const int *N*, float * *x*, const int *incx*, float * *y*, const int *incy*)

  

- void `cblas_scopy`(const int *N*, const float * *x*, const int *incx*, float * *y*, const int *incy*)

  

- void `cblas_saxpy`(const int *N*, const float *alpha*, const float * *x*, const int *incx*, float * *y*, const int *incy*)

  

- void `cblas_dswap`(const int *N*, double * *x*, const int *incx*, double * *y*, const int *incy*)

  

- void `cblas_dcopy`(const int *N*, const double * *x*, const int *incx*, double * *y*, const int *incy*)

  

- void `cblas_daxpy`(const int *N*, const double *alpha*, const double * *x*, const int *incx*, double * *y*, const int *incy*)

  

- void `cblas_cswap`(const int *N*, void * *x*, const int *incx*, void * *y*, const int *incy*)

  

- void `cblas_ccopy`(const int *N*, const void * *x*, const int *incx*, void * *y*, const int *incy*)

  

- void `cblas_caxpy`(const int *N*, const void * *alpha*, const void * *x*, const int *incx*, void * *y*, const int *incy*)

  

- void `cblas_zswap`(const int *N*, void * *x*, const int *incx*, void * *y*, const int *incy*)

  

- void `cblas_zcopy`(const int *N*, const void * *x*, const int *incx*, void * *y*, const int *incy*)

  

- void `cblas_zaxpy`(const int *N*, const void * *alpha*, const void * *x*, const int *incx*, void * *y*, const int *incy*)

  

- void `cblas_srotg`(float * *a*, float * *b*, float * *c*, float * *s*)

  

- void `cblas_srotmg`(float * *d1*, float * *d2*, float * *b1*, const float *b2*, float * *P*)

  

- void `cblas_srot`(const int *N*, float * *x*, const int *incx*, float * *y*, const int *incy*, const float *c*, const float *s*)

  

- void `cblas_srotm`(const int *N*, float * *x*, const int *incx*, float * *y*, const int *incy*, const float * *P*)

  

- void `cblas_drotg`(double * *a*, double * *b*, double * *c*, double * *s*)

  

- void `cblas_drotmg`(double * *d1*, double * *d2*, double * *b1*, const double *b2*, double * *P*)

  

- void `cblas_drot`(const int *N*, double * *x*, const int *incx*, double * *y*, const int *incy*, const double *c*, const double *s*)

  

- void `cblas_drotm`(const int *N*, double * *x*, const int *incx*, double * *y*, const int *incy*, const double * *P*)

  

- void `cblas_sscal`(const int *N*, const float *alpha*, float * *x*, const int *incx*)

  

- void `cblas_dscal`(const int *N*, const double *alpha*, double * *x*, const int *incx*)

  

- void `cblas_cscal`(const int *N*, const void * *alpha*, void * *x*, const int *incx*)

  

- void `cblas_zscal`(const int *N*, const void * *alpha*, void * *x*, const int *incx*)

  

- void `cblas_csscal`(const int *N*, const float *alpha*, void * *x*, const int *incx*)

  

- void `cblas_zdscal`(const int *N*, const double *alpha*, void * *x*, const int *incx*)

  

## Level 2

- void `cblas_sgemv`(const enum CBLAS_ORDER *order*, const enum CBLAS_TRANSPOSE *TransA*, const int *M*, const int *N*, const float *alpha*, const float * *A*, const int *lda*, const float * *x*, const int *incx*, const float *beta*, float * *y*, const int *incy*)

  

- void `cblas_sgbmv`(const enum CBLAS_ORDER *order*, const enum CBLAS_TRANSPOSE *TransA*, const int *M*, const int *N*, const int *KL*, const int *KU*, const float *alpha*, const float * *A*, const int *lda*, const float * *x*, const int *incx*, const float *beta*, float * *y*, const int *incy*)

  

- void `cblas_strmv`(const enum CBLAS_ORDER *order*, const enum CBLAS_UPLO *Uplo*, const enum CBLAS_TRANSPOSE *TransA*, const enum CBLAS_DIAG *Diag*, const int *N*, const float * *A*, const int *lda*, float * *x*, const int *incx*)

  

- void `cblas_stbmv`(const enum CBLAS_ORDER *order*, const enum CBLAS_UPLO *Uplo*, const enum CBLAS_TRANSPOSE *TransA*, const enum CBLAS_DIAG *Diag*, const int *N*, const int *K*, const float * *A*, const int *lda*, float * *x*, const int *incx*)

  

- void `cblas_stpmv`(const enum CBLAS_ORDER *order*, const enum CBLAS_UPLO *Uplo*, const enum CBLAS_TRANSPOSE *TransA*, const enum CBLAS_DIAG *Diag*, const int *N*, const float * *Ap*, float * *x*, const int *incx*)

  

- void `cblas_strsv`(const enum CBLAS_ORDER *order*, const enum CBLAS_UPLO *Uplo*, const enum CBLAS_TRANSPOSE *TransA*, const enum CBLAS_DIAG *Diag*, const int *N*, const float * *A*, const int *lda*, float * *x*, const int *incx*)

  

- void `cblas_stbsv`(const enum CBLAS_ORDER *order*, const enum CBLAS_UPLO *Uplo*, const enum CBLAS_TRANSPOSE *TransA*, const enum CBLAS_DIAG *Diag*, const int *N*, const int *K*, const float * *A*, const int *lda*, float * *x*, const int *incx*)

  

- void `cblas_stpsv`(const enum CBLAS_ORDER *order*, const enum CBLAS_UPLO *Uplo*, const enum CBLAS_TRANSPOSE *TransA*, const enum CBLAS_DIAG *Diag*, const int *N*, const float * *Ap*, float * *x*, const int *incx*)

  

- void `cblas_dgemv`(const enum CBLAS_ORDER *order*, const enum CBLAS_TRANSPOSE *TransA*, const int *M*, const int *N*, const double *alpha*, const double * *A*, const int *lda*, const double * *x*, const int *incx*, const double *beta*, double * *y*, const int *incy*)

  

- void `cblas_dgbmv`(const enum CBLAS_ORDER *order*, const enum CBLAS_TRANSPOSE *TransA*, const int *M*, const int *N*, const int *KL*, const int *KU*, const double *alpha*, const double * *A*, const int *lda*, const double * *x*, const int *incx*, const double *beta*, double * *y*, const int *incy*)

  

- void `cblas_dtrmv`(const enum CBLAS_ORDER *order*, const enum CBLAS_UPLO *Uplo*, const enum CBLAS_TRANSPOSE *TransA*, const enum CBLAS_DIAG *Diag*, const int *N*, const double * *A*, const int *lda*, double * *x*, const int *incx*)

  

- void `cblas_dtbmv`(const enum CBLAS_ORDER *order*, const enum CBLAS_UPLO *Uplo*, const enum CBLAS_TRANSPOSE *TransA*, const enum CBLAS_DIAG *Diag*, const int *N*, const int *K*, const double * *A*, const int *lda*, double * *x*, const int *incx*)

  

- void `cblas_dtpmv`(const enum CBLAS_ORDER *order*, const enum CBLAS_UPLO *Uplo*, const enum CBLAS_TRANSPOSE *TransA*, const enum CBLAS_DIAG *Diag*, const int *N*, const double * *Ap*, double * *x*, const int *incx*)

  

- void `cblas_dtrsv`(const enum CBLAS_ORDER *order*, const enum CBLAS_UPLO *Uplo*, const enum CBLAS_TRANSPOSE *TransA*, const enum CBLAS_DIAG *Diag*, const int *N*, const double * *A*, const int *lda*, double * *x*, const int *incx*)

  

- void `cblas_dtbsv`(const enum CBLAS_ORDER *order*, const enum CBLAS_UPLO *Uplo*, const enum CBLAS_TRANSPOSE *TransA*, const enum CBLAS_DIAG *Diag*, const int *N*, const int *K*, const double * *A*, const int *lda*, double * *x*, const int *incx*)

  

- void `cblas_dtpsv`(const enum CBLAS_ORDER *order*, const enum CBLAS_UPLO *Uplo*, const enum CBLAS_TRANSPOSE *TransA*, const enum CBLAS_DIAG *Diag*, const int *N*, const double * *Ap*, double * *x*, const int *incx*)

  

- void `cblas_cgemv`(const enum CBLAS_ORDER *order*, const enum CBLAS_TRANSPOSE *TransA*, const int *M*, const int *N*, const void * *alpha*, const void * *A*, const int *lda*, const void * *x*, const int *incx*, const void * *beta*, void * *y*, const int *incy*)

  

- void `cblas_cgbmv`(const enum CBLAS_ORDER *order*, const enum CBLAS_TRANSPOSE *TransA*, const int *M*, const int *N*, const int *KL*, const int *KU*, const void * *alpha*, const void * *A*, const int *lda*, const void * *x*, const int *incx*, const void * *beta*, void * *y*, const int *incy*)

  

- void `cblas_ctrmv`(const enum CBLAS_ORDER *order*, const enum CBLAS_UPLO *Uplo*, const enum CBLAS_TRANSPOSE *TransA*, const enum CBLAS_DIAG *Diag*, const int *N*, const void * *A*, const int *lda*, void * *x*, const int *incx*)

  

- void `cblas_ctbmv`(const enum CBLAS_ORDER *order*, const enum CBLAS_UPLO *Uplo*, const enum CBLAS_TRANSPOSE *TransA*, const enum CBLAS_DIAG *Diag*, const int *N*, const int *K*, const void * *A*, const int *lda*, void * *x*, const int *incx*)

  

- void `cblas_ctpmv`(const enum CBLAS_ORDER *order*, const enum CBLAS_UPLO *Uplo*, const enum CBLAS_TRANSPOSE *TransA*, const enum CBLAS_DIAG *Diag*, const int *N*, const void * *Ap*, void * *x*, const int *incx*)

  

- void `cblas_ctrsv`(const enum CBLAS_ORDER *order*, const enum CBLAS_UPLO *Uplo*, const enum CBLAS_TRANSPOSE *TransA*, const enum CBLAS_DIAG *Diag*, const int *N*, const void * *A*, const int *lda*, void * *x*, const int *incx*)

  

- void `cblas_ctbsv`(const enum CBLAS_ORDER *order*, const enum CBLAS_UPLO *Uplo*, const enum CBLAS_TRANSPOSE *TransA*, const enum CBLAS_DIAG *Diag*, const int *N*, const int *K*, const void * *A*, const int *lda*, void * *x*, const int *incx*)

  

- void `cblas_ctpsv`(const enum CBLAS_ORDER *order*, const enum CBLAS_UPLO *Uplo*, const enum CBLAS_TRANSPOSE *TransA*, const enum CBLAS_DIAG *Diag*, const int *N*, const void * *Ap*, void * *x*, const int *incx*)

  

- void `cblas_zgemv`(const enum CBLAS_ORDER *order*, const enum CBLAS_TRANSPOSE *TransA*, const int *M*, const int *N*, const void * *alpha*, const void * *A*, const int *lda*, const void * *x*, const int *incx*, const void * *beta*, void * *y*, const int *incy*)

  

- void `cblas_zgbmv`(const enum CBLAS_ORDER *order*, const enum CBLAS_TRANSPOSE *TransA*, const int *M*, const int *N*, const int *KL*, const int *KU*, const void * *alpha*, const void * *A*, const int *lda*, const void * *x*, const int *incx*, const void * *beta*, void * *y*, const int *incy*)

  

- void `cblas_ztrmv`(const enum CBLAS_ORDER *order*, const enum CBLAS_UPLO *Uplo*, const enum CBLAS_TRANSPOSE *TransA*, const enum CBLAS_DIAG *Diag*, const int *N*, const void * *A*, const int *lda*, void * *x*, const int *incx*)

  

- void `cblas_ztbmv`(const enum CBLAS_ORDER *order*, const enum CBLAS_UPLO *Uplo*, const enum CBLAS_TRANSPOSE *TransA*, const enum CBLAS_DIAG *Diag*, const int *N*, const int *K*, const void * *A*, const int *lda*, void * *x*, const int *incx*)

  

- void `cblas_ztpmv`(const enum CBLAS_ORDER *order*, const enum CBLAS_UPLO *Uplo*, const enum CBLAS_TRANSPOSE *TransA*, const enum CBLAS_DIAG *Diag*, const int *N*, const void * *Ap*, void * *x*, const int *incx*)

  

- void `cblas_ztrsv`(const enum CBLAS_ORDER *order*, const enum CBLAS_UPLO *Uplo*, const enum CBLAS_TRANSPOSE *TransA*, const enum CBLAS_DIAG *Diag*, const int *N*, const void * *A*, const int *lda*, void * *x*, const int *incx*)

  

- void `cblas_ztbsv`(const enum CBLAS_ORDER *order*, const enum CBLAS_UPLO *Uplo*, const enum CBLAS_TRANSPOSE *TransA*, const enum CBLAS_DIAG *Diag*, const int *N*, const int *K*, const void * *A*, const int *lda*, void * *x*, const int *incx*)

  

- void `cblas_ztpsv`(const enum CBLAS_ORDER *order*, const enum CBLAS_UPLO *Uplo*, const enum CBLAS_TRANSPOSE *TransA*, const enum CBLAS_DIAG *Diag*, const int *N*, const void * *Ap*, void * *x*, const int *incx*)

  

- void `cblas_ssymv`(const enum CBLAS_ORDER *order*, const enum CBLAS_UPLO *Uplo*, const int *N*, const float *alpha*, const float * *A*, const int *lda*, const float * *x*, const int *incx*, const float *beta*, float * *y*, const int *incy*)

  

- void `cblas_ssbmv`(const enum CBLAS_ORDER *order*, const enum CBLAS_UPLO *Uplo*, const int *N*, const int *K*, const float *alpha*, const float * *A*, const int *lda*, const float * *x*, const int *incx*, const float *beta*, float * *y*, const int *incy*)

  

- void `cblas_sspmv`(const enum CBLAS_ORDER *order*, const enum CBLAS_UPLO *Uplo*, const int *N*, const float *alpha*, const float * *Ap*, const float * *x*, const int *incx*, const float *beta*, float * *y*, const int *incy*)

  

- void `cblas_sger`(const enum CBLAS_ORDER *order*, const int *M*, const int *N*, const float *alpha*, const float * *x*, const int *incx*, const float * *y*, const int *incy*, float * *A*, const int *lda*)

  

- void `cblas_ssyr`(const enum CBLAS_ORDER *order*, const enum CBLAS_UPLO *Uplo*, const int *N*, const float *alpha*, const float * *x*, const int *incx*, float * *A*, const int *lda*)

  

- void `cblas_sspr`(const enum CBLAS_ORDER *order*, const enum CBLAS_UPLO *Uplo*, const int *N*, const float *alpha*, const float * *x*, const int *incx*, float * *Ap*)

  

- void `cblas_ssyr2`(const enum CBLAS_ORDER *order*, const enum CBLAS_UPLO *Uplo*, const int *N*, const float *alpha*, const float * *x*, const int *incx*, const float * *y*, const int *incy*, float * *A*, const int *lda*)

  

- void `cblas_sspr2`(const enum CBLAS_ORDER *order*, const enum CBLAS_UPLO *Uplo*, const int *N*, const float *alpha*, const float * *x*, const int *incx*, const float * *y*, const int *incy*, float * *A*)

  

- void `cblas_dsymv`(const enum CBLAS_ORDER *order*, const enum CBLAS_UPLO *Uplo*, const int *N*, const double *alpha*, const double * *A*, const int *lda*, const double * *x*, const int *incx*, const double *beta*, double * *y*, const int *incy*)

  

- void `cblas_dsbmv`(const enum CBLAS_ORDER *order*, const enum CBLAS_UPLO *Uplo*, const int *N*, const int *K*, const double *alpha*, const double * *A*, const int *lda*, const double * *x*, const int *incx*, const double *beta*, double * *y*, const int *incy*)

  

- void `cblas_dspmv`(const enum CBLAS_ORDER *order*, const enum CBLAS_UPLO *Uplo*, const int *N*, const double *alpha*, const double * *Ap*, const double * *x*, const int *incx*, const double *beta*, double * *y*, const int *incy*)

  

- void `cblas_dger`(const enum CBLAS_ORDER *order*, const int *M*, const int *N*, const double *alpha*, const double * *x*, const int *incx*, const double * *y*, const int *incy*, double * *A*, const int *lda*)

  

- void `cblas_dsyr`(const enum CBLAS_ORDER *order*, const enum CBLAS_UPLO *Uplo*, const int *N*, const double *alpha*, const double * *x*, const int *incx*, double * *A*, const int *lda*)

  

- void `cblas_dspr`(const enum CBLAS_ORDER *order*, const enum CBLAS_UPLO *Uplo*, const int *N*, const double *alpha*, const double * *x*, const int *incx*, double * *Ap*)

  

- void `cblas_dsyr2`(const enum CBLAS_ORDER *order*, const enum CBLAS_UPLO *Uplo*, const int *N*, const double *alpha*, const double * *x*, const int *incx*, const double * *y*, const int *incy*, double * *A*, const int *lda*)

  

- void `cblas_dspr2`(const enum CBLAS_ORDER *order*, const enum CBLAS_UPLO *Uplo*, const int *N*, const double *alpha*, const double * *x*, const int *incx*, const double * *y*, const int *incy*, double * *A*)

  

- void `cblas_chemv`(const enum CBLAS_ORDER *order*, const enum CBLAS_UPLO *Uplo*, const int *N*, const void * *alpha*, const void * *A*, const int *lda*, const void * *x*, const int *incx*, const void * *beta*, void * *y*, const int *incy*)

  

- void `cblas_chbmv`(const enum CBLAS_ORDER *order*, const enum CBLAS_UPLO *Uplo*, const int *N*, const int *K*, const void * *alpha*, const void * *A*, const int *lda*, const void * *x*, const int *incx*, const void * *beta*, void * *y*, const int *incy*)

  

- void `cblas_chpmv`(const enum CBLAS_ORDER *order*, const enum CBLAS_UPLO *Uplo*, const int *N*, const void * *alpha*, const void * *Ap*, const void * *x*, const int *incx*, const void * *beta*, void * *y*, const int *incy*)

  

- void `cblas_cgeru`(const enum CBLAS_ORDER *order*, const int *M*, const int *N*, const void * *alpha*, const void * *x*, const int *incx*, const void * *y*, const int *incy*, void * *A*, const int *lda*)

  

- void `cblas_cgerc`(const enum CBLAS_ORDER *order*, const int *M*, const int *N*, const void * *alpha*, const void * *x*, const int *incx*, const void * *y*, const int *incy*, void * *A*, const int *lda*)

  

- void `cblas_cher`(const enum CBLAS_ORDER *order*, const enum CBLAS_UPLO *Uplo*, const int *N*, const float *alpha*, const void * *x*, const int *incx*, void * *A*, const int *lda*)

  

- void `cblas_chpr`(const enum CBLAS_ORDER *order*, const enum CBLAS_UPLO *Uplo*, const int *N*, const float *alpha*, const void * *x*, const int *incx*, void * *A*)

  

- void `cblas_cher2`(const enum CBLAS_ORDER *order*, const enum CBLAS_UPLO *Uplo*, const int *N*, const void * *alpha*, const void * *x*, const int *incx*, const void * *y*, const int *incy*, void * *A*, const int *lda*)

  

- void `cblas_chpr2`(const enum CBLAS_ORDER *order*, const enum CBLAS_UPLO *Uplo*, const int *N*, const void * *alpha*, const void * *x*, const int *incx*, const void * *y*, const int *incy*, void * *Ap*)

  

- void `cblas_zhemv`(const enum CBLAS_ORDER *order*, const enum CBLAS_UPLO *Uplo*, const int *N*, const void * *alpha*, const void * *A*, const int *lda*, const void * *x*, const int *incx*, const void * *beta*, void * *y*, const int *incy*)

  

- void `cblas_zhbmv`(const enum CBLAS_ORDER *order*, const enum CBLAS_UPLO *Uplo*, const int *N*, const int *K*, const void * *alpha*, const void * *A*, const int *lda*, const void * *x*, const int *incx*, const void * *beta*, void * *y*, const int *incy*)

  

- void `cblas_zhpmv`(const enum CBLAS_ORDER *order*, const enum CBLAS_UPLO *Uplo*, const int *N*, const void * *alpha*, const void * *Ap*, const void * *x*, const int *incx*, const void * *beta*, void * *y*, const int *incy*)

  

- void `cblas_zgeru`(const enum CBLAS_ORDER *order*, const int *M*, const int *N*, const void * *alpha*, const void * *x*, const int *incx*, const void * *y*, const int *incy*, void * *A*, const int *lda*)

  

- void `cblas_zgerc`(const enum CBLAS_ORDER *order*, const int *M*, const int *N*, const void * *alpha*, const void * *x*, const int *incx*, const void * *y*, const int *incy*, void * *A*, const int *lda*)

  

- void `cblas_zher`(const enum CBLAS_ORDER *order*, const enum CBLAS_UPLO *Uplo*, const int *N*, const double *alpha*, const void * *x*, const int *incx*, void * *A*, const int *lda*)

  

- void `cblas_zhpr`(const enum CBLAS_ORDER *order*, const enum CBLAS_UPLO *Uplo*, const int *N*, const double *alpha*, const void * *x*, const int *incx*, void * *A*)

  

- void `cblas_zher2`(const enum CBLAS_ORDER *order*, const enum CBLAS_UPLO *Uplo*, const int *N*, const void * *alpha*, const void * *x*, const int *incx*, const void * *y*, const int *incy*, void * *A*, const int *lda*)

  

- void `cblas_zhpr2`(const enum CBLAS_ORDER *order*, const enum CBLAS_UPLO *Uplo*, const int *N*, const void * *alpha*, const void * *x*, const int *incx*, const void * *y*, const int *incy*, void * *Ap*)

  

## Level 3

- void `cblas_sgemm`(const enum CBLAS_ORDER *Order*, const enum CBLAS_TRANSPOSE *TransA*, const enum CBLAS_TRANSPOSE *TransB*, const int *M*, const int *N*, const int *K*, const float *alpha*, const float * *A*, const int *lda*, const float * *B*, const int *ldb*, const float *beta*, float * *C*, const int *ldc*)

  

- void `cblas_ssymm`(const enum CBLAS_ORDER *Order*, const enum CBLAS_SIDE *Side*, const enum CBLAS_UPLO *Uplo*, const int *M*, const int *N*, const float *alpha*, const float * *A*, const int *lda*, const float * *B*, const int *ldb*, const float *beta*, float * *C*, const int *ldc*)

  

- void `cblas_ssyrk`(const enum CBLAS_ORDER *Order*, const enum CBLAS_UPLO *Uplo*, const enum CBLAS_TRANSPOSE *Trans*, const int *N*, const int *K*, const float *alpha*, const float * *A*, const int *lda*, const float *beta*, float * *C*, const int *ldc*)

  

- void `cblas_ssyr2k`(const enum CBLAS_ORDER *Order*, const enum CBLAS_UPLO *Uplo*, const enum CBLAS_TRANSPOSE *Trans*, const int *N*, const int *K*, const float *alpha*, const float * *A*, const int *lda*, const float * *B*, const int *ldb*, const float *beta*, float * *C*, const int *ldc*)

  

- void `cblas_strmm`(const enum CBLAS_ORDER *Order*, const enum CBLAS_SIDE *Side*, const enum CBLAS_UPLO *Uplo*, const enum CBLAS_TRANSPOSE *TransA*, const enum CBLAS_DIAG *Diag*, const int *M*, const int *N*, const float *alpha*, const float * *A*, const int *lda*, float * *B*, const int *ldb*)

  

- void `cblas_strsm`(const enum CBLAS_ORDER *Order*, const enum CBLAS_SIDE *Side*, const enum CBLAS_UPLO *Uplo*, const enum CBLAS_TRANSPOSE *TransA*, const enum CBLAS_DIAG *Diag*, const int *M*, const int *N*, const float *alpha*, const float * *A*, const int *lda*, float * *B*, const int *ldb*)

  

- void `cblas_dgemm`(const enum CBLAS_ORDER *Order*, const enum CBLAS_TRANSPOSE *TransA*, const enum CBLAS_TRANSPOSE *TransB*, const int *M*, const int *N*, const int *K*, const double *alpha*, const double * *A*, const int *lda*, const double * *B*, const int *ldb*, const double *beta*, double * *C*, const int *ldc*)

  

- void `cblas_dsymm`(const enum CBLAS_ORDER *Order*, const enum CBLAS_SIDE *Side*, const enum CBLAS_UPLO *Uplo*, const int *M*, const int *N*, const double *alpha*, const double * *A*, const int *lda*, const double * *B*, const int *ldb*, const double *beta*, double * *C*, const int *ldc*)

  

- void `cblas_dsyrk`(const enum CBLAS_ORDER *Order*, const enum CBLAS_UPLO *Uplo*, const enum CBLAS_TRANSPOSE *Trans*, const int *N*, const int *K*, const double *alpha*, const double * *A*, const int *lda*, const double *beta*, double * *C*, const int *ldc*)

  

- void `cblas_dsyr2k`(const enum CBLAS_ORDER *Order*, const enum CBLAS_UPLO *Uplo*, const enum CBLAS_TRANSPOSE *Trans*, const int *N*, const int *K*, const double *alpha*, const double * *A*, const int *lda*, const double * *B*, const int *ldb*, const double *beta*, double * *C*, const int *ldc*)

  

- void `cblas_dtrmm`(const enum CBLAS_ORDER *Order*, const enum CBLAS_SIDE *Side*, const enum CBLAS_UPLO *Uplo*, const enum CBLAS_TRANSPOSE *TransA*, const enum CBLAS_DIAG *Diag*, const int *M*, const int *N*, const double *alpha*, const double * *A*, const int *lda*, double * *B*, const int *ldb*)

  

- void `cblas_dtrsm`(const enum CBLAS_ORDER *Order*, const enum CBLAS_SIDE *Side*, const enum CBLAS_UPLO *Uplo*, const enum CBLAS_TRANSPOSE *TransA*, const enum CBLAS_DIAG *Diag*, const int *M*, const int *N*, const double *alpha*, const double * *A*, const int *lda*, double * *B*, const int *ldb*)

  

- void `cblas_cgemm`(const enum CBLAS_ORDER *Order*, const enum CBLAS_TRANSPOSE *TransA*, const enum CBLAS_TRANSPOSE *TransB*, const int *M*, const int *N*, const int *K*, const void * *alpha*, const void * *A*, const int *lda*, const void * *B*, const int *ldb*, const void * *beta*, void * *C*, const int *ldc*)

  

- void `cblas_csymm`(const enum CBLAS_ORDER *Order*, const enum CBLAS_SIDE *Side*, const enum CBLAS_UPLO *Uplo*, const int *M*, const int *N*, const void * *alpha*, const void * *A*, const int *lda*, const void * *B*, const int *ldb*, const void * *beta*, void * *C*, const int *ldc*)

  

- void `cblas_csyrk`(const enum CBLAS_ORDER *Order*, const enum CBLAS_UPLO *Uplo*, const enum CBLAS_TRANSPOSE *Trans*, const int *N*, const int *K*, const void * *alpha*, const void * *A*, const int *lda*, const void * *beta*, void * *C*, const int *ldc*)

  

- void `cblas_csyr2k`(const enum CBLAS_ORDER *Order*, const enum CBLAS_UPLO *Uplo*, const enum CBLAS_TRANSPOSE *Trans*, const int *N*, const int *K*, const void * *alpha*, const void * *A*, const int *lda*, const void * *B*, const int *ldb*, const void * *beta*, void * *C*, const int *ldc*)

  

- void `cblas_ctrmm`(const enum CBLAS_ORDER *Order*, const enum CBLAS_SIDE *Side*, const enum CBLAS_UPLO *Uplo*, const enum CBLAS_TRANSPOSE *TransA*, const enum CBLAS_DIAG *Diag*, const int *M*, const int *N*, const void * *alpha*, const void * *A*, const int *lda*, void * *B*, const int *ldb*)

  

- void `cblas_ctrsm`(const enum CBLAS_ORDER *Order*, const enum CBLAS_SIDE *Side*, const enum CBLAS_UPLO *Uplo*, const enum CBLAS_TRANSPOSE *TransA*, const enum CBLAS_DIAG *Diag*, const int *M*, const int *N*, const void * *alpha*, const void * *A*, const int *lda*, void * *B*, const int *ldb*)

  

- void `cblas_zgemm`(const enum CBLAS_ORDER *Order*, const enum CBLAS_TRANSPOSE *TransA*, const enum CBLAS_TRANSPOSE *TransB*, const int *M*, const int *N*, const int *K*, const void * *alpha*, const void * *A*, const int *lda*, const void * *B*, const int *ldb*, const void * *beta*, void * *C*, const int *ldc*)

  

- void `cblas_zsymm`(const enum CBLAS_ORDER *Order*, const enum CBLAS_SIDE *Side*, const enum CBLAS_UPLO *Uplo*, const int *M*, const int *N*, const void * *alpha*, const void * *A*, const int *lda*, const void * *B*, const int *ldb*, const void * *beta*, void * *C*, const int *ldc*)

  

- void `cblas_zsyrk`(const enum CBLAS_ORDER *Order*, const enum CBLAS_UPLO *Uplo*, const enum CBLAS_TRANSPOSE *Trans*, const int *N*, const int *K*, const void * *alpha*, const void * *A*, const int *lda*, const void * *beta*, void * *C*, const int *ldc*)

  

- void `cblas_zsyr2k`(const enum CBLAS_ORDER *Order*, const enum CBLAS_UPLO *Uplo*, const enum CBLAS_TRANSPOSE *Trans*, const int *N*, const int *K*, const void * *alpha*, const void * *A*, const int *lda*, const void * *B*, const int *ldb*, const void * *beta*, void * *C*, const int *ldc*)

  

- void `cblas_ztrmm`(const enum CBLAS_ORDER *Order*, const enum CBLAS_SIDE *Side*, const enum CBLAS_UPLO *Uplo*, const enum CBLAS_TRANSPOSE *TransA*, const enum CBLAS_DIAG *Diag*, const int *M*, const int *N*, const void * *alpha*, const void * *A*, const int *lda*, void * *B*, const int *ldb*)

  

- void `cblas_ztrsm`(const enum CBLAS_ORDER *Order*, const enum CBLAS_SIDE *Side*, const enum CBLAS_UPLO *Uplo*, const enum CBLAS_TRANSPOSE *TransA*, const enum CBLAS_DIAG *Diag*, const int *M*, const int *N*, const void * *alpha*, const void * *A*, const int *lda*, void * *B*, const int *ldb*)

  

- void `cblas_chemm`(const enum CBLAS_ORDER *Order*, const enum CBLAS_SIDE *Side*, const enum CBLAS_UPLO *Uplo*, const int *M*, const int *N*, const void * *alpha*, const void * *A*, const int *lda*, const void * *B*, const int *ldb*, const void * *beta*, void * *C*, const int *ldc*)

  

- void `cblas_cherk`(const enum CBLAS_ORDER *Order*, const enum CBLAS_UPLO *Uplo*, const enum CBLAS_TRANSPOSE *Trans*, const int *N*, const int *K*, const float *alpha*, const void * *A*, const int *lda*, const float *beta*, void * *C*, const int *ldc*)

  

- void `cblas_cher2k`(const enum CBLAS_ORDER *Order*, const enum CBLAS_UPLO *Uplo*, const enum CBLAS_TRANSPOSE *Trans*, const int *N*, const int *K*, const void * *alpha*, const void * *A*, const int *lda*, const void * *B*, const int *ldb*, const float *beta*, void * *C*, const int *ldc*)

  

- void `cblas_zhemm`(const enum CBLAS_ORDER *Order*, const enum CBLAS_SIDE *Side*, const enum CBLAS_UPLO *Uplo*, const int *M*, const int *N*, const void * *alpha*, const void * *A*, const int *lda*, const void * *B*, const int *ldb*, const void * *beta*, void * *C*, const int *ldc*)

  

- void `cblas_zherk`(const enum CBLAS_ORDER *Order*, const enum CBLAS_UPLO *Uplo*, const enum CBLAS_TRANSPOSE *Trans*, const int *N*, const int *K*, const double *alpha*, const void * *A*, const int *lda*, const double *beta*, void * *C*, const int *ldc*)

  

- void `cblas_zher2k`(const enum CBLAS_ORDER *Order*, const enum CBLAS_UPLO *Uplo*, const enum CBLAS_TRANSPOSE *Trans*, const int *N*, const int *K*, const void * *alpha*, const void * *A*, const int *lda*, const void * *B*, const int *ldb*, const double *beta*, void * *C*, const int *ldc*)

  

- void `cblas_xerbla`(int *p*, const char * *rout*, const char * *form*, ...)

  

## Examples

The following program computes the product of two matrices using the Level-3 BLAS function SGEMM,

![\left(   \begin{array}{ccc}     0.11 & 0.12 & 0.13\\     0.21 & 0.22 & 0.23   \end{array} \right) \left(   \begin{array}{cc}     1011 & 1012\\     1021 & 1022\\     1031 & 1032   \end{array} \right) = \left(   \begin{array}{cc}     367.76 & 368.12\\     674.06 & 674.72   \end{array} \right)](https://www.gnu.org/software/gsl/doc/html/_images/math/46b30ee68fdca1cd520b5b4149cf80ea6cf2b0dc.png)

The matrices are stored in row major order but could be stored in column major order if the first argument of the call to [`cblas_sgemm()`](https://www.gnu.org/software/gsl/doc/html/cblas.html#c.cblas_sgemm) was changed to `CblasColMajor`.

```
#include <stdio.h>
#include <gsl/gsl_cblas.h>

int
main (void)
{
  int lda = 3;

  float A[] = { 0.11, 0.12, 0.13,
                0.21, 0.22, 0.23 };

  int ldb = 2;

  float B[] = { 1011, 1012,
                1021, 1022,
                1031, 1032 };

  int ldc = 2;

  float C[] = { 0.00, 0.00,
                0.00, 0.00 };

  /* Compute C = A B */

  cblas_sgemm (CblasRowMajor,
               CblasNoTrans, CblasNoTrans, 2, 2, 3,
               1.0, A, lda, B, ldb, 0.0, C, ldc);

  printf ("[ %g, %g\n", C[0], C[1]);
  printf ("  %g, %g ]\n", C[2], C[3]);

  return 0;
}
```

To compile the program use the following command line:

```
$ gcc -Wall demo.c -lgslcblas
```

There is no need to link with the main library `-lgsl` in this case as the CBLAS library is an independent unit. Here is the output from the program,

```
[ 367.76, 368.12
  674.06, 674.72 ]
```