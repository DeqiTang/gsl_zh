# BLAS Support

The Basic Linear Algebra Subprograms (BLAS) define a set of fundamental operations on vectors and matrices which can be used to create optimized higher-level linear algebra functionality.

The library provides a low-level layer which corresponds directly to the C-language BLAS standard, referred to here as “CBLAS”, and a higher-level interface for operations on GSL vectors and matrices. Users who are interested in simple operations on GSL vector and matrix objects should use the high-level layer described in this chapter. The functions are declared in the file `gsl_blas.h`and should satisfy the needs of most users.

Note that GSL matrices are implemented using dense-storage so the interface only includes the corresponding dense-storage BLAS functions. The full BLAS functionality for band-format and packed-format matrices is available through the low-level CBLAS interface. Similarly, GSL vectors are restricted to positive strides, whereas the low-level CBLAS interface supports negative strides as specified in the BLAS standard [[1\]](https://www.gnu.org/software/gsl/doc/html/blas.html#f1).

The interface for the `gsl_cblas` layer is specified in the file `gsl_cblas.h`. This interface corresponds to the BLAS Technical Forum’s standard for the C interface to legacy BLAS implementations. Users who have access to other conforming CBLAS implementations can use these in place of the version provided by the library. Note that users who have only a Fortran BLAS library can use a CBLAS conformant wrapper to convert it into a CBLAS library. A reference CBLAS wrapper for legacy Fortran implementations exists as part of the CBLAS standard and can be obtained from Netlib. The complete set of CBLAS functions is listed in an [appendix](https://www.gnu.org/software/gsl/doc/html/cblas.html#chap-cblas).

There are three levels of BLAS operations,

| **Level 1** | Vector operations, e.g. ![y = \alpha x + y](https://www.gnu.org/software/gsl/doc/html/_images/math/bf4c0bf09422b831ce1c179126efb61a83cc6bb7.png) |
| ----------- | ------------------------------------------------------------ |
| **Level 2** | Matrix-vector operations, e.g. ![y = \alpha A x + \beta y](https://www.gnu.org/software/gsl/doc/html/_images/math/624871f488f273648511a3e51c9fd2ceb3c4f362.png) |
| **Level 3** | Matrix-matrix operations, e.g. ![C = \alpha A B + C](https://www.gnu.org/software/gsl/doc/html/_images/math/f416e48e9915bc751d298e7de209cd0370dd8a7f.png) |

Each routine has a name which specifies the operation, the type of matrices involved and their precisions. Some of the most common operations and their names are given below,

| **DOT**  | scalar product, ![x^T y](https://www.gnu.org/software/gsl/doc/html/_images/math/e872907b8ee554d41f9375630dff5c9a4de74524.png) |
| -------- | ------------------------------------------------------------ |
| **AXPY** | vector sum, ![\alpha x + y](https://www.gnu.org/software/gsl/doc/html/_images/math/d48df2ffd4af5d676901723fe3843dbb3c01eead.png) |
| **MV**   | matrix-vector product, ![A x](https://www.gnu.org/software/gsl/doc/html/_images/math/5fc2bb820b35f1db304cfb47fab3deb258274405.png) |
| **SV**   | matrix-vector solve, ![inv(A) x](https://www.gnu.org/software/gsl/doc/html/_images/math/c0b168a003f9d5b3bd47ed613a306e2c9f60c9a7.png) |
| **MM**   | matrix-matrix product, ![A B](https://www.gnu.org/software/gsl/doc/html/_images/math/818ae2ad823347cb78411eacd060db363a1f585e.png) |
| **SM**   | matrix-matrix solve, ![inv(A) B](https://www.gnu.org/software/gsl/doc/html/_images/math/563e7599b23cae3cb657e96b38bbb9033d1609d6.png) |

The types of matrices are,

| **GE** | general           |
| ------ | ----------------- |
| **GB** | general band      |
| **SY** | symmetric         |
| **SB** | symmetric band    |
| **SP** | symmetric packed  |
| **HE** | hermitian         |
| **HB** | hermitian band    |
| **HP** | hermitian packed  |
| **TR** | triangular        |
| **TB** | triangular band   |
| **TP** | triangular packed |

Each operation is defined for four precisions,

| **S** | single real    |
| ----- | -------------- |
| **D** | double real    |
| **C** | single complex |
| **Z** | double complex |

Thus, for example, the name SGEMM stands for “single-precision general matrix-matrix multiply” and ZGEMM stands for “double-precision complex matrix-matrix multiply”.

Note that the vector and matrix arguments to BLAS functions must not be aliased, as the results are undefined when the underlying arrays overlap ([Aliasing of arrays](https://www.gnu.org/software/gsl/doc/html/usage.html#aliasing-of-arrays)).

## GSL BLAS Interface

GSL provides dense vector and matrix objects, based on the relevant built-in types. The library provides an interface to the BLAS operations which apply to these objects. The interface to this functionality is given in the file `gsl_blas.h`.

### Level 1



- int `gsl_blas_sdsdot`(float *alpha*, const gsl_vector_float * *x*, const gsl_vector_float * *y*, float * *result*)

  This function computes the sum ![\alpha + x^T y](https://www.gnu.org/software/gsl/doc/html/_images/math/be907e2c347edc891f15185a2549b3ff6d6d4f7e.png) for the vectors `x` and `y`, returning the result in `result`.

- int `gsl_blas_sdot`(const gsl_vector_float * *x*, const gsl_vector_float * *y*, float * *result*)

- int `gsl_blas_dsdot`(const gsl_vector_float * *x*, const gsl_vector_float * *y*, double * *result*)

- int `gsl_blas_ddot`(const [gsl_vector](https://www.gnu.org/software/gsl/doc/html/vectors.html#c.gsl_vector) * *x*, const [gsl_vector](https://www.gnu.org/software/gsl/doc/html/vectors.html#c.gsl_vector) * *y*, double * *result*)

  These functions compute the scalar product ![x^T y](https://www.gnu.org/software/gsl/doc/html/_images/math/e872907b8ee554d41f9375630dff5c9a4de74524.png) for the vectors `x` and `y`, returning the result in `result`.

- int `gsl_blas_cdotu`(const gsl_vector_complex_float * *x*, const gsl_vector_complex_float * *y*, gsl_complex_float * *dotu*)

- int `gsl_blas_zdotu`(const gsl_vector_complex * *x*, const gsl_vector_complex * *y*, gsl_complex * *dotu*)

  These functions compute the complex scalar product ![x^T y](https://www.gnu.org/software/gsl/doc/html/_images/math/e872907b8ee554d41f9375630dff5c9a4de74524.png) for the vectors `x` and `y`, returning the result in `dotu`

- int `gsl_blas_cdotc`(const gsl_vector_complex_float * *x*, const gsl_vector_complex_float * *y*, gsl_complex_float * *dotc*)

- int `gsl_blas_zdotc`(const gsl_vector_complex * *x*, const gsl_vector_complex * *y*, gsl_complex * *dotc*)

  These functions compute the complex conjugate scalar product ![x^H y](https://www.gnu.org/software/gsl/doc/html/_images/math/36e5c098c5123a44a8cf5294aa4f9a6bd25eb2ec.png) for the vectors `x` and `y`, returning the result in `dotc`



- float `gsl_blas_snrm2`(const gsl_vector_float * *x*)

- double `gsl_blas_dnrm2`(const [gsl_vector](https://www.gnu.org/software/gsl/doc/html/vectors.html#c.gsl_vector) * *x*)

  These functions compute the Euclidean norm ![||x||_2 = \sqrt{\sum x_i^2}](https://www.gnu.org/software/gsl/doc/html/_images/math/4dca4f3fc2a45ba778cec427713e5d7909666334.png) of the vector `x`.

- float `gsl_blas_scnrm2`(const gsl_vector_complex_float * *x*)

- double `gsl_blas_dznrm2`(const gsl_vector_complex * *x*)

  These functions compute the Euclidean norm of the complex vector `x`,![||x||_2 = \sqrt{\sum (\Re(x_i)^2 + \Im(x_i)^2)}.](https://www.gnu.org/software/gsl/doc/html/_images/math/b0aa7b9e8d8cd4e9d2edd46982cc69d2421cd496.png)



- float `gsl_blas_sasum`(const gsl_vector_float * *x*)

- double `gsl_blas_dasum`(const [gsl_vector](https://www.gnu.org/software/gsl/doc/html/vectors.html#c.gsl_vector) * *x*)

  These functions compute the absolute sum ![\sum |x_i|](https://www.gnu.org/software/gsl/doc/html/_images/math/7f5a721479f3559dfc62f991acb564069dbd6053.png) of the elements of the vector `x`.

- float `gsl_blas_scasum`(const gsl_vector_complex_float * *x*)

- double `gsl_blas_dzasum`(const gsl_vector_complex * *x*)

  These functions compute the sum of the magnitudes of the real and imaginary parts of the complex vector `x`, ![\sum \left( |\Re(x_i)| + |\Im(x_i)| \right)](https://www.gnu.org/software/gsl/doc/html/_images/math/f2acf458b95f05810a941efd9b2b4e9e3849214a.png).



- CBLAS_INDEX_t `gsl_blas_isamax`(const gsl_vector_float * *x*)

- CBLAS_INDEX_t `gsl_blas_idamax`(const [gsl_vector](https://www.gnu.org/software/gsl/doc/html/vectors.html#c.gsl_vector) * *x*)

- CBLAS_INDEX_t `gsl_blas_icamax`(const gsl_vector_complex_float * *x*)

- CBLAS_INDEX_t `gsl_blas_izamax`(const gsl_vector_complex * *x*)

  These functions return the index of the largest element of the vector `x`. The largest element is determined by its absolute magnitude for real vectors and by the sum of the magnitudes of the real and imaginary parts ![|\Re(x_i)| + |\Im(x_i)|](https://www.gnu.org/software/gsl/doc/html/_images/math/5be8fe92f25d71bc5cbfaed8544feed6f9c7a6c8.png) for complex vectors. If the largest value occurs several times then the index of the first occurrence is returned.



- int `gsl_blas_sswap`(gsl_vector_float * *x*, gsl_vector_float * *y*)

- int `gsl_blas_dswap`([gsl_vector](https://www.gnu.org/software/gsl/doc/html/vectors.html#c.gsl_vector) * *x*, [gsl_vector](https://www.gnu.org/software/gsl/doc/html/vectors.html#c.gsl_vector) * *y*)

- int `gsl_blas_cswap`(gsl_vector_complex_float * *x*, gsl_vector_complex_float * *y*)

- int `gsl_blas_zswap`(gsl_vector_complex * *x*, gsl_vector_complex * *y*)

  These functions exchange the elements of the vectors `x` and `y`.



- int `gsl_blas_scopy`(const gsl_vector_float * *x*, gsl_vector_float * *y*)

- int `gsl_blas_dcopy`(const [gsl_vector](https://www.gnu.org/software/gsl/doc/html/vectors.html#c.gsl_vector) * *x*, [gsl_vector](https://www.gnu.org/software/gsl/doc/html/vectors.html#c.gsl_vector) * *y*)

- int `gsl_blas_ccopy`(const gsl_vector_complex_float * *x*, gsl_vector_complex_float * *y*)

- int `gsl_blas_zcopy`(const gsl_vector_complex * *x*, gsl_vector_complex * *y*)

  These functions copy the elements of the vector `x` into the vector `y`.



- int `gsl_blas_saxpy`(float *alpha*, const gsl_vector_float * *x*, gsl_vector_float * *y*)

- int `gsl_blas_daxpy`(double *alpha*, const [gsl_vector](https://www.gnu.org/software/gsl/doc/html/vectors.html#c.gsl_vector) * *x*, [gsl_vector](https://www.gnu.org/software/gsl/doc/html/vectors.html#c.gsl_vector) * *y*)

- int `gsl_blas_caxpy`(const gsl_complex_float *alpha*, const gsl_vector_complex_float * *x*, gsl_vector_complex_float * *y*)

- int `gsl_blas_zaxpy`(const gsl_complex *alpha*, const gsl_vector_complex * *x*, gsl_vector_complex * *y*)

  These functions compute the sum ![y = \alpha x + y](https://www.gnu.org/software/gsl/doc/html/_images/math/bf4c0bf09422b831ce1c179126efb61a83cc6bb7.png) for the vectors `x` and `y`.



- void `gsl_blas_sscal`(float *alpha*, gsl_vector_float * *x*)

- void `gsl_blas_dscal`(double *alpha*, [gsl_vector](https://www.gnu.org/software/gsl/doc/html/vectors.html#c.gsl_vector) * *x*)

- void `gsl_blas_cscal`(const gsl_complex_float *alpha*, gsl_vector_complex_float * *x*)

- void `gsl_blas_zscal`(const gsl_complex *alpha*, gsl_vector_complex * *x*)

- void `gsl_blas_csscal`(float *alpha*, gsl_vector_complex_float * *x*)

- void `gsl_blas_zdscal`(double *alpha*, gsl_vector_complex * *x*)

  These functions rescale the vector `x` by the multiplicative factor [`alpha`](https://www.gnu.org/software/gsl/doc/html/montecarlo.html#c.alpha).



- int `gsl_blas_srotg`(float *a[]*, float *b[]*, float *c[]*, float *s[]*)

- int `gsl_blas_drotg`(double *a[]*, double *b[]*, double *c[]*, double *s[]*)

  These functions compute a Givens rotation ![(c,s)](https://www.gnu.org/software/gsl/doc/html/_images/math/ebe89dba9e4888c3a244712cf0040eaac9371f4e.png) which zeroes the vector ![(a,b)](https://www.gnu.org/software/gsl/doc/html/_images/math/2dc21aa88e97bfdf1a42cc3b66f9429479b9f536.png),![\left( \begin{matrix}    c & s \\   -s & c \end{matrix} \right) \left( \begin{matrix}   a \\   b \end{matrix} \right) = \left( \begin{matrix}   r' \\   0 \end{matrix} \right)](https://www.gnu.org/software/gsl/doc/html/_images/math/5beb387433ff1670b03ef26c091db7bb6d05f2fc.png)The variables `a` and `b` are overwritten by the routine.

- int `gsl_blas_srot`(gsl_vector_float * *x*, gsl_vector_float * *y*, float *c*, float *s*)

- int `gsl_blas_drot`([gsl_vector](https://www.gnu.org/software/gsl/doc/html/vectors.html#c.gsl_vector) * *x*, [gsl_vector](https://www.gnu.org/software/gsl/doc/html/vectors.html#c.gsl_vector) * *y*, const double *c*, const double *s*)

  These functions apply a Givens rotation ![(x', y') = (c x + s y, -s x + c y)](https://www.gnu.org/software/gsl/doc/html/_images/math/7dbc8f683f3f5c9375cd42d6bd0394e302bdaf30.png) to the vectors `x`, `y`.



- int `gsl_blas_srotmg`(float *d1[]*, float *d2[]*, float *b1[]*, float *b2*, float *P[]*)

- int `gsl_blas_drotmg`(double *d1[]*, double *d2[]*, double *b1[]*, double *b2*, double *P[]*)

  These functions compute a modified Givens transformation. The modified Givens transformation is defined in the original Level-1 BLAS specification, given in the references.

- int `gsl_blas_srotm`(gsl_vector_float * *x*, gsl_vector_float * *y*, const float *P[]*)

- int `gsl_blas_drotm`([gsl_vector](https://www.gnu.org/software/gsl/doc/html/vectors.html#c.gsl_vector) * *x*, [gsl_vector](https://www.gnu.org/software/gsl/doc/html/vectors.html#c.gsl_vector) * *y*, const double *P[]*)

  These functions apply a modified Givens transformation.

### Level 2



- int `gsl_blas_sgemv`(CBLAS_TRANSPOSE_t *TransA*, float *alpha*, const gsl_matrix_float * *A*, const gsl_vector_float * *x*, float *beta*, gsl_vector_float * *y*)

- int `gsl_blas_dgemv`(CBLAS_TRANSPOSE_t *TransA*, double *alpha*, const [gsl_matrix](https://www.gnu.org/software/gsl/doc/html/vectors.html#c.gsl_matrix) * *A*, const [gsl_vector](https://www.gnu.org/software/gsl/doc/html/vectors.html#c.gsl_vector)* *x*, double *beta*, [gsl_vector](https://www.gnu.org/software/gsl/doc/html/vectors.html#c.gsl_vector) * *y*)

- int `gsl_blas_cgemv`(CBLAS_TRANSPOSE_t *TransA*, const gsl_complex_float *alpha*, const gsl_matrix_complex_float * *A*, const gsl_vector_complex_float * *x*, const gsl_complex_float *beta*, gsl_vector_complex_float * *y*)

- int `gsl_blas_zgemv`(CBLAS_TRANSPOSE_t *TransA*, const gsl_complex *alpha*, const gsl_matrix_complex * *A*, const gsl_vector_complex * *x*, const gsl_complex *beta*, gsl_vector_complex * *y*)

  These functions compute the matrix-vector product and sum ![y = \alpha op(A) x + \beta y](https://www.gnu.org/software/gsl/doc/html/_images/math/a5e4a3187a069aa2f7fe0faf5c4524edb0d649e9.png), where ![op(A) = A](https://www.gnu.org/software/gsl/doc/html/_images/math/1483f4ba0ce5c3258bc773466e019a3cfdc04d80.png), ![A^T](https://www.gnu.org/software/gsl/doc/html/_images/math/9fd5ac3afd5fccee5c540ba4086bf7c6cda0cd4f.png), ![A^H](https://www.gnu.org/software/gsl/doc/html/_images/math/9f3adf650c736e0fd66806080be9d8b8d090fdf1.png) for `TransA` = `CblasNoTrans`, `CblasTrans`, `CblasConjTrans`.



- int `gsl_blas_strmv`(CBLAS_UPLO_t *Uplo*, CBLAS_TRANSPOSE_t *TransA*, CBLAS_DIAG_t *Diag*, const gsl_matrix_float * *A*, gsl_vector_float * *x*)

- int `gsl_blas_dtrmv`(CBLAS_UPLO_t *Uplo*, CBLAS_TRANSPOSE_t *TransA*, CBLAS_DIAG_t *Diag*, const [gsl_matrix](https://www.gnu.org/software/gsl/doc/html/vectors.html#c.gsl_matrix) * *A*, [gsl_vector](https://www.gnu.org/software/gsl/doc/html/vectors.html#c.gsl_vector) * *x*)

- int `gsl_blas_ctrmv`(CBLAS_UPLO_t *Uplo*, CBLAS_TRANSPOSE_t *TransA*, CBLAS_DIAG_t *Diag*, const gsl_matrix_complex_float * *A*, gsl_vector_complex_float * *x*)

- int `gsl_blas_ztrmv`(CBLAS_UPLO_t *Uplo*, CBLAS_TRANSPOSE_t *TransA*, CBLAS_DIAG_t *Diag*, const gsl_matrix_complex * *A*, gsl_vector_complex * *x*)

  These functions compute the matrix-vector product ![x = op(A) x](https://www.gnu.org/software/gsl/doc/html/_images/math/a362a8e7c63d56bfc8034c73a09c00070950c159.png) for the triangular matrix `A`, where ![op(A) = A](https://www.gnu.org/software/gsl/doc/html/_images/math/1483f4ba0ce5c3258bc773466e019a3cfdc04d80.png), ![A^T](https://www.gnu.org/software/gsl/doc/html/_images/math/9fd5ac3afd5fccee5c540ba4086bf7c6cda0cd4f.png), ![A^H](https://www.gnu.org/software/gsl/doc/html/_images/math/9f3adf650c736e0fd66806080be9d8b8d090fdf1.png) for `TransA` = `CblasNoTrans`, `CblasTrans`, `CblasConjTrans`. When`Uplo` is `CblasUpper` then the upper triangle of `A` is used, and when `Uplo` is `CblasLower` then the lower triangle of `A` is used. If `Diag` is `CblasNonUnit` then the diagonal of the matrix is used, but if `Diag` is `CblasUnit` then the diagonal elements of the matrix `A` are taken as unity and are not referenced.



- int `gsl_blas_strsv`(CBLAS_UPLO_t *Uplo*, CBLAS_TRANSPOSE_t *TransA*, CBLAS_DIAG_t *Diag*, const gsl_matrix_float * *A*, gsl_vector_float * *x*)

- int `gsl_blas_dtrsv`(CBLAS_UPLO_t *Uplo*, CBLAS_TRANSPOSE_t *TransA*, CBLAS_DIAG_t *Diag*, const [gsl_matrix](https://www.gnu.org/software/gsl/doc/html/vectors.html#c.gsl_matrix) * *A*, [gsl_vector](https://www.gnu.org/software/gsl/doc/html/vectors.html#c.gsl_vector) * *x*)

- int `gsl_blas_ctrsv`(CBLAS_UPLO_t *Uplo*, CBLAS_TRANSPOSE_t *TransA*, CBLAS_DIAG_t *Diag*, const gsl_matrix_complex_float * *A*, gsl_vector_complex_float * *x*)

- int `gsl_blas_ztrsv`(CBLAS_UPLO_t *Uplo*, CBLAS_TRANSPOSE_t *TransA*, CBLAS_DIAG_t *Diag*, const gsl_matrix_complex * *A*, gsl_vector_complex * *x*)

  These functions compute ![inv(op(A)) x](https://www.gnu.org/software/gsl/doc/html/_images/math/39ab3e012448277024df92520e538a9397b8f0f6.png) for `x`, where ![op(A) = A](https://www.gnu.org/software/gsl/doc/html/_images/math/1483f4ba0ce5c3258bc773466e019a3cfdc04d80.png), ![A^T](https://www.gnu.org/software/gsl/doc/html/_images/math/9fd5ac3afd5fccee5c540ba4086bf7c6cda0cd4f.png), ![A^H](https://www.gnu.org/software/gsl/doc/html/_images/math/9f3adf650c736e0fd66806080be9d8b8d090fdf1.png) for `TransA` =`CblasNoTrans`, `CblasTrans`, `CblasConjTrans`. When `Uplo` is `CblasUpper` then the upper triangle of `A` is used, and when `Uplo` is `CblasLower` then the lower triangle of `A` is used. If `Diag` is `CblasNonUnit` then the diagonal of the matrix is used, but if `Diag` is `CblasUnit` then the diagonal elements of the matrix `A` are taken as unity and are not referenced.



- int `gsl_blas_ssymv`(CBLAS_UPLO_t *Uplo*, float *alpha*, const gsl_matrix_float * *A*, const gsl_vector_float * *x*, float *beta*, gsl_vector_float * *y*)

- int `gsl_blas_dsymv`(CBLAS_UPLO_t *Uplo*, double *alpha*, const [gsl_matrix](https://www.gnu.org/software/gsl/doc/html/vectors.html#c.gsl_matrix) * *A*, const [gsl_vector](https://www.gnu.org/software/gsl/doc/html/vectors.html#c.gsl_vector) * *x*, double *beta*, [gsl_vector](https://www.gnu.org/software/gsl/doc/html/vectors.html#c.gsl_vector) * *y*)

  These functions compute the matrix-vector product and sum ![y = \alpha A x + \beta y](https://www.gnu.org/software/gsl/doc/html/_images/math/624871f488f273648511a3e51c9fd2ceb3c4f362.png) for the symmetric matrix `A`. Since the matrix `A` is symmetric only its upper half or lower half need to be stored. When `Uplo` is `CblasUpper` then the upper triangle and diagonal of `A` are used, and when `Uplo` is `CblasLower` then the lower triangle and diagonal of `A` are used.



- int `gsl_blas_chemv`(CBLAS_UPLO_t *Uplo*, const gsl_complex_float *alpha*, const gsl_matrix_complex_float * *A*, const gsl_vector_complex_float * *x*, const gsl_complex_float *beta*, gsl_vector_complex_float * *y*)

- int `gsl_blas_zhemv`(CBLAS_UPLO_t *Uplo*, const gsl_complex *alpha*, const gsl_matrix_complex * *A*, const gsl_vector_complex * *x*, const gsl_complex *beta*, gsl_vector_complex * *y*)

  These functions compute the matrix-vector product and sum ![y = \alpha A x + \beta y](https://www.gnu.org/software/gsl/doc/html/_images/math/624871f488f273648511a3e51c9fd2ceb3c4f362.png) for the hermitian matrix `A`. Since the matrix `A` is hermitian only its upper half or lower half need to be stored. When `Uplo` is `CblasUpper` then the upper triangle and diagonal of `A` are used, and when `Uplo` is `CblasLower` then the lower triangle and diagonal of `A` are used. The imaginary elements of the diagonal are automatically assumed to be zero and are not referenced.



- int `gsl_blas_sger`(float *alpha*, const gsl_vector_float * *x*, const gsl_vector_float * *y*, gsl_matrix_float * *A*)

- int `gsl_blas_dger`(double *alpha*, const [gsl_vector](https://www.gnu.org/software/gsl/doc/html/vectors.html#c.gsl_vector) * *x*, const [gsl_vector](https://www.gnu.org/software/gsl/doc/html/vectors.html#c.gsl_vector) * *y*, [gsl_matrix](https://www.gnu.org/software/gsl/doc/html/vectors.html#c.gsl_matrix) * *A*)

- int `gsl_blas_cgeru`(const gsl_complex_float *alpha*, const gsl_vector_complex_float * *x*, const gsl_vector_complex_float * *y*, gsl_matrix_complex_float * *A*)

- int `gsl_blas_zgeru`(const gsl_complex *alpha*, const gsl_vector_complex * *x*, const gsl_vector_complex * *y*, gsl_matrix_complex * *A*)

  These functions compute the rank-1 update ![A = \alpha x y^T + A](https://www.gnu.org/software/gsl/doc/html/_images/math/e8814216355e71bc23b5d4afb85296d07df4879d.png) of the matrix `A`.



- int `gsl_blas_cgerc`(const gsl_complex_float *alpha*, const gsl_vector_complex_float * *x*, const gsl_vector_complex_float * *y*, gsl_matrix_complex_float * *A*)

- int `gsl_blas_zgerc`(const gsl_complex *alpha*, const gsl_vector_complex * *x*, const gsl_vector_complex * *y*, gsl_matrix_complex * *A*)

  These functions compute the conjugate rank-1 update ![A = \alpha x y^H + A](https://www.gnu.org/software/gsl/doc/html/_images/math/f8611621c1349de247ebb04bf4bc07ed442d4944.png) of the matrix `A`.



- int `gsl_blas_ssyr`(CBLAS_UPLO_t *Uplo*, float *alpha*, const gsl_vector_float * *x*, gsl_matrix_float * *A*)

- int `gsl_blas_dsyr`(CBLAS_UPLO_t *Uplo*, double *alpha*, const [gsl_vector](https://www.gnu.org/software/gsl/doc/html/vectors.html#c.gsl_vector) * *x*, [gsl_matrix](https://www.gnu.org/software/gsl/doc/html/vectors.html#c.gsl_matrix) * *A*)

  These functions compute the symmetric rank-1 update ![A = \alpha x x^T + A](https://www.gnu.org/software/gsl/doc/html/_images/math/78d1d26d1e5af1006cf93134633e4e3f0b4d2c2b.png) of the symmetric matrix `A`. Since the matrix `A` is symmetric only its upper half or lower half need to be stored. When `Uplo` is `CblasUpper` then the upper triangle and diagonal of `A` are used, and when `Uplo` is `CblasLower` then the lower triangle and diagonal of `A` are used.



- int `gsl_blas_cher`(CBLAS_UPLO_t *Uplo*, float *alpha*, const gsl_vector_complex_float * *x*, gsl_matrix_complex_float * *A*)

- int `gsl_blas_zher`(CBLAS_UPLO_t *Uplo*, double *alpha*, const gsl_vector_complex * *x*, gsl_matrix_complex * *A*)

  These functions compute the hermitian rank-1 update ![A = \alpha x x^H + A](https://www.gnu.org/software/gsl/doc/html/_images/math/1b9ebd549dab508e56d26f2a977211ebc72eeafd.png) of the hermitian matrix `A`. Since the matrix `A` is hermitian only its upper half or lower half need to be stored. When`Uplo` is `CblasUpper` then the upper triangle and diagonal of `A` are used, and when `Uplo` is `CblasLower` then the lower triangle and diagonal of `A` are used. The imaginary elements of the diagonal are automatically set to zero.



- int `gsl_blas_ssyr2`(CBLAS_UPLO_t *Uplo*, float *alpha*, const gsl_vector_float * *x*, const gsl_vector_float * *y*, gsl_matrix_float * *A*)

- int `gsl_blas_dsyr2`(CBLAS_UPLO_t *Uplo*, double *alpha*, const [gsl_vector](https://www.gnu.org/software/gsl/doc/html/vectors.html#c.gsl_vector) * *x*, const [gsl_vector](https://www.gnu.org/software/gsl/doc/html/vectors.html#c.gsl_vector) * *y*, [gsl_matrix](https://www.gnu.org/software/gsl/doc/html/vectors.html#c.gsl_matrix) * *A*)

  These functions compute the symmetric rank-2 update ![A = \alpha x y^T + \alpha y x^T + A](https://www.gnu.org/software/gsl/doc/html/_images/math/28a38491faba2f200b769afac1c33efcf1fa1a78.png) of the symmetric matrix `A`. Since the matrix `A` is symmetric only its upper half or lower half need to be stored. When `Uplo` is `CblasUpper` then the upper triangle and diagonal of `A` are used, and when `Uplo` is `CblasLower` then the lower triangle and diagonal of `A` are used.



- int `gsl_blas_cher2`(CBLAS_UPLO_t *Uplo*, const gsl_complex_float *alpha*, const gsl_vector_complex_float * *x*, const gsl_vector_complex_float * *y*, gsl_matrix_complex_float * *A*)

- int `gsl_blas_zher2`(CBLAS_UPLO_t *Uplo*, const gsl_complex *alpha*, const gsl_vector_complex * *x*, const gsl_vector_complex * *y*, gsl_matrix_complex * *A*)

  These functions compute the hermitian rank-2 update ![A = \alpha x y^H + \alpha^* y x^H + A](https://www.gnu.org/software/gsl/doc/html/_images/math/1e688f7ae922888ecd1d564bed4b40bb460fda4c.png) of the hermitian matrix `A`. Since the matrix `A` is hermitian only its upper half or lower half need to be stored. When `Uplo` is `CblasUpper` then the upper triangle and diagonal of `A` are used, and when `Uplo` is `CblasLower` then the lower triangle and diagonal of `A` are used. The imaginary elements of the diagonal are automatically set to zero.

### Level 3



- int `gsl_blas_sgemm`(CBLAS_TRANSPOSE_t *TransA*, CBLAS_TRANSPOSE_t *TransB*, float *alpha*, const gsl_matrix_float * *A*, const gsl_matrix_float * *B*, float *beta*, gsl_matrix_float * *C*)

- int `gsl_blas_dgemm`(CBLAS_TRANSPOSE_t *TransA*, CBLAS_TRANSPOSE_t *TransB*, double *alpha*, const [gsl_matrix](https://www.gnu.org/software/gsl/doc/html/vectors.html#c.gsl_matrix) * *A*, const [gsl_matrix](https://www.gnu.org/software/gsl/doc/html/vectors.html#c.gsl_matrix) * *B*, double *beta*, [gsl_matrix](https://www.gnu.org/software/gsl/doc/html/vectors.html#c.gsl_matrix) * *C*)

- int `gsl_blas_cgemm`(CBLAS_TRANSPOSE_t *TransA*, CBLAS_TRANSPOSE_t *TransB*, const gsl_complex_float *alpha*, const gsl_matrix_complex_float * *A*, const gsl_matrix_complex_float * *B*, const gsl_complex_float *beta*, gsl_matrix_complex_float * *C*)

- int `gsl_blas_zgemm`(CBLAS_TRANSPOSE_t *TransA*, CBLAS_TRANSPOSE_t *TransB*, const gsl_complex *alpha*, const gsl_matrix_complex * *A*, const gsl_matrix_complex * *B*, const gsl_complex *beta*, gsl_matrix_complex * *C*)

  These functions compute the matrix-matrix product and sum ![C = \alpha op(A) op(B) + \beta C](https://www.gnu.org/software/gsl/doc/html/_images/math/f1d617ce30a6240e5dc43578a226eccf14cd8ad8.png) where ![op(A) = A](https://www.gnu.org/software/gsl/doc/html/_images/math/1483f4ba0ce5c3258bc773466e019a3cfdc04d80.png), ![A^T](https://www.gnu.org/software/gsl/doc/html/_images/math/9fd5ac3afd5fccee5c540ba4086bf7c6cda0cd4f.png), ![A^H](https://www.gnu.org/software/gsl/doc/html/_images/math/9f3adf650c736e0fd66806080be9d8b8d090fdf1.png) for `TransA` = `CblasNoTrans`, `CblasTrans`, `CblasConjTrans` and similarly for the parameter `TransB`.



- int `gsl_blas_ssymm`(CBLAS_SIDE_t *Side*, CBLAS_UPLO_t *Uplo*, float *alpha*, const gsl_matrix_float * *A*, const gsl_matrix_float * *B*, float *beta*, gsl_matrix_float * *C*)

- int `gsl_blas_dsymm`(CBLAS_SIDE_t *Side*, CBLAS_UPLO_t *Uplo*, double *alpha*, const [gsl_matrix](https://www.gnu.org/software/gsl/doc/html/vectors.html#c.gsl_matrix) * *A*, const [gsl_matrix](https://www.gnu.org/software/gsl/doc/html/vectors.html#c.gsl_matrix) * *B*, double *beta*, [gsl_matrix](https://www.gnu.org/software/gsl/doc/html/vectors.html#c.gsl_matrix) * *C*)

- int `gsl_blas_csymm`(CBLAS_SIDE_t *Side*, CBLAS_UPLO_t *Uplo*, const gsl_complex_float *alpha*, const gsl_matrix_complex_float * *A*, const gsl_matrix_complex_float * *B*, const gsl_complex_float *beta*, gsl_matrix_complex_float * *C*)

- int `gsl_blas_zsymm`(CBLAS_SIDE_t *Side*, CBLAS_UPLO_t *Uplo*, const gsl_complex *alpha*, const gsl_matrix_complex * *A*, const gsl_matrix_complex * *B*, const gsl_complex *beta*, gsl_matrix_complex * *C*)

  These functions compute the matrix-matrix product and sum ![C = \alpha A B + \beta C](https://www.gnu.org/software/gsl/doc/html/_images/math/e65821bfd5d334138d12f19fbc14055b4808e4a8.png) for `Side` is `CblasLeft` and ![C = \alpha B A + \beta C](https://www.gnu.org/software/gsl/doc/html/_images/math/860474428f0a483215677ed6991e1dd65be68eeb.png) for `Side` is `CblasRight`, where the matrix `A` is symmetric. When `Uplo` is `CblasUpper` then the upper triangle and diagonal of `A` are used, and when `Uplo` is `CblasLower` then the lower triangle and diagonal of `A` are used.



- int `gsl_blas_chemm`(CBLAS_SIDE_t *Side*, CBLAS_UPLO_t *Uplo*, const gsl_complex_float *alpha*, const gsl_matrix_complex_float * *A*, const gsl_matrix_complex_float * *B*, const gsl_complex_float *beta*, gsl_matrix_complex_float * *C*)

- int `gsl_blas_zhemm`(CBLAS_SIDE_t *Side*, CBLAS_UPLO_t *Uplo*, const gsl_complex *alpha*, const gsl_matrix_complex * *A*, const gsl_matrix_complex * *B*, const gsl_complex *beta*, gsl_matrix_complex * *C*)

  These functions compute the matrix-matrix product and sum ![C = \alpha A B + \beta C](https://www.gnu.org/software/gsl/doc/html/_images/math/e65821bfd5d334138d12f19fbc14055b4808e4a8.png) for `Side` is `CblasLeft` and ![C = \alpha B A + \beta C](https://www.gnu.org/software/gsl/doc/html/_images/math/860474428f0a483215677ed6991e1dd65be68eeb.png) for `Side` is `CblasRight`, where the matrix `A` is hermitian. When `Uplo` is `CblasUpper` then the upper triangle and diagonal of `A` are used, and when `Uplo` is `CblasLower` then the lower triangle and diagonal of `A` are used. The imaginary elements of the diagonal are automatically set to zero.



- int `gsl_blas_strmm`(CBLAS_SIDE_t *Side*, CBLAS_UPLO_t *Uplo*, CBLAS_TRANSPOSE_t *TransA*, CBLAS_DIAG_t *Diag*, float *alpha*, const gsl_matrix_float * *A*, gsl_matrix_float * *B*)

- int `gsl_blas_dtrmm`(CBLAS_SIDE_t *Side*, CBLAS_UPLO_t *Uplo*, CBLAS_TRANSPOSE_t *TransA*, CBLAS_DIAG_t *Diag*, double *alpha*, const [gsl_matrix](https://www.gnu.org/software/gsl/doc/html/vectors.html#c.gsl_matrix) * *A*, [gsl_matrix](https://www.gnu.org/software/gsl/doc/html/vectors.html#c.gsl_matrix) * *B*)

- int `gsl_blas_ctrmm`(CBLAS_SIDE_t *Side*, CBLAS_UPLO_t *Uplo*, CBLAS_TRANSPOSE_t *TransA*, CBLAS_DIAG_t *Diag*, const gsl_complex_float *alpha*, const gsl_matrix_complex_float * *A*, gsl_matrix_complex_float * *B*)

- int `gsl_blas_ztrmm`(CBLAS_SIDE_t *Side*, CBLAS_UPLO_t *Uplo*, CBLAS_TRANSPOSE_t *TransA*, CBLAS_DIAG_t *Diag*, const gsl_complex *alpha*, const gsl_matrix_complex * *A*, gsl_matrix_complex * *B*)

  These functions compute the matrix-matrix product ![B = \alpha op(A) B](https://www.gnu.org/software/gsl/doc/html/_images/math/dbfbc5f8691ead41598032122cdfa469586c7885.png) for `Side` is `CblasLeft` and ![B = \alpha B op(A)](https://www.gnu.org/software/gsl/doc/html/_images/math/30f8ff1d58c1e276f8be9f8a173c9e08236b5b45.png) for `Side` is `CblasRight`. The matrix `A` is triangular and ![op(A) = A](https://www.gnu.org/software/gsl/doc/html/_images/math/1483f4ba0ce5c3258bc773466e019a3cfdc04d80.png), ![A^T](https://www.gnu.org/software/gsl/doc/html/_images/math/9fd5ac3afd5fccee5c540ba4086bf7c6cda0cd4f.png), ![A^H](https://www.gnu.org/software/gsl/doc/html/_images/math/9f3adf650c736e0fd66806080be9d8b8d090fdf1.png)for `TransA` = `CblasNoTrans`, `CblasTrans`, `CblasConjTrans`. When `Uplo` is `CblasUpper` then the upper triangle of `A` is used, and when `Uplo` is `CblasLower` then the lower triangle of `A` is used. If `Diag` is `CblasNonUnit` then the diagonal of `A` is used, but if `Diag` is `CblasUnit` then the diagonal elements of the matrix `A` are taken as unity and are not referenced.



- int `gsl_blas_strsm`(CBLAS_SIDE_t *Side*, CBLAS_UPLO_t *Uplo*, CBLAS_TRANSPOSE_t *TransA*, CBLAS_DIAG_t *Diag*, float *alpha*, const gsl_matrix_float * *A*, gsl_matrix_float * *B*)

- int `gsl_blas_dtrsm`(CBLAS_SIDE_t *Side*, CBLAS_UPLO_t *Uplo*, CBLAS_TRANSPOSE_t *TransA*, CBLAS_DIAG_t *Diag*, double *alpha*, const [gsl_matrix](https://www.gnu.org/software/gsl/doc/html/vectors.html#c.gsl_matrix) * *A*, [gsl_matrix](https://www.gnu.org/software/gsl/doc/html/vectors.html#c.gsl_matrix) * *B*)

- int `gsl_blas_ctrsm`(CBLAS_SIDE_t *Side*, CBLAS_UPLO_t *Uplo*, CBLAS_TRANSPOSE_t *TransA*, CBLAS_DIAG_t *Diag*, const gsl_complex_float *alpha*, const gsl_matrix_complex_float * *A*, gsl_matrix_complex_float * *B*)

- int `gsl_blas_ztrsm`(CBLAS_SIDE_t *Side*, CBLAS_UPLO_t *Uplo*, CBLAS_TRANSPOSE_t *TransA*, CBLAS_DIAG_t *Diag*, const gsl_complex *alpha*, const gsl_matrix_complex * *A*, gsl_matrix_complex * *B*)

  These functions compute the inverse-matrix matrix product ![B = \alpha op(inv(A))B](https://www.gnu.org/software/gsl/doc/html/_images/math/c740c4807cd59897acab516d2954e74f10f130b3.png) for `Side` is`CblasLeft` and ![B = \alpha B op(inv(A))](https://www.gnu.org/software/gsl/doc/html/_images/math/c38c56b7f23009164fc706c4df0e6a973ec10671.png) for `Side` is `CblasRight`. The matrix `A` is triangular and ![op(A) = A](https://www.gnu.org/software/gsl/doc/html/_images/math/1483f4ba0ce5c3258bc773466e019a3cfdc04d80.png), ![A^T](https://www.gnu.org/software/gsl/doc/html/_images/math/9fd5ac3afd5fccee5c540ba4086bf7c6cda0cd4f.png), ![A^H](https://www.gnu.org/software/gsl/doc/html/_images/math/9f3adf650c736e0fd66806080be9d8b8d090fdf1.png) for `TransA` = `CblasNoTrans`, `CblasTrans`, `CblasConjTrans`. When `Uplo` is `CblasUpper` then the upper triangle of `A` is used, and when `Uplo` is `CblasLower` then the lower triangle of `A` is used. If `Diag` is `CblasNonUnit` then the diagonal of `A` is used, but if `Diag` is `CblasUnit` then the diagonal elements of the matrix `A` are taken as unity and are not referenced.



- int `gsl_blas_ssyrk`(CBLAS_UPLO_t *Uplo*, CBLAS_TRANSPOSE_t *Trans*, float *alpha*, const gsl_matrix_float * *A*, float *beta*, gsl_matrix_float * *C*)

- int `gsl_blas_dsyrk`(CBLAS_UPLO_t *Uplo*, CBLAS_TRANSPOSE_t *Trans*, double *alpha*, const [gsl_matrix](https://www.gnu.org/software/gsl/doc/html/vectors.html#c.gsl_matrix)* *A*, double *beta*, [gsl_matrix](https://www.gnu.org/software/gsl/doc/html/vectors.html#c.gsl_matrix) * *C*)

- int `gsl_blas_csyrk`(CBLAS_UPLO_t *Uplo*, CBLAS_TRANSPOSE_t *Trans*, const gsl_complex_float *alpha*, const gsl_matrix_complex_float * *A*, const gsl_complex_float *beta*, gsl_matrix_complex_float * *C*)

- int `gsl_blas_zsyrk`(CBLAS_UPLO_t *Uplo*, CBLAS_TRANSPOSE_t *Trans*, const gsl_complex *alpha*, const gsl_matrix_complex * *A*, const gsl_complex *beta*, gsl_matrix_complex * *C*)

  These functions compute a rank-k update of the symmetric matrix `C`, ![C = \alpha A A^T + \beta C](https://www.gnu.org/software/gsl/doc/html/_images/math/77ce32eccbbd4a69fa8b4a98494168aa95ccc6ef.png) when `Trans` is `CblasNoTrans` and ![C = \alpha A^T A + \beta C](https://www.gnu.org/software/gsl/doc/html/_images/math/4f68a8e443706646ec098e8a489b3e28285b4e4d.png) when `Trans` is `CblasTrans`. Since the matrix `C` is symmetric only its upper half or lower half need to be stored. When `Uplo` is `CblasUpper`then the upper triangle and diagonal of `C` are used, and when `Uplo` is `CblasLower` then the lower triangle and diagonal of `C` are used.



- int `gsl_blas_cherk`(CBLAS_UPLO_t *Uplo*, CBLAS_TRANSPOSE_t *Trans*, float *alpha*, const gsl_matrix_complex_float * *A*, float *beta*, gsl_matrix_complex_float * *C*)

- int `gsl_blas_zherk`(CBLAS_UPLO_t *Uplo*, CBLAS_TRANSPOSE_t *Trans*, double *alpha*, const gsl_matrix_complex * *A*, double *beta*, gsl_matrix_complex * *C*)

  These functions compute a rank-k update of the hermitian matrix `C`, ![C = \alpha A A^H + \beta C](https://www.gnu.org/software/gsl/doc/html/_images/math/d6bdfff2f9d1684d281bc62d919b8061da690240.png) when `Trans` is `CblasNoTrans` and ![C = \alpha A^H A + \beta C](https://www.gnu.org/software/gsl/doc/html/_images/math/c0729cd259eb53131087365741284914d45e5b13.png) when `Trans` is `CblasConjTrans`. Since the matrix `C` is hermitian only its upper half or lower half need to be stored. When `Uplo` is`CblasUpper` then the upper triangle and diagonal of `C` are used, and when `Uplo` is `CblasLower`then the lower triangle and diagonal of `C` are used. The imaginary elements of the diagonal are automatically set to zero.



- int `gsl_blas_ssyr2k`(CBLAS_UPLO_t *Uplo*, CBLAS_TRANSPOSE_t *Trans*, float *alpha*, const gsl_matrix_float * *A*, const gsl_matrix_float * *B*, float *beta*, gsl_matrix_float * *C*)

- int `gsl_blas_dsyr2k`(CBLAS_UPLO_t *Uplo*, CBLAS_TRANSPOSE_t *Trans*, double *alpha*, const [gsl_matrix](https://www.gnu.org/software/gsl/doc/html/vectors.html#c.gsl_matrix)* *A*, const [gsl_matrix](https://www.gnu.org/software/gsl/doc/html/vectors.html#c.gsl_matrix) * *B*, double *beta*, [gsl_matrix](https://www.gnu.org/software/gsl/doc/html/vectors.html#c.gsl_matrix) * *C*)

- int `gsl_blas_csyr2k`(CBLAS_UPLO_t *Uplo*, CBLAS_TRANSPOSE_t *Trans*, const gsl_complex_float *alpha*, const gsl_matrix_complex_float * *A*, const gsl_matrix_complex_float * *B*, const gsl_complex_float *beta*, gsl_matrix_complex_float * *C*)

- int `gsl_blas_zsyr2k`(CBLAS_UPLO_t *Uplo*, CBLAS_TRANSPOSE_t *Trans*, const gsl_complex *alpha*, const gsl_matrix_complex * *A*, const gsl_matrix_complex * *B*, const gsl_complex *beta*, gsl_matrix_complex * *C*)

  These functions compute a rank-2k update of the symmetric matrix `C`, ![C = \alpha A B^T + \alpha B A^T + \beta C](https://www.gnu.org/software/gsl/doc/html/_images/math/053ba1f9432b036f0e14ce1e8b2fc51c5d329ef3.png) when `Trans` is `CblasNoTrans` and ![C = \alpha A^T B + \alpha B^T A + \beta C](https://www.gnu.org/software/gsl/doc/html/_images/math/e4e735aad28bc1139b86415a500b78cc8b0eee1c.png)when `Trans` is `CblasTrans`. Since the matrix `C` is symmetric only its upper half or lower half need to be stored. When `Uplo` is `CblasUpper` then the upper triangle and diagonal of `C` are used, and when `Uplo` is `CblasLower` then the lower triangle and diagonal of `C` are used.



- int `gsl_blas_cher2k`(CBLAS_UPLO_t *Uplo*, CBLAS_TRANSPOSE_t *Trans*, const gsl_complex_float *alpha*, const gsl_matrix_complex_float * *A*, const gsl_matrix_complex_float * *B*, float *beta*, gsl_matrix_complex_float * *C*)

- int `gsl_blas_zher2k`(CBLAS_UPLO_t *Uplo*, CBLAS_TRANSPOSE_t *Trans*, const gsl_complex *alpha*, const gsl_matrix_complex * *A*, const gsl_matrix_complex * *B*, double *beta*, gsl_matrix_complex * *C*)

  These functions compute a rank-2k update of the hermitian matrix `C`, ![C = \alpha A B^H + \alpha^* B A^H + \beta C](https://www.gnu.org/software/gsl/doc/html/_images/math/beaae7206d62aa5cb21b09ed2c0e1c891dc5a9fc.png) when `Trans` is `CblasNoTrans` and ![C = \alpha A^H B + \alpha^* B^H A + \beta C](https://www.gnu.org/software/gsl/doc/html/_images/math/a045c75a24674e0ae5c1bc22c287d7d64288b5ad.png) when `Trans` is `CblasConjTrans`. Since the matrix `C` is hermitian only its upper half or lower half need to be stored. When `Uplo` is `CblasUpper` then the upper triangle and diagonal of `C` are used, and when `Uplo` is `CblasLower` then the lower triangle and diagonal of `C` are used. The imaginary elements of the diagonal are automatically set to zero.

## Examples

The following program computes the product of two matrices using the Level-3 BLAS function DGEMM,

![\left( \begin{matrix}   0.11&0.12&0.13 \\   0.21&0.22&0.23 \end{matrix} \right) \left( \begin{matrix}   1011&1012 \\   1021&1022 \\   1031&1031 \end{matrix} \right) = \left( \begin{matrix}   367.76&368.12 \\   674.06&674.72 \end{matrix} \right)](https://www.gnu.org/software/gsl/doc/html/_images/math/295744c265e2375eea35cb5df906a7b8100f81a5.png)

The matrices are stored in row major order, according to the C convention for arrays.

```
#include <stdio.h>
#include <gsl/gsl_blas.h>

int
main (void)
{
  double a[] = { 0.11, 0.12, 0.13,
                 0.21, 0.22, 0.23 };

  double b[] = { 1011, 1012,
                 1021, 1022,
                 1031, 1032 };

  double c[] = { 0.00, 0.00,
                 0.00, 0.00 };

  gsl_matrix_view A = gsl_matrix_view_array(a, 2, 3);
  gsl_matrix_view B = gsl_matrix_view_array(b, 3, 2);
  gsl_matrix_view C = gsl_matrix_view_array(c, 2, 2);

  /* Compute C = A B */

  gsl_blas_dgemm (CblasNoTrans, CblasNoTrans,
                  1.0, &A.matrix, &B.matrix,
                  0.0, &C.matrix);

  printf ("[ %g, %g\n", c[0], c[1]);
  printf ("  %g, %g ]\n", c[2], c[3]);

  return 0;
}
```

Here is the output from the program,

```
[ 367.76, 368.12
  674.06, 674.72 ]
```



## References and Further Reading

Information on the BLAS standards, including both the legacy and updated interface standards, is available online from the BLAS Homepage and BLAS Technical Forum web-site.

- BLAS Homepage, <http://www.netlib.org/blas/>
- BLAS Technical Forum, <http://www.netlib.org/blas/blast-forum/>

The following papers contain the specifications for Level 1, Level 2 and Level 3 BLAS.

- C. Lawson, R. Hanson, D. Kincaid, F. Krogh, “Basic Linear Algebra Subprograms for Fortran Usage”, ACM Transactions on Mathematical Software, Vol.: 5 (1979), Pages 308–325.
- J.J. Dongarra, J. DuCroz, S. Hammarling, R. Hanson, “An Extended Set of Fortran Basic Linear Algebra Subprograms”, ACM Transactions on Mathematical Software, Vol.: 14, No.: 1 (1988), Pages 1–32.
- J.J. Dongarra, I. Duff, J. DuCroz, S. Hammarling, “A Set of Level 3 Basic Linear Algebra Subprograms”, ACM Transactions on Mathematical Software, Vol.: 16 (1990), Pages 1–28.

Postscript versions of the latter two papers are available from <http://www.netlib.org/blas/>. A CBLAS wrapper for Fortran BLAS libraries is available from the same location.

### Footnotes

[[1]](https://www.gnu.org/software/gsl/doc/html/blas.html#id1)In the low-level CBLAS interface, a negative stride accesses the vector elements in reverse order, i.e. the ![i](https://www.gnu.org/software/gsl/doc/html/_images/math/1cc632900aae8b2837a1f383619c0ad753be7d29.png)-th element is given by ![(N-i)*|incx|](https://www.gnu.org/software/gsl/doc/html/_images/math/06be90af20cf3efc17d4b4cdd99f56502f4faa62.png) for ![incx < 0](https://www.gnu.org/software/gsl/doc/html/_images/math/db34591c31985b279a480df200850f9965de7a15.png).