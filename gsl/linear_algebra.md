# Linear Algebra 线性代数
This chapter describes functions for solving linear systems. The library provides linear algebra operations which operate directly on the [`gsl_vector`](https://www.gnu.org/software/gsl/doc/html/vectors.html#c.gsl_vector) and [`gsl_matrix`](https://www.gnu.org/software/gsl/doc/html/vectors.html#c.gsl_matrix) objects. These routines use the standard algorithms from Golub & Van Loan’s *Matrix Computations* with Level-1 and Level-2 BLAS calls for efficiency.
本章描述用于求解线性系统的函数。本库提供了线性代数操作，其直接作用在[`gsl_vector`](https://www.gnu.org/softwares/gsl/doc/html/vectors.html#c.gsl_vector)和[`gsl_matrix`](https://www.gnu.org/software/gsl/doc/html/vectors.html#c.gsl_matrix)对象上。这些程序使用来自Golub和Van Loan的*Matrix Computation* 的标准算法，通过1级和2级BLAS调用来寻求效率。

The functions described in this chapter are declared in the header file `gsl_linalg.h`.
本章描述的函数被声明在头文件`gsl_linalg.h`中。


## LU Decomposition LU分解

A general ![N](https://www.gnu.org/software/gsl/doc/html/_images/math/99f8a3c97c19212ca3eee7d00dddc8d68b2d1c27.png)-by-![N](https://www.gnu.org/software/gsl/doc/html/_images/math/99f8a3c97c19212ca3eee7d00dddc8d68b2d1c27.png) square matrix ![A](https://www.gnu.org/software/gsl/doc/html/_images/math/3daf97c82cdab413e28d010c06770f224d64dc5a.png) has an ![LU](https://www.gnu.org/software/gsl/doc/html/_images/math/dba5651ef03e528e7fd6d4ae7afda2eaa9836116.png) decomposition into upper and lower triangular matrices,

![P A = L U](https://www.gnu.org/software/gsl/doc/html/_images/math/a7abeeef102edfbffeb4f88fa5eca9805ff5bd14.png)

一个通用的$N-by-N$方阵$A$具有一个$$LU$$分解为一个上三角阵和一个下三角阵，

$PA=LU$

where ![P](https://www.gnu.org/software/gsl/doc/html/_images/math/d1a0a8c0064fe114c19c8fd71ffec7a98480f707.png) is a permutation matrix, ![L](https://www.gnu.org/software/gsl/doc/html/_images/math/82359a7c56cc02d70a5ed0e98df631577f7d172b.png) is unit lower triangular matrix and ![U](https://www.gnu.org/software/gsl/doc/html/_images/math/3d0e8e66abfde1273e99e86e8a090a7073d842e3.png) is upper triangular matrix. For square matrices this decomposition can be used to convert the linear system ![A x = b](https://www.gnu.org/software/gsl/doc/html/_images/math/bf1b0fe2f3273671dbd86939aa2ca8992c570e65.png) into a pair of triangular systems (![L y = P b](https://www.gnu.org/software/gsl/doc/html/_images/math/bf9541228f059a57b6efe8bafb08e1cedad98a58.png), ![U x = y](https://www.gnu.org/software/gsl/doc/html/_images/math/a349be500c98b2fe20c54c2f8c0d9f9d78b3b980.png)), which can be solved by forward and back-substitution. Note that the ![LU](https://www.gnu.org/software/gsl/doc/html/_images/math/dba5651ef03e528e7fd6d4ae7afda2eaa9836116.png) decomposition is valid for singular matrices.

其中$$P$$是一个置换矩阵，$$L$$是一个单位下三角矩阵而$$U$$是一个上三角矩阵。对于方阵，这个分解能够被用于转换线性系统$$Ax=b$$为一对三角系统$$(Ly=Pb, Ux=y)$$，其能够由前向和后向代换法求解。注意$$LU$$分解仅对奇异阵有效。

- int `gsl_linalg_LU_decomp`([gsl_matrix](https://www.gnu.org/software/gsl/doc/html/vectors.html#c.gsl_matrix) * *A*, [gsl_permutation](https://www.gnu.org/software/gsl/doc/html/permutation.html#c.gsl_permutation) * *p*, int * *signum*)

- int `gsl_linalg_complex_LU_decomp`(gsl_matrix_complex * *A*, [gsl_permutation](https://www.gnu.org/software/gsl/doc/html/permutation.html#c.gsl_permutation) * *p*, int * *signum*)

  These functions factorize the square matrix `A` into the ![LU](https://www.gnu.org/software/gsl/doc/html/_images/math/dba5651ef03e528e7fd6d4ae7afda2eaa9836116.png) decomposition ![PA = LU](https://www.gnu.org/software/gsl/doc/html/_images/math/d3e1ef726a15da8060f714833cf273f27b1e0ca4.png). On output the diagonal and upper triangular part of the input matrix `A` contain the matrix ![U](https://www.gnu.org/software/gsl/doc/html/_images/math/3d0e8e66abfde1273e99e86e8a090a7073d842e3.png). The lower triangular part of the input matrix (excluding the diagonal) contains ![L](https://www.gnu.org/software/gsl/doc/html/_images/math/82359a7c56cc02d70a5ed0e98df631577f7d172b.png). The diagonal elements of ![L](https://www.gnu.org/software/gsl/doc/html/_images/math/82359a7c56cc02d70a5ed0e98df631577f7d172b.png) are unity, and are not stored.The permutation matrix ![P](https://www.gnu.org/software/gsl/doc/html/_images/math/d1a0a8c0064fe114c19c8fd71ffec7a98480f707.png) is encoded in the permutation `p` on output. The ![j](https://www.gnu.org/software/gsl/doc/html/_images/math/7bc2ad053edeecad39b867d60495ac2caa4555dc.png)-th column of the matrix ![P](https://www.gnu.org/software/gsl/doc/html/_images/math/d1a0a8c0064fe114c19c8fd71ffec7a98480f707.png) is given by the ![k](https://www.gnu.org/software/gsl/doc/html/_images/math/90645b5c2abc13d7de104fd87abeb2f19406e5f8.png)-th column of the identity matrix, where ![k = p_j](https://www.gnu.org/software/gsl/doc/html/_images/math/180cd8182731bf565d61a5a2d33ede4726b8e5e9.png) the ![j](https://www.gnu.org/software/gsl/doc/html/_images/math/7bc2ad053edeecad39b867d60495ac2caa4555dc.png)-th element of the permutation vector. The sign of the permutation is given by `signum`. It has the value ![(-1)^n](https://www.gnu.org/software/gsl/doc/html/_images/math/11f2b5e69df0c211c2b82609c97389c1f9d8df14.png), where ![n](https://www.gnu.org/software/gsl/doc/html/_images/math/a24554f1502cf204e7aa24a6c064962c3504de48.png) is the number of interchanges in the permutation.The algorithm used in the decomposition is Gaussian Elimination with partial pivoting (Golub & Van Loan, *Matrix Computations*, Algorithm 3.4.1).

  这些函数

- int `gsl_linalg_LU_solve`(const [gsl_matrix](https://www.gnu.org/software/gsl/doc/html/vectors.html#c.gsl_matrix) * *LU*, const [gsl_permutation](https://www.gnu.org/software/gsl/doc/html/permutation.html#c.gsl_permutation) * *p*, const [gsl_vector](https://www.gnu.org/software/gsl/doc/html/vectors.html#c.gsl_vector) * *b*, [gsl_vector](https://www.gnu.org/software/gsl/doc/html/vectors.html#c.gsl_vector) * *x*)

- int `gsl_linalg_complex_LU_solve`(const gsl_matrix_complex * *LU*, const [gsl_permutation](https://www.gnu.org/software/gsl/doc/html/permutation.html#c.gsl_permutation) * *p*, const gsl_vector_complex * *b*, gsl_vector_complex * *x*)

  These functions solve the square system ![A x = b](https://www.gnu.org/software/gsl/doc/html/_images/math/bf1b0fe2f3273671dbd86939aa2ca8992c570e65.png) using the ![LU](https://www.gnu.org/software/gsl/doc/html/_images/math/dba5651ef03e528e7fd6d4ae7afda2eaa9836116.png) decomposition of ![A](https://www.gnu.org/software/gsl/doc/html/_images/math/3daf97c82cdab413e28d010c06770f224d64dc5a.png) into (`LU`, `p`) given by [`gsl_linalg_LU_decomp()`](https://www.gnu.org/software/gsl/doc/html/linalg.html#c.gsl_linalg_LU_decomp) or [`gsl_linalg_complex_LU_decomp()`](https://www.gnu.org/software/gsl/doc/html/linalg.html#c.gsl_linalg_complex_LU_decomp) as input.

- int `gsl_linalg_LU_svx`(const [gsl_matrix](https://www.gnu.org/software/gsl/doc/html/vectors.html#c.gsl_matrix) * *LU*, const [gsl_permutation](https://www.gnu.org/software/gsl/doc/html/permutation.html#c.gsl_permutation) * *p*, [gsl_vector](https://www.gnu.org/software/gsl/doc/html/vectors.html#c.gsl_vector) * *x*)

- int `gsl_linalg_complex_LU_svx`(const gsl_matrix_complex * *LU*, const [gsl_permutation](https://www.gnu.org/software/gsl/doc/html/permutation.html#c.gsl_permutation) * *p*, gsl_vector_complex * *x*)

  These functions solve the square system ![A x = b](https://www.gnu.org/software/gsl/doc/html/_images/math/bf1b0fe2f3273671dbd86939aa2ca8992c570e65.png) in-place using the precomputed ![LU](https://www.gnu.org/software/gsl/doc/html/_images/math/dba5651ef03e528e7fd6d4ae7afda2eaa9836116.png)decomposition of ![A](https://www.gnu.org/software/gsl/doc/html/_images/math/3daf97c82cdab413e28d010c06770f224d64dc5a.png) into (`LU`, `p`). On input `x` should contain the right-hand side ![b](https://www.gnu.org/software/gsl/doc/html/_images/math/88082f6488798774b9cd95a59c1ab4c6cafe909f.png), which is replaced by the solution on output.



- int `gsl_linalg_LU_refine`(const [gsl_matrix](https://www.gnu.org/software/gsl/doc/html/vectors.html#c.gsl_matrix) * *A*, const [gsl_matrix](https://www.gnu.org/software/gsl/doc/html/vectors.html#c.gsl_matrix) * *LU*, const [gsl_permutation](https://www.gnu.org/software/gsl/doc/html/permutation.html#c.gsl_permutation) * *p*, const [gsl_vector](https://www.gnu.org/software/gsl/doc/html/vectors.html#c.gsl_vector) * *b*, [gsl_vector](https://www.gnu.org/software/gsl/doc/html/vectors.html#c.gsl_vector) * *x*, [gsl_vector](https://www.gnu.org/software/gsl/doc/html/vectors.html#c.gsl_vector) * *work*)

- int `gsl_linalg_complex_LU_refine`(const gsl_matrix_complex * *A*, const gsl_matrix_complex * *LU*, const [gsl_permutation](https://www.gnu.org/software/gsl/doc/html/permutation.html#c.gsl_permutation) * *p*, const gsl_vector_complex * *b*, gsl_vector_complex * *x*, gsl_vector_complex * *work*)

  These functions apply an iterative improvement to `x`, the solution of ![A x = b](https://www.gnu.org/software/gsl/doc/html/_images/math/bf1b0fe2f3273671dbd86939aa2ca8992c570e65.png), from the precomputed ![LU](https://www.gnu.org/software/gsl/doc/html/_images/math/dba5651ef03e528e7fd6d4ae7afda2eaa9836116.png) decomposition of ![A](https://www.gnu.org/software/gsl/doc/html/_images/math/3daf97c82cdab413e28d010c06770f224d64dc5a.png) into (`LU`, `p`). Additional workspace of length `N` is required in `work`.



- int `gsl_linalg_LU_invert`(const [gsl_matrix](https://www.gnu.org/software/gsl/doc/html/vectors.html#c.gsl_matrix) * *LU*, const [gsl_permutation](https://www.gnu.org/software/gsl/doc/html/permutation.html#c.gsl_permutation) * *p*, [gsl_matrix](https://www.gnu.org/software/gsl/doc/html/vectors.html#c.gsl_matrix) * *inverse*)

- int `gsl_linalg_complex_LU_invert`(const gsl_matrix_complex * *LU*, const [gsl_permutation](https://www.gnu.org/software/gsl/doc/html/permutation.html#c.gsl_permutation) * *p*, gsl_matrix_complex * *inverse*)

  These functions compute the inverse of a matrix ![A](https://www.gnu.org/software/gsl/doc/html/_images/math/3daf97c82cdab413e28d010c06770f224d64dc5a.png) from its ![LU](https://www.gnu.org/software/gsl/doc/html/_images/math/dba5651ef03e528e7fd6d4ae7afda2eaa9836116.png) decomposition (`LU`, `p`), storing the result in the matrix `inverse`. The inverse is computed by solving the system ![A x = b](https://www.gnu.org/software/gsl/doc/html/_images/math/bf1b0fe2f3273671dbd86939aa2ca8992c570e65.png) for each column of the identity matrix. It is preferable to avoid direct use of the inverse whenever possible, as the linear solver functions can obtain the same result more efficiently and reliably (consult any introductory textbook on numerical linear algebra for details).



- double `gsl_linalg_LU_det`([gsl_matrix](https://www.gnu.org/software/gsl/doc/html/vectors.html#c.gsl_matrix) * *LU*, int *signum*)

- gsl_complex `gsl_linalg_complex_LU_det`(gsl_matrix_complex * *LU*, int *signum*)

  These functions compute the determinant of a matrix ![A](https://www.gnu.org/software/gsl/doc/html/_images/math/3daf97c82cdab413e28d010c06770f224d64dc5a.png) from its ![LU](https://www.gnu.org/software/gsl/doc/html/_images/math/dba5651ef03e528e7fd6d4ae7afda2eaa9836116.png) decomposition, `LU`. The determinant is computed as the product of the diagonal elements of ![U](https://www.gnu.org/software/gsl/doc/html/_images/math/3d0e8e66abfde1273e99e86e8a090a7073d842e3.png) and the sign of the row permutation `signum`.



- double `gsl_linalg_LU_lndet`([gsl_matrix](https://www.gnu.org/software/gsl/doc/html/vectors.html#c.gsl_matrix) * *LU*)

- double `gsl_linalg_complex_LU_lndet`(gsl_matrix_complex * *LU*)

  These functions compute the logarithm of the absolute value of the determinant of a matrix ![A](https://www.gnu.org/software/gsl/doc/html/_images/math/3daf97c82cdab413e28d010c06770f224d64dc5a.png), ![\ln|\det(A)|](https://www.gnu.org/software/gsl/doc/html/_images/math/8d50e5648f292131c671961bd0abb57bd668ba48.png), from its ![LU](https://www.gnu.org/software/gsl/doc/html/_images/math/dba5651ef03e528e7fd6d4ae7afda2eaa9836116.png) decomposition, `LU`. This function may be useful if the direct computation of the determinant would overflow or underflow.



- int `gsl_linalg_LU_sgndet`([gsl_matrix](https://www.gnu.org/software/gsl/doc/html/vectors.html#c.gsl_matrix) * *LU*, int *signum*)

- gsl_complex `gsl_linalg_complex_LU_sgndet`(gsl_matrix_complex * *LU*, int *signum*)

  These functions compute the sign or phase factor of the determinant of a matrix ![A](https://www.gnu.org/software/gsl/doc/html/_images/math/3daf97c82cdab413e28d010c06770f224d64dc5a.png), ![\det(A)/|\det(A)|](https://www.gnu.org/software/gsl/doc/html/_images/math/0fb15790b84a146ca5d077362cb952bc6b755656.png), from its ![LU](https://www.gnu.org/software/gsl/doc/html/_images/math/dba5651ef03e528e7fd6d4ae7afda2eaa9836116.png) decomposition, `LU`.



## QR Decomposition

A general rectangular ![M](https://www.gnu.org/software/gsl/doc/html/_images/math/85ca40542b59eef42eec1488ce2784e307962a5f.png)-by-![N](https://www.gnu.org/software/gsl/doc/html/_images/math/99f8a3c97c19212ca3eee7d00dddc8d68b2d1c27.png) matrix ![A](https://www.gnu.org/software/gsl/doc/html/_images/math/3daf97c82cdab413e28d010c06770f224d64dc5a.png) has a ![QR](https://www.gnu.org/software/gsl/doc/html/_images/math/467c97789b8b4e6d169789aa163ba055c58e9223.png) decomposition into the product of an orthogonal ![M](https://www.gnu.org/software/gsl/doc/html/_images/math/85ca40542b59eef42eec1488ce2784e307962a5f.png)-by-![M](https://www.gnu.org/software/gsl/doc/html/_images/math/85ca40542b59eef42eec1488ce2784e307962a5f.png) square matrix ![Q](https://www.gnu.org/software/gsl/doc/html/_images/math/7fcf2b077fb63ef0b6eb5c350673266b94327847.png) (where ![Q^T Q = I](https://www.gnu.org/software/gsl/doc/html/_images/math/6eecea648f2f23c81d354f00be66df644ccdfaec.png)) and an ![M](https://www.gnu.org/software/gsl/doc/html/_images/math/85ca40542b59eef42eec1488ce2784e307962a5f.png)-by-![N](https://www.gnu.org/software/gsl/doc/html/_images/math/99f8a3c97c19212ca3eee7d00dddc8d68b2d1c27.png) right-triangular matrix ![R](https://www.gnu.org/software/gsl/doc/html/_images/math/c9b67a086a723a5eec85a146053dd0a2d2378bac.png),

![A = Q R](https://www.gnu.org/software/gsl/doc/html/_images/math/d676906375e55f0faf904af4652e652c743efa02.png)

This decomposition can be used to convert the linear system ![A x = b](https://www.gnu.org/software/gsl/doc/html/_images/math/bf1b0fe2f3273671dbd86939aa2ca8992c570e65.png) into the triangular system ![R x = Q^T b](https://www.gnu.org/software/gsl/doc/html/_images/math/b3fd8e1be426224ac9aa267b39a9de422b0f4c53.png), which can be solved by back-substitution. Another use of the ![QR](https://www.gnu.org/software/gsl/doc/html/_images/math/467c97789b8b4e6d169789aa163ba055c58e9223.png) decomposition is to compute an orthonormal basis for a set of vectors. The first ![N](https://www.gnu.org/software/gsl/doc/html/_images/math/99f8a3c97c19212ca3eee7d00dddc8d68b2d1c27.png) columns of ![Q](https://www.gnu.org/software/gsl/doc/html/_images/math/7fcf2b077fb63ef0b6eb5c350673266b94327847.png) form an orthonormal basis for the range of ![A](https://www.gnu.org/software/gsl/doc/html/_images/math/3daf97c82cdab413e28d010c06770f224d64dc5a.png), ![ran(A)](https://www.gnu.org/software/gsl/doc/html/_images/math/76d3fd64ee04f68953955d1372dfbcd8e0d38a18.png), when ![A](https://www.gnu.org/software/gsl/doc/html/_images/math/3daf97c82cdab413e28d010c06770f224d64dc5a.png) has full column rank.

- int `gsl_linalg_QR_decomp`([gsl_matrix](https://www.gnu.org/software/gsl/doc/html/vectors.html#c.gsl_matrix) * *A*, [gsl_vector](https://www.gnu.org/software/gsl/doc/html/vectors.html#c.gsl_vector) * *tau*)

  This function factorizes the ![M](https://www.gnu.org/software/gsl/doc/html/_images/math/85ca40542b59eef42eec1488ce2784e307962a5f.png)-by-![N](https://www.gnu.org/software/gsl/doc/html/_images/math/99f8a3c97c19212ca3eee7d00dddc8d68b2d1c27.png) matrix `A` into the ![QR](https://www.gnu.org/software/gsl/doc/html/_images/math/467c97789b8b4e6d169789aa163ba055c58e9223.png) decomposition ![A = Q R](https://www.gnu.org/software/gsl/doc/html/_images/math/b402ff0e8fbca3329adbdf5e61d0d3e030be41cb.png). On output the diagonal and upper triangular part of the input matrix contain the matrix ![R](https://www.gnu.org/software/gsl/doc/html/_images/math/c9b67a086a723a5eec85a146053dd0a2d2378bac.png). The vector `tau` and the columns of the lower triangular part of the matrix `A` contain the Householder coefficients and Householder vectors which encode the orthogonal matrix `Q`. The vector `tau` must be of length ![k=\min(M,N)](https://www.gnu.org/software/gsl/doc/html/_images/math/cec598e7cfbb40428221e784f7d5bacd36bc8046.png). The matrix ![Q](https://www.gnu.org/software/gsl/doc/html/_images/math/7fcf2b077fb63ef0b6eb5c350673266b94327847.png) is related to these components by, ![Q = Q_k ... Q_2 Q_1](https://www.gnu.org/software/gsl/doc/html/_images/math/394c1fcb62e3c6b327ff0d104bb5917fa5f19b84.png) where ![Q_i = I - \tau_i v_i v_i^T](https://www.gnu.org/software/gsl/doc/html/_images/math/724f9e9a84326f1f8d512c16f51ec3ba5ff89b63.png) and ![v_i](https://www.gnu.org/software/gsl/doc/html/_images/math/78498e145a0277e0022f7d72bb541ef552325f62.png) is the Householder vector ![v_i = (0,...,1,A(i+1,i),A(i+2,i),...,A(m,i))](https://www.gnu.org/software/gsl/doc/html/_images/math/878186d588411d9510a013b62cb4ef2b6c55a11e.png). This is the same storage scheme as used by LAPACK.The algorithm used to perform the decomposition is Householder QR (Golub & Van Loan, “Matrix Computations”, Algorithm 5.2.1).

- int `gsl_linalg_QR_solve`(const [gsl_matrix](https://www.gnu.org/software/gsl/doc/html/vectors.html#c.gsl_matrix) * *QR*, const [gsl_vector](https://www.gnu.org/software/gsl/doc/html/vectors.html#c.gsl_vector) * *tau*, const [gsl_vector](https://www.gnu.org/software/gsl/doc/html/vectors.html#c.gsl_vector) * *b*, [gsl_vector](https://www.gnu.org/software/gsl/doc/html/vectors.html#c.gsl_vector) * *x*)

  This function solves the square system ![A x = b](https://www.gnu.org/software/gsl/doc/html/_images/math/bf1b0fe2f3273671dbd86939aa2ca8992c570e65.png) using the ![QR](https://www.gnu.org/software/gsl/doc/html/_images/math/467c97789b8b4e6d169789aa163ba055c58e9223.png) decomposition of ![A](https://www.gnu.org/software/gsl/doc/html/_images/math/3daf97c82cdab413e28d010c06770f224d64dc5a.png) held in (`QR`, `tau`) which must have been computed previously with [`gsl_linalg_QR_decomp()`](https://www.gnu.org/software/gsl/doc/html/linalg.html#c.gsl_linalg_QR_decomp). The least-squares solution for rectangular systems can be found using [`gsl_linalg_QR_lssolve()`](https://www.gnu.org/software/gsl/doc/html/linalg.html#c.gsl_linalg_QR_lssolve).

- int `gsl_linalg_QR_svx`(const [gsl_matrix](https://www.gnu.org/software/gsl/doc/html/vectors.html#c.gsl_matrix) * *QR*, const [gsl_vector](https://www.gnu.org/software/gsl/doc/html/vectors.html#c.gsl_vector) * *tau*, [gsl_vector](https://www.gnu.org/software/gsl/doc/html/vectors.html#c.gsl_vector) * *x*)

  This function solves the square system ![A x = b](https://www.gnu.org/software/gsl/doc/html/_images/math/bf1b0fe2f3273671dbd86939aa2ca8992c570e65.png) in-place using the ![QR](https://www.gnu.org/software/gsl/doc/html/_images/math/467c97789b8b4e6d169789aa163ba055c58e9223.png) decomposition of ![A](https://www.gnu.org/software/gsl/doc/html/_images/math/3daf97c82cdab413e28d010c06770f224d64dc5a.png)held in (`QR`, `tau`) which must have been computed previously by [`gsl_linalg_QR_decomp()`](https://www.gnu.org/software/gsl/doc/html/linalg.html#c.gsl_linalg_QR_decomp). On input `x` should contain the right-hand side ![b](https://www.gnu.org/software/gsl/doc/html/_images/math/88082f6488798774b9cd95a59c1ab4c6cafe909f.png), which is replaced by the solution on output.

- int `gsl_linalg_QR_lssolve`(const [gsl_matrix](https://www.gnu.org/software/gsl/doc/html/vectors.html#c.gsl_matrix) * *QR*, const [gsl_vector](https://www.gnu.org/software/gsl/doc/html/vectors.html#c.gsl_vector) * *tau*, const [gsl_vector](https://www.gnu.org/software/gsl/doc/html/vectors.html#c.gsl_vector) * *b*, [gsl_vector](https://www.gnu.org/software/gsl/doc/html/vectors.html#c.gsl_vector) * *x*, [gsl_vector](https://www.gnu.org/software/gsl/doc/html/vectors.html#c.gsl_vector) * *residual*)

  This function finds the least squares solution to the overdetermined system ![A x = b](https://www.gnu.org/software/gsl/doc/html/_images/math/bf1b0fe2f3273671dbd86939aa2ca8992c570e65.png) where the matrix `A` has more rows than columns. The least squares solution minimizes the Euclidean norm of the residual, ![||Ax - b||](https://www.gnu.org/software/gsl/doc/html/_images/math/9c0a20cca26cd490f519bfa5af388bb9f7a7e877.png).The routine requires as input the ![QR](https://www.gnu.org/software/gsl/doc/html/_images/math/467c97789b8b4e6d169789aa163ba055c58e9223.png) decomposition of ![A](https://www.gnu.org/software/gsl/doc/html/_images/math/3daf97c82cdab413e28d010c06770f224d64dc5a.png) into (`QR`, `tau`) given by [`gsl_linalg_QR_decomp()`](https://www.gnu.org/software/gsl/doc/html/linalg.html#c.gsl_linalg_QR_decomp). The solution is returned in `x`. The residual is computed as a by-product and stored in `residual`.

- int `gsl_linalg_QR_QTvec`(const [gsl_matrix](https://www.gnu.org/software/gsl/doc/html/vectors.html#c.gsl_matrix) * *QR*, const [gsl_vector](https://www.gnu.org/software/gsl/doc/html/vectors.html#c.gsl_vector) * *tau*, [gsl_vector](https://www.gnu.org/software/gsl/doc/html/vectors.html#c.gsl_vector) * *v*)

  This function applies the matrix ![Q^T](https://www.gnu.org/software/gsl/doc/html/_images/math/7f8ec85a5eb1167b604a5766a3c0dc526d37caaf.png) encoded in the decomposition (`QR`, `tau`) to the vector `v`, storing the result ![Q^T v](https://www.gnu.org/software/gsl/doc/html/_images/math/556b080d1d11975b54e4cc840b491f2d31e5d6c4.png) in `v`. The matrix multiplication is carried out directly using the encoding of the Householder vectors without needing to form the full matrix ![Q^T](https://www.gnu.org/software/gsl/doc/html/_images/math/7f8ec85a5eb1167b604a5766a3c0dc526d37caaf.png).

- int `gsl_linalg_QR_Qvec`(const [gsl_matrix](https://www.gnu.org/software/gsl/doc/html/vectors.html#c.gsl_matrix) * *QR*, const [gsl_vector](https://www.gnu.org/software/gsl/doc/html/vectors.html#c.gsl_vector) * *tau*, [gsl_vector](https://www.gnu.org/software/gsl/doc/html/vectors.html#c.gsl_vector) * *v*)

  This function applies the matrix ![Q](https://www.gnu.org/software/gsl/doc/html/_images/math/7fcf2b077fb63ef0b6eb5c350673266b94327847.png) encoded in the decomposition (`QR`, `tau`) to the vector `v`, storing the result ![Q v](https://www.gnu.org/software/gsl/doc/html/_images/math/e12836afcb39ec36a2bf15b9947b6f6f9446b9c9.png) in `v`. The matrix multiplication is carried out directly using the encoding of the Householder vectors without needing to form the full matrix ![Q](https://www.gnu.org/software/gsl/doc/html/_images/math/7fcf2b077fb63ef0b6eb5c350673266b94327847.png).

- int `gsl_linalg_QR_QTmat`(const [gsl_matrix](https://www.gnu.org/software/gsl/doc/html/vectors.html#c.gsl_matrix) * *QR*, const [gsl_vector](https://www.gnu.org/software/gsl/doc/html/vectors.html#c.gsl_vector) * *tau*, [gsl_matrix](https://www.gnu.org/software/gsl/doc/html/vectors.html#c.gsl_matrix) * *A*)

  This function applies the matrix ![Q^T](https://www.gnu.org/software/gsl/doc/html/_images/math/7f8ec85a5eb1167b604a5766a3c0dc526d37caaf.png) encoded in the decomposition (`QR`, `tau`) to the matrix `A`, storing the result ![Q^T A](https://www.gnu.org/software/gsl/doc/html/_images/math/28fe4dda452d99e9a1fc8d26acdcae568a123222.png) in `A`. The matrix multiplication is carried out directly using the encoding of the Householder vectors without needing to form the full matrix ![Q^T](https://www.gnu.org/software/gsl/doc/html/_images/math/7f8ec85a5eb1167b604a5766a3c0dc526d37caaf.png).

- int `gsl_linalg_QR_Rsolve`(const [gsl_matrix](https://www.gnu.org/software/gsl/doc/html/vectors.html#c.gsl_matrix) * *QR*, const [gsl_vector](https://www.gnu.org/software/gsl/doc/html/vectors.html#c.gsl_vector) * *b*, [gsl_vector](https://www.gnu.org/software/gsl/doc/html/vectors.html#c.gsl_vector) * *x*)

  This function solves the triangular system ![R x = b](https://www.gnu.org/software/gsl/doc/html/_images/math/e62de3680c5d6ae54c038051beec7f1bf1e3391b.png) for `x`. It may be useful if the product ![b' = Q^T b](https://www.gnu.org/software/gsl/doc/html/_images/math/7378e0691ac3c8650eb038dc3b45643cf42e132b.png) has already been computed using [`gsl_linalg_QR_QTvec()`](https://www.gnu.org/software/gsl/doc/html/linalg.html#c.gsl_linalg_QR_QTvec).

- int `gsl_linalg_QR_Rsvx`(const [gsl_matrix](https://www.gnu.org/software/gsl/doc/html/vectors.html#c.gsl_matrix) * *QR*, [gsl_vector](https://www.gnu.org/software/gsl/doc/html/vectors.html#c.gsl_vector) * *x*)

  This function solves the triangular system ![R x = b](https://www.gnu.org/software/gsl/doc/html/_images/math/e62de3680c5d6ae54c038051beec7f1bf1e3391b.png) for `x` in-place. On input `x` should contain the right-hand side ![b](https://www.gnu.org/software/gsl/doc/html/_images/math/88082f6488798774b9cd95a59c1ab4c6cafe909f.png) and is replaced by the solution on output. This function may be useful if the product ![b' = Q^T b](https://www.gnu.org/software/gsl/doc/html/_images/math/7378e0691ac3c8650eb038dc3b45643cf42e132b.png) has already been computed using [`gsl_linalg_QR_QTvec()`](https://www.gnu.org/software/gsl/doc/html/linalg.html#c.gsl_linalg_QR_QTvec).

- int `gsl_linalg_QR_unpack`(const [gsl_matrix](https://www.gnu.org/software/gsl/doc/html/vectors.html#c.gsl_matrix) * *QR*, const [gsl_vector](https://www.gnu.org/software/gsl/doc/html/vectors.html#c.gsl_vector) * *tau*, [gsl_matrix](https://www.gnu.org/software/gsl/doc/html/vectors.html#c.gsl_matrix) * *Q*, [gsl_matrix](https://www.gnu.org/software/gsl/doc/html/vectors.html#c.gsl_matrix) * *R*)

  This function unpacks the encoded ![QR](https://www.gnu.org/software/gsl/doc/html/_images/math/467c97789b8b4e6d169789aa163ba055c58e9223.png) decomposition (`QR`, `tau`) into the matrices `Q` and `R`, where `Q` is ![M](https://www.gnu.org/software/gsl/doc/html/_images/math/85ca40542b59eef42eec1488ce2784e307962a5f.png)-by-![M](https://www.gnu.org/software/gsl/doc/html/_images/math/85ca40542b59eef42eec1488ce2784e307962a5f.png) and `R` is ![M](https://www.gnu.org/software/gsl/doc/html/_images/math/85ca40542b59eef42eec1488ce2784e307962a5f.png)-by-![N](https://www.gnu.org/software/gsl/doc/html/_images/math/99f8a3c97c19212ca3eee7d00dddc8d68b2d1c27.png).

- int `gsl_linalg_QR_QRsolve`([gsl_matrix](https://www.gnu.org/software/gsl/doc/html/vectors.html#c.gsl_matrix) * *Q*, [gsl_matrix](https://www.gnu.org/software/gsl/doc/html/vectors.html#c.gsl_matrix) * *R*, const [gsl_vector](https://www.gnu.org/software/gsl/doc/html/vectors.html#c.gsl_vector) * *b*, [gsl_vector](https://www.gnu.org/software/gsl/doc/html/vectors.html#c.gsl_vector) * *x*)

  This function solves the system ![R x = Q^T b](https://www.gnu.org/software/gsl/doc/html/_images/math/b3fd8e1be426224ac9aa267b39a9de422b0f4c53.png) for `x`. It can be used when the ![QR](https://www.gnu.org/software/gsl/doc/html/_images/math/467c97789b8b4e6d169789aa163ba055c58e9223.png) decomposition of a matrix is available in unpacked form as (`Q`, `R`).

- int `gsl_linalg_QR_update`([gsl_matrix](https://www.gnu.org/software/gsl/doc/html/vectors.html#c.gsl_matrix) * *Q*, [gsl_matrix](https://www.gnu.org/software/gsl/doc/html/vectors.html#c.gsl_matrix) * *R*, [gsl_vector](https://www.gnu.org/software/gsl/doc/html/vectors.html#c.gsl_vector) * *w*, const [gsl_vector](https://www.gnu.org/software/gsl/doc/html/vectors.html#c.gsl_vector) * *v*)

  This function performs a rank-1 update ![w v^T](https://www.gnu.org/software/gsl/doc/html/_images/math/fbea7e56ef7a5485bbe20c50aa1792c98c0da281.png) of the ![QR](https://www.gnu.org/software/gsl/doc/html/_images/math/467c97789b8b4e6d169789aa163ba055c58e9223.png) decomposition (`Q`, `R`). The update is given by ![Q'R' = Q (R + w v^T)](https://www.gnu.org/software/gsl/doc/html/_images/math/3e8b52fb2f7fe3b1e6d250c43e1f94ed39f87327.png) where the output matrices ![Q](https://www.gnu.org/software/gsl/doc/html/_images/math/7fcf2b077fb63ef0b6eb5c350673266b94327847.png) and ![R](https://www.gnu.org/software/gsl/doc/html/_images/math/c9b67a086a723a5eec85a146053dd0a2d2378bac.png) are also orthogonal and right triangular. Note that `w` is destroyed by the update.

- int `gsl_linalg_R_solve`(const [gsl_matrix](https://www.gnu.org/software/gsl/doc/html/vectors.html#c.gsl_matrix) * *R*, const [gsl_vector](https://www.gnu.org/software/gsl/doc/html/vectors.html#c.gsl_vector) * *b*, [gsl_vector](https://www.gnu.org/software/gsl/doc/html/vectors.html#c.gsl_vector) * *x*)

  This function solves the triangular system ![R x = b](https://www.gnu.org/software/gsl/doc/html/_images/math/e62de3680c5d6ae54c038051beec7f1bf1e3391b.png) for the ![N](https://www.gnu.org/software/gsl/doc/html/_images/math/99f8a3c97c19212ca3eee7d00dddc8d68b2d1c27.png)-by-![N](https://www.gnu.org/software/gsl/doc/html/_images/math/99f8a3c97c19212ca3eee7d00dddc8d68b2d1c27.png) matrix `R`.

- int `gsl_linalg_R_svx`(const [gsl_matrix](https://www.gnu.org/software/gsl/doc/html/vectors.html#c.gsl_matrix) * *R*, [gsl_vector](https://www.gnu.org/software/gsl/doc/html/vectors.html#c.gsl_vector) * *x*)

  This function solves the triangular system ![R x = b](https://www.gnu.org/software/gsl/doc/html/_images/math/e62de3680c5d6ae54c038051beec7f1bf1e3391b.png) in-place. On input `x` should contain the right-hand side ![b](https://www.gnu.org/software/gsl/doc/html/_images/math/88082f6488798774b9cd95a59c1ab4c6cafe909f.png), which is replaced by the solution on output.



## QR Decomposition with Column Pivoting

The ![QR](https://www.gnu.org/software/gsl/doc/html/_images/math/467c97789b8b4e6d169789aa163ba055c58e9223.png) decomposition of an ![M](https://www.gnu.org/software/gsl/doc/html/_images/math/85ca40542b59eef42eec1488ce2784e307962a5f.png)-by-![N](https://www.gnu.org/software/gsl/doc/html/_images/math/99f8a3c97c19212ca3eee7d00dddc8d68b2d1c27.png) matrix ![A](https://www.gnu.org/software/gsl/doc/html/_images/math/3daf97c82cdab413e28d010c06770f224d64dc5a.png) can be extended to the rank deficient case by introducing a column permutation ![P](https://www.gnu.org/software/gsl/doc/html/_images/math/d1a0a8c0064fe114c19c8fd71ffec7a98480f707.png),

![A P = Q R](https://www.gnu.org/software/gsl/doc/html/_images/math/2bb1c5e3d46678d736c81d687bf868c1bf271163.png)

The first ![r](https://www.gnu.org/software/gsl/doc/html/_images/math/31c0f9d105c9cb974d61acbed45d11aa9e323d05.png) columns of ![Q](https://www.gnu.org/software/gsl/doc/html/_images/math/7fcf2b077fb63ef0b6eb5c350673266b94327847.png) form an orthonormal basis for the range of ![A](https://www.gnu.org/software/gsl/doc/html/_images/math/3daf97c82cdab413e28d010c06770f224d64dc5a.png) for a matrix with column rank ![r](https://www.gnu.org/software/gsl/doc/html/_images/math/31c0f9d105c9cb974d61acbed45d11aa9e323d05.png). This decomposition can also be used to convert the linear system ![A x = b](https://www.gnu.org/software/gsl/doc/html/_images/math/bf1b0fe2f3273671dbd86939aa2ca8992c570e65.png) into the triangular system ![R y = Q^T b, x = P y](https://www.gnu.org/software/gsl/doc/html/_images/math/5820cea79cc44f01d14df9d9066f5805d7e08adc.png), which can be solved by back-substitution and permutation. We denote the ![QR](https://www.gnu.org/software/gsl/doc/html/_images/math/467c97789b8b4e6d169789aa163ba055c58e9223.png) decomposition with column pivoting by ![QRP^T](https://www.gnu.org/software/gsl/doc/html/_images/math/0b15e5a9faef611e7f649cd93b9eb499aa83b70e.png) since ![A = Q R P^T](https://www.gnu.org/software/gsl/doc/html/_images/math/5cb43efb6444c61ba8f4955c10234384581e8faf.png). When ![A](https://www.gnu.org/software/gsl/doc/html/_images/math/3daf97c82cdab413e28d010c06770f224d64dc5a.png) is rank deficient with ![r = {\rm rank}(A)](https://www.gnu.org/software/gsl/doc/html/_images/math/dc287936ed1af111c95caf85e0bb987c73fbf634.png), the matrix ![R](https://www.gnu.org/software/gsl/doc/html/_images/math/c9b67a086a723a5eec85a146053dd0a2d2378bac.png) can be partitioned as

![R = \left( \begin{matrix}   R_{11} & R_{12} \\   0 & R_{22} \end{matrix} \right) \approx \left( \begin{matrix}   R_{11} & R_{12} \\   0 & 0 \end{matrix} \right)](https://www.gnu.org/software/gsl/doc/html/_images/math/8e58052c62bdf314fc3f1e22facd612118b434ad.png)

where ![R_{11}](https://www.gnu.org/software/gsl/doc/html/_images/math/9e923c840bb3f4f86d2e6203f886a587a40024cd.png) is ![r](https://www.gnu.org/software/gsl/doc/html/_images/math/31c0f9d105c9cb974d61acbed45d11aa9e323d05.png)-by-![r](https://www.gnu.org/software/gsl/doc/html/_images/math/31c0f9d105c9cb974d61acbed45d11aa9e323d05.png) and nonsingular. In this case, a *basic* least squares solution for the overdetermined system ![A x = b](https://www.gnu.org/software/gsl/doc/html/_images/math/bf1b0fe2f3273671dbd86939aa2ca8992c570e65.png) can be obtained as

![x = P \left( \begin{matrix}   R_{11}^{-1} c_1 \\   0 \end{matrix} \right)](https://www.gnu.org/software/gsl/doc/html/_images/math/95a51eeeecb37eaf9ee9d5fc66118c5b5d29e1e4.png)

where ![c_1](https://www.gnu.org/software/gsl/doc/html/_images/math/a015265e0623d5bd749e304b3a2cd1563611f055.png) consists of the first ![r](https://www.gnu.org/software/gsl/doc/html/_images/math/31c0f9d105c9cb974d61acbed45d11aa9e323d05.png) elements of ![Q^T b](https://www.gnu.org/software/gsl/doc/html/_images/math/74a156c62feffb1f2ce6765084bc63fac5c152ac.png). This basic solution is not guaranteed to be the minimum norm solution unless ![R_{12} = 0](https://www.gnu.org/software/gsl/doc/html/_images/math/4b8de8e7cd474d346343f49243fb5170a3648ed6.png) (see [Complete Orthogonal Decomposition](https://www.gnu.org/software/gsl/doc/html/linalg.html#cod)).

- int `gsl_linalg_QRPT_decomp`([gsl_matrix](https://www.gnu.org/software/gsl/doc/html/vectors.html#c.gsl_matrix) * *A*, [gsl_vector](https://www.gnu.org/software/gsl/doc/html/vectors.html#c.gsl_vector) * *tau*, [gsl_permutation](https://www.gnu.org/software/gsl/doc/html/permutation.html#c.gsl_permutation) * *p*, int * *signum*, [gsl_vector](https://www.gnu.org/software/gsl/doc/html/vectors.html#c.gsl_vector) * *norm*)

  This function factorizes the ![M](https://www.gnu.org/software/gsl/doc/html/_images/math/85ca40542b59eef42eec1488ce2784e307962a5f.png)-by-![N](https://www.gnu.org/software/gsl/doc/html/_images/math/99f8a3c97c19212ca3eee7d00dddc8d68b2d1c27.png) matrix `A` into the ![QRP^T](https://www.gnu.org/software/gsl/doc/html/_images/math/0b15e5a9faef611e7f649cd93b9eb499aa83b70e.png) decomposition ![A = Q R P^T](https://www.gnu.org/software/gsl/doc/html/_images/math/5cb43efb6444c61ba8f4955c10234384581e8faf.png). On output the diagonal and upper triangular part of the input matrix contain the matrix ![R](https://www.gnu.org/software/gsl/doc/html/_images/math/c9b67a086a723a5eec85a146053dd0a2d2378bac.png). The permutation matrix ![P](https://www.gnu.org/software/gsl/doc/html/_images/math/d1a0a8c0064fe114c19c8fd71ffec7a98480f707.png) is stored in the permutation `p`. The sign of the permutation is given by`signum`. It has the value ![(-1)^n](https://www.gnu.org/software/gsl/doc/html/_images/math/11f2b5e69df0c211c2b82609c97389c1f9d8df14.png), where ![n](https://www.gnu.org/software/gsl/doc/html/_images/math/a24554f1502cf204e7aa24a6c064962c3504de48.png) is the number of interchanges in the permutation. The vector `tau` and the columns of the lower triangular part of the matrix `A` contain the Householder coefficients and vectors which encode the orthogonal matrix `Q`. The vector `tau`must be of length ![k=\min(M,N)](https://www.gnu.org/software/gsl/doc/html/_images/math/cec598e7cfbb40428221e784f7d5bacd36bc8046.png). The matrix ![Q](https://www.gnu.org/software/gsl/doc/html/_images/math/7fcf2b077fb63ef0b6eb5c350673266b94327847.png) is related to these components by, ![Q = Q_k ... Q_2 Q_1](https://www.gnu.org/software/gsl/doc/html/_images/math/394c1fcb62e3c6b327ff0d104bb5917fa5f19b84.png) where ![Q_i = I - \tau_i v_i v_i^T](https://www.gnu.org/software/gsl/doc/html/_images/math/724f9e9a84326f1f8d512c16f51ec3ba5ff89b63.png) and ![v_i](https://www.gnu.org/software/gsl/doc/html/_images/math/78498e145a0277e0022f7d72bb541ef552325f62.png) is the Householder vector![v_i = (0,...,1,A(i+1,i),A(i+2,i),...,A(m,i))](https://www.gnu.org/software/gsl/doc/html/_images/math/797595dbdb24bbb428501fd7d7bfab49848559eb.png)This is the same storage scheme as used by LAPACK. The vector `norm` is a workspace of length `N` used for column pivoting.The algorithm used to perform the decomposition is Householder QR with column pivoting (Golub & Van Loan, “Matrix Computations”, Algorithm 5.4.1).

- int `gsl_linalg_QRPT_decomp2`(const [gsl_matrix](https://www.gnu.org/software/gsl/doc/html/vectors.html#c.gsl_matrix) * *A*, [gsl_matrix](https://www.gnu.org/software/gsl/doc/html/vectors.html#c.gsl_matrix) * *q*, [gsl_matrix](https://www.gnu.org/software/gsl/doc/html/vectors.html#c.gsl_matrix) * *r*, [gsl_vector](https://www.gnu.org/software/gsl/doc/html/vectors.html#c.gsl_vector)* *tau*, [gsl_permutation](https://www.gnu.org/software/gsl/doc/html/permutation.html#c.gsl_permutation) * *p*, int * *signum*, [gsl_vector](https://www.gnu.org/software/gsl/doc/html/vectors.html#c.gsl_vector) * *norm*)

  This function factorizes the matrix `A` into the decomposition ![A = Q R P^T](https://www.gnu.org/software/gsl/doc/html/_images/math/5cb43efb6444c61ba8f4955c10234384581e8faf.png) without modifying `A` itself and storing the output in the separate matrices `q` and `r`.

- int `gsl_linalg_QRPT_solve`(const [gsl_matrix](https://www.gnu.org/software/gsl/doc/html/vectors.html#c.gsl_matrix) * *QR*, const [gsl_vector](https://www.gnu.org/software/gsl/doc/html/vectors.html#c.gsl_vector) * *tau*, const [gsl_permutation](https://www.gnu.org/software/gsl/doc/html/permutation.html#c.gsl_permutation)* *p*, const [gsl_vector](https://www.gnu.org/software/gsl/doc/html/vectors.html#c.gsl_vector) * *b*, [gsl_vector](https://www.gnu.org/software/gsl/doc/html/vectors.html#c.gsl_vector) * *x*)

  This function solves the square system ![A x = b](https://www.gnu.org/software/gsl/doc/html/_images/math/bf1b0fe2f3273671dbd86939aa2ca8992c570e65.png) using the ![QRP^T](https://www.gnu.org/software/gsl/doc/html/_images/math/0b15e5a9faef611e7f649cd93b9eb499aa83b70e.png) decomposition of ![A](https://www.gnu.org/software/gsl/doc/html/_images/math/3daf97c82cdab413e28d010c06770f224d64dc5a.png) held in (`QR`, `tau`, `p`) which must have been computed previously by [`gsl_linalg_QRPT_decomp()`](https://www.gnu.org/software/gsl/doc/html/linalg.html#c.gsl_linalg_QRPT_decomp).

- int `gsl_linalg_QRPT_svx`(const [gsl_matrix](https://www.gnu.org/software/gsl/doc/html/vectors.html#c.gsl_matrix) * *QR*, const [gsl_vector](https://www.gnu.org/software/gsl/doc/html/vectors.html#c.gsl_vector) * *tau*, const [gsl_permutation](https://www.gnu.org/software/gsl/doc/html/permutation.html#c.gsl_permutation)* *p*, [gsl_vector](https://www.gnu.org/software/gsl/doc/html/vectors.html#c.gsl_vector) * *x*)

  This function solves the square system ![A x = b](https://www.gnu.org/software/gsl/doc/html/_images/math/bf1b0fe2f3273671dbd86939aa2ca8992c570e65.png) in-place using the ![QRP^T](https://www.gnu.org/software/gsl/doc/html/_images/math/0b15e5a9faef611e7f649cd93b9eb499aa83b70e.png) decomposition of ![A](https://www.gnu.org/software/gsl/doc/html/_images/math/3daf97c82cdab413e28d010c06770f224d64dc5a.png)held in (`QR`, `tau`, `p`). On input `x` should contain the right-hand side ![b](https://www.gnu.org/software/gsl/doc/html/_images/math/88082f6488798774b9cd95a59c1ab4c6cafe909f.png), which is replaced by the solution on output.

- int `gsl_linalg_QRPT_lssolve`(const [gsl_matrix](https://www.gnu.org/software/gsl/doc/html/vectors.html#c.gsl_matrix) * *QR*, const [gsl_vector](https://www.gnu.org/software/gsl/doc/html/vectors.html#c.gsl_vector) * *tau*, const [gsl_permutation](https://www.gnu.org/software/gsl/doc/html/permutation.html#c.gsl_permutation) * *p*, const [gsl_vector](https://www.gnu.org/software/gsl/doc/html/vectors.html#c.gsl_vector) * *b*, [gsl_vector](https://www.gnu.org/software/gsl/doc/html/vectors.html#c.gsl_vector) * *x*, [gsl_vector](https://www.gnu.org/software/gsl/doc/html/vectors.html#c.gsl_vector) * *residual*)

  This function finds the least squares solution to the overdetermined system ![A x = b](https://www.gnu.org/software/gsl/doc/html/_images/math/bf1b0fe2f3273671dbd86939aa2ca8992c570e65.png) where the matrix `A` has more rows than columns and is assumed to have full rank. The least squares solution minimizes the Euclidean norm of the residual, ![||b - A x||](https://www.gnu.org/software/gsl/doc/html/_images/math/95e033617d6ff2c39ce8fc4b303faf8fa31ba822.png). The routine requires as input the ![QR](https://www.gnu.org/software/gsl/doc/html/_images/math/467c97789b8b4e6d169789aa163ba055c58e9223.png) decomposition of ![A](https://www.gnu.org/software/gsl/doc/html/_images/math/3daf97c82cdab413e28d010c06770f224d64dc5a.png) into (`QR`, `tau`, `p`) given by [`gsl_linalg_QRPT_decomp()`](https://www.gnu.org/software/gsl/doc/html/linalg.html#c.gsl_linalg_QRPT_decomp). The solution is returned in `x`. The residual is computed as a by-product and stored in `residual`. For rank deficient matrices, [`gsl_linalg_QRPT_lssolve2()`](https://www.gnu.org/software/gsl/doc/html/linalg.html#c.gsl_linalg_QRPT_lssolve2) should be used instead.

- int `gsl_linalg_QRPT_lssolve2`(const [gsl_matrix](https://www.gnu.org/software/gsl/doc/html/vectors.html#c.gsl_matrix) * *QR*, const [gsl_vector](https://www.gnu.org/software/gsl/doc/html/vectors.html#c.gsl_vector) * *tau*, const [gsl_permutation](https://www.gnu.org/software/gsl/doc/html/permutation.html#c.gsl_permutation) * *p*, const [gsl_vector](https://www.gnu.org/software/gsl/doc/html/vectors.html#c.gsl_vector) * *b*, const size_t *rank*, [gsl_vector](https://www.gnu.org/software/gsl/doc/html/vectors.html#c.gsl_vector) * *x*, [gsl_vector](https://www.gnu.org/software/gsl/doc/html/vectors.html#c.gsl_vector) * *residual*)

  This function finds the least squares solution to the overdetermined system ![A x = b](https://www.gnu.org/software/gsl/doc/html/_images/math/bf1b0fe2f3273671dbd86939aa2ca8992c570e65.png) where the matrix `A` has more rows than columns and has rank given by the input `rank`. If the user does not know the rank of ![A](https://www.gnu.org/software/gsl/doc/html/_images/math/3daf97c82cdab413e28d010c06770f224d64dc5a.png), the routine [`gsl_linalg_QRPT_rank()`](https://www.gnu.org/software/gsl/doc/html/linalg.html#c.gsl_linalg_QRPT_rank) can be called to estimate it. The least squares solution is the so-called “basic” solution discussed above and may not be the minimum norm solution. The routine requires as input the ![QR](https://www.gnu.org/software/gsl/doc/html/_images/math/467c97789b8b4e6d169789aa163ba055c58e9223.png) decomposition of ![A](https://www.gnu.org/software/gsl/doc/html/_images/math/3daf97c82cdab413e28d010c06770f224d64dc5a.png) into (`QR`, `tau`, `p`) given by [`gsl_linalg_QRPT_decomp()`](https://www.gnu.org/software/gsl/doc/html/linalg.html#c.gsl_linalg_QRPT_decomp). The solution is returned in `x`. The residual is computed as a by-product and stored in `residual`.

- int `gsl_linalg_QRPT_QRsolve`(const [gsl_matrix](https://www.gnu.org/software/gsl/doc/html/vectors.html#c.gsl_matrix) * *Q*, const [gsl_matrix](https://www.gnu.org/software/gsl/doc/html/vectors.html#c.gsl_matrix) * *R*, const [gsl_permutation](https://www.gnu.org/software/gsl/doc/html/permutation.html#c.gsl_permutation)* *p*, const [gsl_vector](https://www.gnu.org/software/gsl/doc/html/vectors.html#c.gsl_vector) * *b*, [gsl_vector](https://www.gnu.org/software/gsl/doc/html/vectors.html#c.gsl_vector) * *x*)

  This function solves the square system ![R P^T x = Q^T b](https://www.gnu.org/software/gsl/doc/html/_images/math/daad27eb4594dcfc4dcc71a422a518a93cca1583.png) for `x`. It can be used when the ![QR](https://www.gnu.org/software/gsl/doc/html/_images/math/467c97789b8b4e6d169789aa163ba055c58e9223.png)decomposition of a matrix is available in unpacked form as (`Q`, `R`).

- int `gsl_linalg_QRPT_update`([gsl_matrix](https://www.gnu.org/software/gsl/doc/html/vectors.html#c.gsl_matrix) * *Q*, [gsl_matrix](https://www.gnu.org/software/gsl/doc/html/vectors.html#c.gsl_matrix) * *R*, const [gsl_permutation](https://www.gnu.org/software/gsl/doc/html/permutation.html#c.gsl_permutation) * *p*, [gsl_vector](https://www.gnu.org/software/gsl/doc/html/vectors.html#c.gsl_vector)* *w*, const [gsl_vector](https://www.gnu.org/software/gsl/doc/html/vectors.html#c.gsl_vector) * *v*)

  This function performs a rank-1 update ![w v^T](https://www.gnu.org/software/gsl/doc/html/_images/math/fbea7e56ef7a5485bbe20c50aa1792c98c0da281.png) of the ![QRP^T](https://www.gnu.org/software/gsl/doc/html/_images/math/0b15e5a9faef611e7f649cd93b9eb499aa83b70e.png) decomposition (`Q`, `R`, `p`). The update is given by ![Q'R' = Q (R + w v^T P)](https://www.gnu.org/software/gsl/doc/html/_images/math/6af83c25a0926680d9cef13587474edcb64b762c.png) where the output matrices ![Q'](https://www.gnu.org/software/gsl/doc/html/_images/math/421ab82e2284b0d5d3510a97288922b5cb294ec9.png) and ![R'](https://www.gnu.org/software/gsl/doc/html/_images/math/defd36efec5807157008da653d5c30a09c937394.png) are also orthogonal and right triangular. Note that `w` is destroyed by the update. The permutation `p` is not changed.

- int `gsl_linalg_QRPT_Rsolve`(const [gsl_matrix](https://www.gnu.org/software/gsl/doc/html/vectors.html#c.gsl_matrix) * *QR*, const [gsl_permutation](https://www.gnu.org/software/gsl/doc/html/permutation.html#c.gsl_permutation) * *p*, const [gsl_vector](https://www.gnu.org/software/gsl/doc/html/vectors.html#c.gsl_vector)* *b*, [gsl_vector](https://www.gnu.org/software/gsl/doc/html/vectors.html#c.gsl_vector) * *x*)

  This function solves the triangular system ![R P^T x = b](https://www.gnu.org/software/gsl/doc/html/_images/math/330b0ab7da7b3a996e017d424a2cb7bcd0a708ef.png) for the ![N](https://www.gnu.org/software/gsl/doc/html/_images/math/99f8a3c97c19212ca3eee7d00dddc8d68b2d1c27.png)-by-![N](https://www.gnu.org/software/gsl/doc/html/_images/math/99f8a3c97c19212ca3eee7d00dddc8d68b2d1c27.png) matrix ![R](https://www.gnu.org/software/gsl/doc/html/_images/math/c9b67a086a723a5eec85a146053dd0a2d2378bac.png) contained in `QR`.

- int `gsl_linalg_QRPT_Rsvx`(const [gsl_matrix](https://www.gnu.org/software/gsl/doc/html/vectors.html#c.gsl_matrix) * *QR*, const [gsl_permutation](https://www.gnu.org/software/gsl/doc/html/permutation.html#c.gsl_permutation) * *p*, [gsl_vector](https://www.gnu.org/software/gsl/doc/html/vectors.html#c.gsl_vector) * *x*)

  This function solves the triangular system ![R P^T x = b](https://www.gnu.org/software/gsl/doc/html/_images/math/330b0ab7da7b3a996e017d424a2cb7bcd0a708ef.png) in-place for the ![N](https://www.gnu.org/software/gsl/doc/html/_images/math/99f8a3c97c19212ca3eee7d00dddc8d68b2d1c27.png)-by-![N](https://www.gnu.org/software/gsl/doc/html/_images/math/99f8a3c97c19212ca3eee7d00dddc8d68b2d1c27.png) matrix ![R](https://www.gnu.org/software/gsl/doc/html/_images/math/c9b67a086a723a5eec85a146053dd0a2d2378bac.png)contained in `QR`. On input `x` should contain the right-hand side ![b](https://www.gnu.org/software/gsl/doc/html/_images/math/88082f6488798774b9cd95a59c1ab4c6cafe909f.png), which is replaced by the solution on output.

- size_t `gsl_linalg_QRPT_rank`(const [gsl_matrix](https://www.gnu.org/software/gsl/doc/html/vectors.html#c.gsl_matrix) * *QR*, const double *tol*)

  This function estimates the rank of the triangular matrix ![R](https://www.gnu.org/software/gsl/doc/html/_images/math/c9b67a086a723a5eec85a146053dd0a2d2378bac.png) contained in `QR`. The algorithm simply counts the number of diagonal elements of ![R](https://www.gnu.org/software/gsl/doc/html/_images/math/c9b67a086a723a5eec85a146053dd0a2d2378bac.png) whose absolute value is greater than the specified tolerance `tol`. If the input `tol` is negative, a default value of ![20 (M + N) eps(max(|diag(R)|))](https://www.gnu.org/software/gsl/doc/html/_images/math/808d41bc358fd7dbf9b900df55d8773316a49852.png) is used.

- int `gsl_linalg_QRPT_rcond`(const [gsl_matrix](https://www.gnu.org/software/gsl/doc/html/vectors.html#c.gsl_matrix) * *QR*, double * *rcond*, [gsl_vector](https://www.gnu.org/software/gsl/doc/html/vectors.html#c.gsl_vector) * *work*)

  This function estimates the reciprocal condition number (using the 1-norm) of the ![R](https://www.gnu.org/software/gsl/doc/html/_images/math/c9b67a086a723a5eec85a146053dd0a2d2378bac.png) factor, stored in the upper triangle of `QR`. The reciprocal condition number estimate, defined as ![1 / (||R||_1 \cdot ||R^{-1}||_1)](https://www.gnu.org/software/gsl/doc/html/_images/math/fa0df26b989ae9d43374bbb084e59ff3e01c0e9e.png), is stored in `rcond`. Additional workspace of size ![3 N](https://www.gnu.org/software/gsl/doc/html/_images/math/1b2b0ffdf792c9f396c727992622338aadd67f74.png) is required in `work`.





## Complete Orthogonal Decomposition

The complete orthogonal decomposition of a ![M](https://www.gnu.org/software/gsl/doc/html/_images/math/85ca40542b59eef42eec1488ce2784e307962a5f.png)-by-![N](https://www.gnu.org/software/gsl/doc/html/_images/math/99f8a3c97c19212ca3eee7d00dddc8d68b2d1c27.png) matrix ![A](https://www.gnu.org/software/gsl/doc/html/_images/math/3daf97c82cdab413e28d010c06770f224d64dc5a.png) is a generalization of the QR decomposition with column pivoting, given by

![A P = Q \left( \begin{matrix}   R_{11} & 0 \\   0 & 0 \end{matrix} \right) Z^T](https://www.gnu.org/software/gsl/doc/html/_images/math/49efaf7443fa880d7a3a1867592dc5fd6944e73e.png)

where ![P](https://www.gnu.org/software/gsl/doc/html/_images/math/d1a0a8c0064fe114c19c8fd71ffec7a98480f707.png) is a ![N](https://www.gnu.org/software/gsl/doc/html/_images/math/99f8a3c97c19212ca3eee7d00dddc8d68b2d1c27.png)-by-![N](https://www.gnu.org/software/gsl/doc/html/_images/math/99f8a3c97c19212ca3eee7d00dddc8d68b2d1c27.png) permutation matrix, ![Q](https://www.gnu.org/software/gsl/doc/html/_images/math/7fcf2b077fb63ef0b6eb5c350673266b94327847.png) is ![M](https://www.gnu.org/software/gsl/doc/html/_images/math/85ca40542b59eef42eec1488ce2784e307962a5f.png)-by-![M](https://www.gnu.org/software/gsl/doc/html/_images/math/85ca40542b59eef42eec1488ce2784e307962a5f.png) orthogonal, ![R_{11}](https://www.gnu.org/software/gsl/doc/html/_images/math/9e923c840bb3f4f86d2e6203f886a587a40024cd.png) is ![r](https://www.gnu.org/software/gsl/doc/html/_images/math/31c0f9d105c9cb974d61acbed45d11aa9e323d05.png)-by-![r](https://www.gnu.org/software/gsl/doc/html/_images/math/31c0f9d105c9cb974d61acbed45d11aa9e323d05.png) upper triangular, with ![r = {\rm rank}(A)](https://www.gnu.org/software/gsl/doc/html/_images/math/dc287936ed1af111c95caf85e0bb987c73fbf634.png), and ![Z](https://www.gnu.org/software/gsl/doc/html/_images/math/1c5cb7f4fbed6ab88f22a1444e3dc273fd24860a.png) is ![N](https://www.gnu.org/software/gsl/doc/html/_images/math/99f8a3c97c19212ca3eee7d00dddc8d68b2d1c27.png)-by-![N](https://www.gnu.org/software/gsl/doc/html/_images/math/99f8a3c97c19212ca3eee7d00dddc8d68b2d1c27.png) orthogonal. If ![A](https://www.gnu.org/software/gsl/doc/html/_images/math/3daf97c82cdab413e28d010c06770f224d64dc5a.png) has full rank, then ![R_{11} = R](https://www.gnu.org/software/gsl/doc/html/_images/math/61895703174509421882cef5c008fc82dc564816.png), ![Z = I](https://www.gnu.org/software/gsl/doc/html/_images/math/fee79fee0e3f0cf8a329d9f7f14cfff0d662021e.png) and this reduces to the QR decomposition with column pivoting.

For a rank deficient least squares problem, ![\min_x{|| b - Ax||^2}](https://www.gnu.org/software/gsl/doc/html/_images/math/23a0daf31338d0937589e561366f3ed1cf51e726.png), the solution vector ![x](https://www.gnu.org/software/gsl/doc/html/_images/math/1c5715c8edd0bbc330baa0f0494f1f7951619754.png) is not unique. However if we further require that ![||x||^2](https://www.gnu.org/software/gsl/doc/html/_images/math/4047fbb464a213c6b25a3d5871315f3f5a4e3189.png) is minimized, then the complete orthogonal decomposition gives us the ability to compute the unique minimum norm solution, which is given by

![x = P Z \left( \begin{matrix}   R_{11}^{-1} c_1 \\   0 \end{matrix} \right)](https://www.gnu.org/software/gsl/doc/html/_images/math/f2841d79068a307c6b27b90d1b51a550da67d5c4.png)

and the vector ![c_1](https://www.gnu.org/software/gsl/doc/html/_images/math/a015265e0623d5bd749e304b3a2cd1563611f055.png) is the first ![r](https://www.gnu.org/software/gsl/doc/html/_images/math/31c0f9d105c9cb974d61acbed45d11aa9e323d05.png) elements of ![Q^T b](https://www.gnu.org/software/gsl/doc/html/_images/math/74a156c62feffb1f2ce6765084bc63fac5c152ac.png).

The COD also enables a straightforward solution of regularized least squares problems in Tikhonov standard form, written as

![\min_x ||b - A x||^2 + \lambda^2 ||x||^2](https://www.gnu.org/software/gsl/doc/html/_images/math/fb22a169fd53274ab0508258be74a66d7aff4e2e.png)

where ![\lambda > 0](https://www.gnu.org/software/gsl/doc/html/_images/math/6e9b95f084b43d5b271bb562bb305dd15e963ad1.png) is a regularization parameter which represents a tradeoff between minimizing the residual norm ![||b - A x||](https://www.gnu.org/software/gsl/doc/html/_images/math/95e033617d6ff2c39ce8fc4b303faf8fa31ba822.png) and the solution norm ![||x||](https://www.gnu.org/software/gsl/doc/html/_images/math/04251f512a3974b1d34ec448936a9d7d8734a982.png). For this system, the solution is given by

![x = P Z \left( \begin{matrix}   y_1 \\   0 \end{matrix} \right)](https://www.gnu.org/software/gsl/doc/html/_images/math/cbf56262639736fb5e63ca351c8f7b5f0ef2b175.png)

where ![y_1](https://www.gnu.org/software/gsl/doc/html/_images/math/490f3f4cf44a8fca8b9258bcbadba4fdc5cfc2b2.png) is a vector of length ![r](https://www.gnu.org/software/gsl/doc/html/_images/math/31c0f9d105c9cb974d61acbed45d11aa9e323d05.png) which is found by solving

![\left( \begin{matrix}   R_{11} \\   \lambda I_r \end{matrix} \right) y_1 = \left( \begin{matrix}   c_1 \\   0 \end{matrix} \right)](https://www.gnu.org/software/gsl/doc/html/_images/math/95a614fa7a9e56a23ae17be6453ddea184ec6676.png)

and ![c_1](https://www.gnu.org/software/gsl/doc/html/_images/math/a015265e0623d5bd749e304b3a2cd1563611f055.png) is defined above. The equation above can be solved efficiently for different values of ![\lambda](https://www.gnu.org/software/gsl/doc/html/_images/math/458d03720fd1c320b48535a3bfc70fa5c1711b78.png)using QR factorizations of the left hand side matrix.

- int `gsl_linalg_COD_decomp`([gsl_matrix](https://www.gnu.org/software/gsl/doc/html/vectors.html#c.gsl_matrix) * *A*, [gsl_vector](https://www.gnu.org/software/gsl/doc/html/vectors.html#c.gsl_vector) * *tau_Q*, [gsl_vector](https://www.gnu.org/software/gsl/doc/html/vectors.html#c.gsl_vector) * *tau_Z*, [gsl_permutation](https://www.gnu.org/software/gsl/doc/html/permutation.html#c.gsl_permutation) * *p*, size_t * *rank*, [gsl_vector](https://www.gnu.org/software/gsl/doc/html/vectors.html#c.gsl_vector) * *work*)

- int `gsl_linalg_COD_decomp_e`([gsl_matrix](https://www.gnu.org/software/gsl/doc/html/vectors.html#c.gsl_matrix) * *A*, [gsl_vector](https://www.gnu.org/software/gsl/doc/html/vectors.html#c.gsl_vector) * *tau_Q*, [gsl_vector](https://www.gnu.org/software/gsl/doc/html/vectors.html#c.gsl_vector) * *tau_Z*, [gsl_permutation](https://www.gnu.org/software/gsl/doc/html/permutation.html#c.gsl_permutation) * *p*, double *tol*, size_t * *rank*, [gsl_vector](https://www.gnu.org/software/gsl/doc/html/vectors.html#c.gsl_vector) * *work*)

  These functions factor the ![M](https://www.gnu.org/software/gsl/doc/html/_images/math/85ca40542b59eef42eec1488ce2784e307962a5f.png)-by-![N](https://www.gnu.org/software/gsl/doc/html/_images/math/99f8a3c97c19212ca3eee7d00dddc8d68b2d1c27.png) matrix `A` into the decomposition ![A = Q R Z P^T](https://www.gnu.org/software/gsl/doc/html/_images/math/ae21d260dd3811a142059b4e50d74306fa916e7c.png). The rank of `A` is computed as the number of diagonal elements of ![R](https://www.gnu.org/software/gsl/doc/html/_images/math/c9b67a086a723a5eec85a146053dd0a2d2378bac.png) greater than the tolerance `tol`and output in `rank`. If `tol` is not specified, a default value is used (see [`gsl_linalg_QRPT_rank()`](https://www.gnu.org/software/gsl/doc/html/linalg.html#c.gsl_linalg_QRPT_rank)). On output, the permutation matrix ![P](https://www.gnu.org/software/gsl/doc/html/_images/math/d1a0a8c0064fe114c19c8fd71ffec7a98480f707.png) is stored in `p`. The matrix ![R_{11}](https://www.gnu.org/software/gsl/doc/html/_images/math/9e923c840bb3f4f86d2e6203f886a587a40024cd.png)is stored in the upper `rank`-by-`rank` block of `A`. The matrices ![Q](https://www.gnu.org/software/gsl/doc/html/_images/math/7fcf2b077fb63ef0b6eb5c350673266b94327847.png) and ![Z](https://www.gnu.org/software/gsl/doc/html/_images/math/1c5cb7f4fbed6ab88f22a1444e3dc273fd24860a.png) are encoded in packed storage in `A` on output. The vectors `tau_Q` and `tau_Z` contain the Householder scalars corresponding to the matrices ![Q](https://www.gnu.org/software/gsl/doc/html/_images/math/7fcf2b077fb63ef0b6eb5c350673266b94327847.png) and ![Z](https://www.gnu.org/software/gsl/doc/html/_images/math/1c5cb7f4fbed6ab88f22a1444e3dc273fd24860a.png) respectively and must be of length ![k = \min(M,N)](https://www.gnu.org/software/gsl/doc/html/_images/math/31c568301eed089942f705bb0073598abab7abe5.png). The vector `work` is additional workspace of length ![N](https://www.gnu.org/software/gsl/doc/html/_images/math/99f8a3c97c19212ca3eee7d00dddc8d68b2d1c27.png).

- int `gsl_linalg_COD_lssolve`(const [gsl_matrix](https://www.gnu.org/software/gsl/doc/html/vectors.html#c.gsl_matrix) * *QRZT*, const [gsl_vector](https://www.gnu.org/software/gsl/doc/html/vectors.html#c.gsl_vector) * *tau_Q*, const [gsl_vector](https://www.gnu.org/software/gsl/doc/html/vectors.html#c.gsl_vector)* *tau_Z*, const [gsl_permutation](https://www.gnu.org/software/gsl/doc/html/permutation.html#c.gsl_permutation) * *p*, const size_t *rank*, const [gsl_vector](https://www.gnu.org/software/gsl/doc/html/vectors.html#c.gsl_vector) * *b*, [gsl_vector](https://www.gnu.org/software/gsl/doc/html/vectors.html#c.gsl_vector) * *x*, [gsl_vector](https://www.gnu.org/software/gsl/doc/html/vectors.html#c.gsl_vector)* *residual*)

  This function finds the unique minimum norm least squares solution to the overdetermined system ![A x = b](https://www.gnu.org/software/gsl/doc/html/_images/math/bf1b0fe2f3273671dbd86939aa2ca8992c570e65.png) where the matrix `A` has more rows than columns. The least squares solution minimizes the Euclidean norm of the residual, ![||b - A x||](https://www.gnu.org/software/gsl/doc/html/_images/math/95e033617d6ff2c39ce8fc4b303faf8fa31ba822.png) as well as the norm of the solution ![||x||](https://www.gnu.org/software/gsl/doc/html/_images/math/04251f512a3974b1d34ec448936a9d7d8734a982.png). The routine requires as input the ![QRZT](https://www.gnu.org/software/gsl/doc/html/_images/math/124af73ee0acec99d8788164e80e6a5f1554f932.png) decomposition of ![A](https://www.gnu.org/software/gsl/doc/html/_images/math/3daf97c82cdab413e28d010c06770f224d64dc5a.png) into (`QRZT`, `tau_Q`, `tau_Z`, `p`, `rank`) given by [`gsl_linalg_COD_decomp()`](https://www.gnu.org/software/gsl/doc/html/linalg.html#c.gsl_linalg_COD_decomp). The solution is returned in `x`. The residual, ![b - Ax](https://www.gnu.org/software/gsl/doc/html/_images/math/bf719bb5387c6d21defa162d4f1f38e64c1d89dd.png), is computed as a by-product and stored in `residual`.

- int `gsl_linalg_COD_lssolve2`(const double *lambda*, const [gsl_matrix](https://www.gnu.org/software/gsl/doc/html/vectors.html#c.gsl_matrix) * *QRZT*, const [gsl_vector](https://www.gnu.org/software/gsl/doc/html/vectors.html#c.gsl_vector)* *tau_Q*, const [gsl_vector](https://www.gnu.org/software/gsl/doc/html/vectors.html#c.gsl_vector) * *tau_Z*, const [gsl_permutation](https://www.gnu.org/software/gsl/doc/html/permutation.html#c.gsl_permutation) * *p*, const size_t *rank*, const [gsl_vector](https://www.gnu.org/software/gsl/doc/html/vectors.html#c.gsl_vector) * *b*, [gsl_vector](https://www.gnu.org/software/gsl/doc/html/vectors.html#c.gsl_vector) * *x*, [gsl_vector](https://www.gnu.org/software/gsl/doc/html/vectors.html#c.gsl_vector) * *residual*, [gsl_matrix](https://www.gnu.org/software/gsl/doc/html/vectors.html#c.gsl_matrix) * *S*, [gsl_vector](https://www.gnu.org/software/gsl/doc/html/vectors.html#c.gsl_vector) * *work*)

  This function finds the solution to the regularized least squares problem in Tikhonov standard form, ![\min_x ||b - Ax||^2 + \lambda^2 ||x||^2](https://www.gnu.org/software/gsl/doc/html/_images/math/b553ed4154bd1f3d747b1880741c930cebff0bab.png). The routine requires as input the ![QRZT](https://www.gnu.org/software/gsl/doc/html/_images/math/124af73ee0acec99d8788164e80e6a5f1554f932.png) decomposition of ![A](https://www.gnu.org/software/gsl/doc/html/_images/math/3daf97c82cdab413e28d010c06770f224d64dc5a.png) into (`QRZT`, `tau_Q`, `tau_Z`, `p`, `rank`) given by [`gsl_linalg_COD_decomp()`](https://www.gnu.org/software/gsl/doc/html/linalg.html#c.gsl_linalg_COD_decomp). The parameter ![\lambda](https://www.gnu.org/software/gsl/doc/html/_images/math/458d03720fd1c320b48535a3bfc70fa5c1711b78.png)is supplied in `lambda`. The solution is returned in `x`. The residual, ![b - Ax](https://www.gnu.org/software/gsl/doc/html/_images/math/bf719bb5387c6d21defa162d4f1f38e64c1d89dd.png), is stored in `residual` on output. `S` is additional workspace of size `rank`-by-`rank`. `work` is additional workspace of length `rank`.

- int `gsl_linalg_COD_unpack`(const [gsl_matrix](https://www.gnu.org/software/gsl/doc/html/vectors.html#c.gsl_matrix) * *QRZT*, const [gsl_vector](https://www.gnu.org/software/gsl/doc/html/vectors.html#c.gsl_vector) * *tau_Q*, const [gsl_vector](https://www.gnu.org/software/gsl/doc/html/vectors.html#c.gsl_vector)* *tau_Z*, const size_t *rank*, [gsl_matrix](https://www.gnu.org/software/gsl/doc/html/vectors.html#c.gsl_matrix) * *Q*, [gsl_matrix](https://www.gnu.org/software/gsl/doc/html/vectors.html#c.gsl_matrix) * *R*, [gsl_matrix](https://www.gnu.org/software/gsl/doc/html/vectors.html#c.gsl_matrix) * *Z*)

  This function unpacks the encoded ![QRZT](https://www.gnu.org/software/gsl/doc/html/_images/math/124af73ee0acec99d8788164e80e6a5f1554f932.png) decomposition (`QRZT`, `tau_Q`, `tau_Z`, `rank`) into the matrices `Q`, `R`, and `Z`, where `Q` is ![M](https://www.gnu.org/software/gsl/doc/html/_images/math/85ca40542b59eef42eec1488ce2784e307962a5f.png)-by-![M](https://www.gnu.org/software/gsl/doc/html/_images/math/85ca40542b59eef42eec1488ce2784e307962a5f.png), `R` is ![M](https://www.gnu.org/software/gsl/doc/html/_images/math/85ca40542b59eef42eec1488ce2784e307962a5f.png)-by-![N](https://www.gnu.org/software/gsl/doc/html/_images/math/99f8a3c97c19212ca3eee7d00dddc8d68b2d1c27.png), and `Z` is ![N](https://www.gnu.org/software/gsl/doc/html/_images/math/99f8a3c97c19212ca3eee7d00dddc8d68b2d1c27.png)-by-![N](https://www.gnu.org/software/gsl/doc/html/_images/math/99f8a3c97c19212ca3eee7d00dddc8d68b2d1c27.png).

- int `gsl_linalg_COD_matZ`(const [gsl_matrix](https://www.gnu.org/software/gsl/doc/html/vectors.html#c.gsl_matrix) * *QRZT*, const [gsl_vector](https://www.gnu.org/software/gsl/doc/html/vectors.html#c.gsl_vector) * *tau_Z*, const size_t *rank*, [gsl_matrix](https://www.gnu.org/software/gsl/doc/html/vectors.html#c.gsl_matrix) * *A*, [gsl_vector](https://www.gnu.org/software/gsl/doc/html/vectors.html#c.gsl_vector) * *work*)

  This function multiplies the input matrix `A` on the right by `Z`, ![A' = A Z](https://www.gnu.org/software/gsl/doc/html/_images/math/7bce2dcde0e616eb28a34ae8d4e5c7cc009f9155.png) using the encoded ![QRZT](https://www.gnu.org/software/gsl/doc/html/_images/math/124af73ee0acec99d8788164e80e6a5f1554f932.png) decomposition (`QRZT`, `tau_Z`, `rank`). `A` must have ![N](https://www.gnu.org/software/gsl/doc/html/_images/math/99f8a3c97c19212ca3eee7d00dddc8d68b2d1c27.png) columns but may have any number of rows. Additional workspace of length ![M](https://www.gnu.org/software/gsl/doc/html/_images/math/85ca40542b59eef42eec1488ce2784e307962a5f.png) is provided in `work`.



## Singular Value Decomposition

A general rectangular ![M](https://www.gnu.org/software/gsl/doc/html/_images/math/85ca40542b59eef42eec1488ce2784e307962a5f.png)-by-![N](https://www.gnu.org/software/gsl/doc/html/_images/math/99f8a3c97c19212ca3eee7d00dddc8d68b2d1c27.png) matrix ![A](https://www.gnu.org/software/gsl/doc/html/_images/math/3daf97c82cdab413e28d010c06770f224d64dc5a.png) has a singular value decomposition (SVD) into the product of an ![M](https://www.gnu.org/software/gsl/doc/html/_images/math/85ca40542b59eef42eec1488ce2784e307962a5f.png)-by-![N](https://www.gnu.org/software/gsl/doc/html/_images/math/99f8a3c97c19212ca3eee7d00dddc8d68b2d1c27.png) orthogonal matrix ![U](https://www.gnu.org/software/gsl/doc/html/_images/math/3d0e8e66abfde1273e99e86e8a090a7073d842e3.png), an ![N](https://www.gnu.org/software/gsl/doc/html/_images/math/99f8a3c97c19212ca3eee7d00dddc8d68b2d1c27.png)-by-![N](https://www.gnu.org/software/gsl/doc/html/_images/math/99f8a3c97c19212ca3eee7d00dddc8d68b2d1c27.png) diagonal matrix of singular values ![S](https://www.gnu.org/software/gsl/doc/html/_images/math/4687d5918c214ce691dee3322ecf3d2046b28e1f.png) and the transpose of an ![N](https://www.gnu.org/software/gsl/doc/html/_images/math/99f8a3c97c19212ca3eee7d00dddc8d68b2d1c27.png)-by-![N](https://www.gnu.org/software/gsl/doc/html/_images/math/99f8a3c97c19212ca3eee7d00dddc8d68b2d1c27.png) orthogonal square matrix ![V](https://www.gnu.org/software/gsl/doc/html/_images/math/632777989bf30753abbbf45a4a9b07ebf0fd36ad.png),

![A = U S V^T](https://www.gnu.org/software/gsl/doc/html/_images/math/c3b6630f24aad355776b58fd1ef99f4dfc4c1b1d.png)

The singular values ![\sigma_i = S_{ii}](https://www.gnu.org/software/gsl/doc/html/_images/math/188242c14c1725156bff83dd3d263e7f46b30544.png) are all non-negative and are generally chosen to form a non-increasing sequence

![\sigma_1 \ge \sigma_2 \ge ... \ge \sigma_N \ge 0](https://www.gnu.org/software/gsl/doc/html/_images/math/9aba5b6447a953fa667bf397788ebe4f69904239.png)

The singular value decomposition of a matrix has many practical uses. The condition number of the matrix is given by the ratio of the largest singular value to the smallest singular value. The presence of a zero singular value indicates that the matrix is singular. The number of non-zero singular values indicates the rank of the matrix. In practice singular value decomposition of a rank-deficient matrix will not produce exact zeroes for singular values, due to finite numerical precision. Small singular values should be edited by choosing a suitable tolerance.

For a rank-deficient matrix, the null space of ![A](https://www.gnu.org/software/gsl/doc/html/_images/math/3daf97c82cdab413e28d010c06770f224d64dc5a.png) is given by the columns of ![V](https://www.gnu.org/software/gsl/doc/html/_images/math/632777989bf30753abbbf45a4a9b07ebf0fd36ad.png) corresponding to the zero singular values. Similarly, the range of ![A](https://www.gnu.org/software/gsl/doc/html/_images/math/3daf97c82cdab413e28d010c06770f224d64dc5a.png) is given by columns of ![U](https://www.gnu.org/software/gsl/doc/html/_images/math/3d0e8e66abfde1273e99e86e8a090a7073d842e3.png) corresponding to the non-zero singular values.

Note that the routines here compute the “thin” version of the SVD with ![U](https://www.gnu.org/software/gsl/doc/html/_images/math/3d0e8e66abfde1273e99e86e8a090a7073d842e3.png) as ![M](https://www.gnu.org/software/gsl/doc/html/_images/math/85ca40542b59eef42eec1488ce2784e307962a5f.png)-by-![N](https://www.gnu.org/software/gsl/doc/html/_images/math/99f8a3c97c19212ca3eee7d00dddc8d68b2d1c27.png) orthogonal matrix. This allows in-place computation and is the most commonly-used form in practice. Mathematically, the “full” SVD is defined with ![U](https://www.gnu.org/software/gsl/doc/html/_images/math/3d0e8e66abfde1273e99e86e8a090a7073d842e3.png) as an ![M](https://www.gnu.org/software/gsl/doc/html/_images/math/85ca40542b59eef42eec1488ce2784e307962a5f.png)-by-![M](https://www.gnu.org/software/gsl/doc/html/_images/math/85ca40542b59eef42eec1488ce2784e307962a5f.png) orthogonal matrix and ![S](https://www.gnu.org/software/gsl/doc/html/_images/math/4687d5918c214ce691dee3322ecf3d2046b28e1f.png) as an ![M](https://www.gnu.org/software/gsl/doc/html/_images/math/85ca40542b59eef42eec1488ce2784e307962a5f.png)-by-![N](https://www.gnu.org/software/gsl/doc/html/_images/math/99f8a3c97c19212ca3eee7d00dddc8d68b2d1c27.png) diagonal matrix (with additional rows of zeros).

- int `gsl_linalg_SV_decomp`([gsl_matrix](https://www.gnu.org/software/gsl/doc/html/vectors.html#c.gsl_matrix) * *A*, [gsl_matrix](https://www.gnu.org/software/gsl/doc/html/vectors.html#c.gsl_matrix) * *V*, [gsl_vector](https://www.gnu.org/software/gsl/doc/html/vectors.html#c.gsl_vector) * *S*, [gsl_vector](https://www.gnu.org/software/gsl/doc/html/vectors.html#c.gsl_vector) * *work*)

  This function factorizes the ![M](https://www.gnu.org/software/gsl/doc/html/_images/math/85ca40542b59eef42eec1488ce2784e307962a5f.png)-by-![N](https://www.gnu.org/software/gsl/doc/html/_images/math/99f8a3c97c19212ca3eee7d00dddc8d68b2d1c27.png) matrix `A` into the singular value decomposition ![A = U S V^T](https://www.gnu.org/software/gsl/doc/html/_images/math/d26e5eb59422ca146ae6fa8ae0d5ec73b8573868.png) for ![M \ge N](https://www.gnu.org/software/gsl/doc/html/_images/math/736ad7d9a01455eddb311df76669dff8b8af8f4b.png). On output the matrix `A` is replaced by ![U](https://www.gnu.org/software/gsl/doc/html/_images/math/3d0e8e66abfde1273e99e86e8a090a7073d842e3.png). The diagonal elements of the singular value matrix ![S](https://www.gnu.org/software/gsl/doc/html/_images/math/4687d5918c214ce691dee3322ecf3d2046b28e1f.png) are stored in the vector `S`. The singular values are non-negative and form a non-increasing sequence from ![S_1](https://www.gnu.org/software/gsl/doc/html/_images/math/91f16fda59db71e0ee05e4faa4c294a190c6dc94.png) to ![S_N](https://www.gnu.org/software/gsl/doc/html/_images/math/287ff1f99032ab690da44e96df1497403cb6812b.png). The matrix `V` contains the elements of ![V](https://www.gnu.org/software/gsl/doc/html/_images/math/632777989bf30753abbbf45a4a9b07ebf0fd36ad.png)in untransposed form. To form the product ![U S V^T](https://www.gnu.org/software/gsl/doc/html/_images/math/c6adcbc5abf13bdae8f0ac1eb6e47f7092db707d.png) it is necessary to take the transpose of `V`. A workspace of length `N` is required in `work`.This routine uses the Golub-Reinsch SVD algorithm.

- int `gsl_linalg_SV_decomp_mod`([gsl_matrix](https://www.gnu.org/software/gsl/doc/html/vectors.html#c.gsl_matrix) * *A*, [gsl_matrix](https://www.gnu.org/software/gsl/doc/html/vectors.html#c.gsl_matrix) * *X*, [gsl_matrix](https://www.gnu.org/software/gsl/doc/html/vectors.html#c.gsl_matrix) * *V*, [gsl_vector](https://www.gnu.org/software/gsl/doc/html/vectors.html#c.gsl_vector) * *S*, [gsl_vector](https://www.gnu.org/software/gsl/doc/html/vectors.html#c.gsl_vector) * *work*)

  This function computes the SVD using the modified Golub-Reinsch algorithm, which is faster for ![M \gg N](https://www.gnu.org/software/gsl/doc/html/_images/math/89fa2bf27b4ab76c1e34be04de14e3310a997c36.png). It requires the vector `work` of length `N` and the ![N](https://www.gnu.org/software/gsl/doc/html/_images/math/99f8a3c97c19212ca3eee7d00dddc8d68b2d1c27.png)-by-![N](https://www.gnu.org/software/gsl/doc/html/_images/math/99f8a3c97c19212ca3eee7d00dddc8d68b2d1c27.png) matrix `X` as additional working space.



- int `gsl_linalg_SV_decomp_jacobi`([gsl_matrix](https://www.gnu.org/software/gsl/doc/html/vectors.html#c.gsl_matrix) * *A*, [gsl_matrix](https://www.gnu.org/software/gsl/doc/html/vectors.html#c.gsl_matrix) * *V*, [gsl_vector](https://www.gnu.org/software/gsl/doc/html/vectors.html#c.gsl_vector) * *S*)

  This function computes the SVD of the ![M](https://www.gnu.org/software/gsl/doc/html/_images/math/85ca40542b59eef42eec1488ce2784e307962a5f.png)-by-![N](https://www.gnu.org/software/gsl/doc/html/_images/math/99f8a3c97c19212ca3eee7d00dddc8d68b2d1c27.png) matrix `A` using one-sided Jacobi orthogonalization for ![M \ge N](https://www.gnu.org/software/gsl/doc/html/_images/math/736ad7d9a01455eddb311df76669dff8b8af8f4b.png). The Jacobi method can compute singular values to higher relative accuracy than Golub-Reinsch algorithms (see references for details).

- int `gsl_linalg_SV_solve`(const [gsl_matrix](https://www.gnu.org/software/gsl/doc/html/vectors.html#c.gsl_matrix) * *U*, const [gsl_matrix](https://www.gnu.org/software/gsl/doc/html/vectors.html#c.gsl_matrix) * *V*, const [gsl_vector](https://www.gnu.org/software/gsl/doc/html/vectors.html#c.gsl_vector) * *S*, const [gsl_vector](https://www.gnu.org/software/gsl/doc/html/vectors.html#c.gsl_vector) * *b*, [gsl_vector](https://www.gnu.org/software/gsl/doc/html/vectors.html#c.gsl_vector) * *x*)

  This function solves the system ![A x = b](https://www.gnu.org/software/gsl/doc/html/_images/math/bf1b0fe2f3273671dbd86939aa2ca8992c570e65.png) using the singular value decomposition (`U`, `S`, `V`) of ![A](https://www.gnu.org/software/gsl/doc/html/_images/math/3daf97c82cdab413e28d010c06770f224d64dc5a.png) which must have been computed previously with [`gsl_linalg_SV_decomp()`](https://www.gnu.org/software/gsl/doc/html/linalg.html#c.gsl_linalg_SV_decomp).Only non-zero singular values are used in computing the solution. The parts of the solution corresponding to singular values of zero are ignored. Other singular values can be edited out by setting them to zero before calling this function.In the over-determined case where `A` has more rows than columns the system is solved in the least squares sense, returning the solution `x` which minimizes ![||A x - b||_2](https://www.gnu.org/software/gsl/doc/html/_images/math/9ef92fe409cbcfbee25080259183237769ad31eb.png).

- int `gsl_linalg_SV_leverage`(const [gsl_matrix](https://www.gnu.org/software/gsl/doc/html/vectors.html#c.gsl_matrix) * *U*, [gsl_vector](https://www.gnu.org/software/gsl/doc/html/vectors.html#c.gsl_vector) * *h*)

  This function computes the statistical leverage values ![h_i](https://www.gnu.org/software/gsl/doc/html/_images/math/881774e10df53f0390bc0f793f8454163afe3d57.png) of a matrix ![A](https://www.gnu.org/software/gsl/doc/html/_images/math/3daf97c82cdab413e28d010c06770f224d64dc5a.png) using its singular value decomposition (`U`, `S`, `V`) previously computed with [`gsl_linalg_SV_decomp()`](https://www.gnu.org/software/gsl/doc/html/linalg.html#c.gsl_linalg_SV_decomp). ![h_i](https://www.gnu.org/software/gsl/doc/html/_images/math/881774e10df53f0390bc0f793f8454163afe3d57.png) are the diagonal values of the matrix ![A (A^T A)^{-1} A^T](https://www.gnu.org/software/gsl/doc/html/_images/math/5b44ebe09f726ce0692ca95c52a9224faa92c2ec.png) and depend only on the matrix `U` which is the input to this function.





## Cholesky Decomposition

A symmetric, positive definite square matrix ![A](https://www.gnu.org/software/gsl/doc/html/_images/math/3daf97c82cdab413e28d010c06770f224d64dc5a.png) has a Cholesky decomposition into a product of a lower triangular matrix ![L](https://www.gnu.org/software/gsl/doc/html/_images/math/82359a7c56cc02d70a5ed0e98df631577f7d172b.png) and its transpose ![L^T](https://www.gnu.org/software/gsl/doc/html/_images/math/36fd14d4f6a032c0d439613fa95d6e11c291bf09.png),

![A = L L^T](https://www.gnu.org/software/gsl/doc/html/_images/math/ad02cee1d01fac6be3b6dc84a322d01dd973d326.png)

This is sometimes referred to as taking the square-root of a matrix. The Cholesky decomposition can only be carried out when all the eigenvalues of the matrix are positive. This decomposition can be used to convert the linear system ![A x = b](https://www.gnu.org/software/gsl/doc/html/_images/math/bf1b0fe2f3273671dbd86939aa2ca8992c570e65.png) into a pair of triangular systems (![L y = b](https://www.gnu.org/software/gsl/doc/html/_images/math/18e896ebb9553db9b0e20dd3621e7cbabe73b642.png), ![L^T x = y](https://www.gnu.org/software/gsl/doc/html/_images/math/e58fc969cd5c8da4aeb1f2c1abce0e36921eeace.png)), which can be solved by forward and back-substitution.

If the matrix ![A](https://www.gnu.org/software/gsl/doc/html/_images/math/3daf97c82cdab413e28d010c06770f224d64dc5a.png) is near singular, it is sometimes possible to reduce the condition number and recover a more accurate solution vector ![x](https://www.gnu.org/software/gsl/doc/html/_images/math/1c5715c8edd0bbc330baa0f0494f1f7951619754.png) by scaling as

![\left( S A S \right) \left( S^{-1} x \right) = S b](https://www.gnu.org/software/gsl/doc/html/_images/math/9b9feed87846567b926d19b7f041f8bba690d369.png)

where ![S](https://www.gnu.org/software/gsl/doc/html/_images/math/4687d5918c214ce691dee3322ecf3d2046b28e1f.png) is a diagonal matrix whose elements are given by ![S_{ii} = 1/\sqrt{A_{ii}}](https://www.gnu.org/software/gsl/doc/html/_images/math/f205500b66dfb801c12b4c4d1746a0a1ce4abb41.png). This scaling is also known as Jacobi preconditioning. There are routines below to solve both the scaled and unscaled systems.

- int `gsl_linalg_cholesky_decomp1`([gsl_matrix](https://www.gnu.org/software/gsl/doc/html/vectors.html#c.gsl_matrix) * *A*)

- int `gsl_linalg_complex_cholesky_decomp`(gsl_matrix_complex * *A*)

  These functions factorize the symmetric, positive-definite square matrix `A` into the Cholesky decomposition ![A = L L^T](https://www.gnu.org/software/gsl/doc/html/_images/math/90b03e70f53691e2440071a0e4f5da1ae24e268d.png) (or ![A = L L^{\dagger}](https://www.gnu.org/software/gsl/doc/html/_images/math/4f386261cdd7361ebd157936843da613aa1d3cd1.png) for the complex case). On input, the values from the diagonal and lower-triangular part of the matrix `A` are used (the upper triangular part is ignored). On output the diagonal and lower triangular part of the input matrix `A` contain the matrix ![L](https://www.gnu.org/software/gsl/doc/html/_images/math/82359a7c56cc02d70a5ed0e98df631577f7d172b.png), while the upper triangular part is unmodified. If the matrix is not positive-definite then the decomposition will fail, returning the error code [`GSL_EDOM`](https://www.gnu.org/software/gsl/doc/html/err.html#c.GSL_EDOM).When testing whether a matrix is positive-definite, disable the error handler first to avoid triggering an error.

- int `gsl_linalg_cholesky_decomp`([gsl_matrix](https://www.gnu.org/software/gsl/doc/html/vectors.html#c.gsl_matrix) * *A*)

  This function is now deprecated and is provided only for backward compatibility.

- int `gsl_linalg_cholesky_solve`(const [gsl_matrix](https://www.gnu.org/software/gsl/doc/html/vectors.html#c.gsl_matrix) * *cholesky*, const [gsl_vector](https://www.gnu.org/software/gsl/doc/html/vectors.html#c.gsl_vector) * *b*, [gsl_vector](https://www.gnu.org/software/gsl/doc/html/vectors.html#c.gsl_vector) * *x*)

- int `gsl_linalg_complex_cholesky_solve`(const gsl_matrix_complex * *cholesky*, const gsl_vector_complex * *b*, gsl_vector_complex * *x*)

  These functions solve the system ![A x = b](https://www.gnu.org/software/gsl/doc/html/_images/math/bf1b0fe2f3273671dbd86939aa2ca8992c570e65.png) using the Cholesky decomposition of ![A](https://www.gnu.org/software/gsl/doc/html/_images/math/3daf97c82cdab413e28d010c06770f224d64dc5a.png) held in the matrix `cholesky` which must have been previously computed by [`gsl_linalg_cholesky_decomp()`](https://www.gnu.org/software/gsl/doc/html/linalg.html#c.gsl_linalg_cholesky_decomp) or [`gsl_linalg_complex_cholesky_decomp()`](https://www.gnu.org/software/gsl/doc/html/linalg.html#c.gsl_linalg_complex_cholesky_decomp).

- int `gsl_linalg_cholesky_svx`(const [gsl_matrix](https://www.gnu.org/software/gsl/doc/html/vectors.html#c.gsl_matrix) * *cholesky*, [gsl_vector](https://www.gnu.org/software/gsl/doc/html/vectors.html#c.gsl_vector) * *x*)

- int `gsl_linalg_complex_cholesky_svx`(const gsl_matrix_complex * *cholesky*, gsl_vector_complex * *x*)

  These functions solve the system ![A x = b](https://www.gnu.org/software/gsl/doc/html/_images/math/bf1b0fe2f3273671dbd86939aa2ca8992c570e65.png) in-place using the Cholesky decomposition of ![A](https://www.gnu.org/software/gsl/doc/html/_images/math/3daf97c82cdab413e28d010c06770f224d64dc5a.png)held in the matrix `cholesky` which must have been previously computed by[`gsl_linalg_cholesky_decomp()`](https://www.gnu.org/software/gsl/doc/html/linalg.html#c.gsl_linalg_cholesky_decomp) or [`gsl_linalg_complex_cholesky_decomp()`](https://www.gnu.org/software/gsl/doc/html/linalg.html#c.gsl_linalg_complex_cholesky_decomp). On input `x` should contain the right-hand side ![b](https://www.gnu.org/software/gsl/doc/html/_images/math/88082f6488798774b9cd95a59c1ab4c6cafe909f.png), which is replaced by the solution on output.

- int `gsl_linalg_cholesky_invert`([gsl_matrix](https://www.gnu.org/software/gsl/doc/html/vectors.html#c.gsl_matrix) * *cholesky*)

- int `gsl_linalg_complex_cholesky_invert`(gsl_matrix_complex * *cholesky*)

  These functions compute the inverse of a matrix from its Cholesky decomposition `cholesky`, which must have been previously computed by [`gsl_linalg_cholesky_decomp()`](https://www.gnu.org/software/gsl/doc/html/linalg.html#c.gsl_linalg_cholesky_decomp) or[`gsl_linalg_complex_cholesky_decomp()`](https://www.gnu.org/software/gsl/doc/html/linalg.html#c.gsl_linalg_complex_cholesky_decomp). On output, the inverse is stored in-place in `cholesky`.

- int `gsl_linalg_cholesky_decomp2`([gsl_matrix](https://www.gnu.org/software/gsl/doc/html/vectors.html#c.gsl_matrix) * *A*, [gsl_vector](https://www.gnu.org/software/gsl/doc/html/vectors.html#c.gsl_vector) * *S*)

  This function calculates a diagonal scaling transformation ![S](https://www.gnu.org/software/gsl/doc/html/_images/math/4687d5918c214ce691dee3322ecf3d2046b28e1f.png) for the symmetric, positive-definite square matrix `A`, and then computes the Cholesky decomposition ![S A S = L L^T](https://www.gnu.org/software/gsl/doc/html/_images/math/4ec66f1d0357472b0c97cc4a5313f214dfb9a507.png). On input, the values from the diagonal and lower-triangular part of the matrix `A` are used (the upper triangular part is ignored). On output the diagonal and lower triangular part of the input matrix `A` contain the matrix ![L](https://www.gnu.org/software/gsl/doc/html/_images/math/82359a7c56cc02d70a5ed0e98df631577f7d172b.png), while the upper triangular part of the input matrix is overwritten with ![L^T](https://www.gnu.org/software/gsl/doc/html/_images/math/36fd14d4f6a032c0d439613fa95d6e11c291bf09.png) (the diagonal terms being identical for both ![L](https://www.gnu.org/software/gsl/doc/html/_images/math/82359a7c56cc02d70a5ed0e98df631577f7d172b.png) and ![L^T](https://www.gnu.org/software/gsl/doc/html/_images/math/36fd14d4f6a032c0d439613fa95d6e11c291bf09.png)). If the matrix is not positive-definite then the decomposition will fail, returning the error code [`GSL_EDOM`](https://www.gnu.org/software/gsl/doc/html/err.html#c.GSL_EDOM). The diagonal scale factors are stored in `S` on output.When testing whether a matrix is positive-definite, disable the error handler first to avoid triggering an error.

- int `gsl_linalg_cholesky_solve2`(const [gsl_matrix](https://www.gnu.org/software/gsl/doc/html/vectors.html#c.gsl_matrix) * *cholesky*, const [gsl_vector](https://www.gnu.org/software/gsl/doc/html/vectors.html#c.gsl_vector) * *S*, const [gsl_vector](https://www.gnu.org/software/gsl/doc/html/vectors.html#c.gsl_vector) * *b*, [gsl_vector](https://www.gnu.org/software/gsl/doc/html/vectors.html#c.gsl_vector) * *x*)

  This function solves the system ![(S A S) (S^{-1} x) = S b](https://www.gnu.org/software/gsl/doc/html/_images/math/0f3e2ccfdb50053406b8af953f5aac9f5e3dfa68.png) using the Cholesky decomposition of ![S A S](https://www.gnu.org/software/gsl/doc/html/_images/math/c6b6a6dbdca6ac2a12b75bad4fbd13f57720c8d6.png) held in the matrix `cholesky` which must have been previously computed by [`gsl_linalg_cholesky_decomp2()`](https://www.gnu.org/software/gsl/doc/html/linalg.html#c.gsl_linalg_cholesky_decomp2).

- int `gsl_linalg_cholesky_svx2`(const [gsl_matrix](https://www.gnu.org/software/gsl/doc/html/vectors.html#c.gsl_matrix) * *cholesky*, const [gsl_vector](https://www.gnu.org/software/gsl/doc/html/vectors.html#c.gsl_vector) * *S*, [gsl_vector](https://www.gnu.org/software/gsl/doc/html/vectors.html#c.gsl_vector) * *x*)

  This function solves the system ![(S A S) (S^{-1} x) = S b](https://www.gnu.org/software/gsl/doc/html/_images/math/0f3e2ccfdb50053406b8af953f5aac9f5e3dfa68.png) in-place using the Cholesky decomposition of ![S A S](https://www.gnu.org/software/gsl/doc/html/_images/math/c6b6a6dbdca6ac2a12b75bad4fbd13f57720c8d6.png) held in the matrix `cholesky` which must have been previously computed by [`gsl_linalg_cholesky_decomp2()`](https://www.gnu.org/software/gsl/doc/html/linalg.html#c.gsl_linalg_cholesky_decomp2). On input `x` should contain the right-hand side ![b](https://www.gnu.org/software/gsl/doc/html/_images/math/88082f6488798774b9cd95a59c1ab4c6cafe909f.png), which is replaced by the solution on output.

- int `gsl_linalg_cholesky_scale`(const [gsl_matrix](https://www.gnu.org/software/gsl/doc/html/vectors.html#c.gsl_matrix) * *A*, [gsl_vector](https://www.gnu.org/software/gsl/doc/html/vectors.html#c.gsl_vector) * *S*)

  This function calculates a diagonal scaling transformation of the symmetric, positive definite matrix `A`, such that ![S A S](https://www.gnu.org/software/gsl/doc/html/_images/math/c6b6a6dbdca6ac2a12b75bad4fbd13f57720c8d6.png) has a condition number within a factor of ![N](https://www.gnu.org/software/gsl/doc/html/_images/math/99f8a3c97c19212ca3eee7d00dddc8d68b2d1c27.png) of the matrix of smallest possible condition number over all possible diagonal scalings. On output, `S` contains the scale factors, given by ![S_i = 1/\sqrt{A_{ii}}](https://www.gnu.org/software/gsl/doc/html/_images/math/23f0805db32e2cd7e6024cb71a04a546ca94a15f.png). For any ![A_{ii} \le 0](https://www.gnu.org/software/gsl/doc/html/_images/math/bd7f6fa1bf0e5d9d759ec18c7cd09fc254e5f28d.png), the corresponding scale factor ![S_i](https://www.gnu.org/software/gsl/doc/html/_images/math/6aeb172cc40a350fa0bff6b6d2ea84baaf0bfad5.png) is set to ![1](https://www.gnu.org/software/gsl/doc/html/_images/math/585581da152cce359927e7b9eeadef5ca7c3a42d.png).

- int `gsl_linalg_cholesky_scale_apply`([gsl_matrix](https://www.gnu.org/software/gsl/doc/html/vectors.html#c.gsl_matrix) * *A*, const [gsl_vector](https://www.gnu.org/software/gsl/doc/html/vectors.html#c.gsl_vector) * *S*)

  This function applies the scaling transformation `S` to the matrix `A`. On output, `A` is replaced by ![S A S](https://www.gnu.org/software/gsl/doc/html/_images/math/c6b6a6dbdca6ac2a12b75bad4fbd13f57720c8d6.png).

- int `gsl_linalg_cholesky_rcond`(const [gsl_matrix](https://www.gnu.org/software/gsl/doc/html/vectors.html#c.gsl_matrix) * *cholesky*, double * *rcond*, [gsl_vector](https://www.gnu.org/software/gsl/doc/html/vectors.html#c.gsl_vector) * *work*)

  This function estimates the reciprocal condition number (using the 1-norm) of the symmetric positive definite matrix ![A](https://www.gnu.org/software/gsl/doc/html/_images/math/3daf97c82cdab413e28d010c06770f224d64dc5a.png), using its Cholesky decomposition provided in `cholesky`. The reciprocal condition number estimate, defined as ![1 / (||A||_1 \cdot ||A^{-1}||_1)](https://www.gnu.org/software/gsl/doc/html/_images/math/ad2d6cfcbfb6f15efb4fc25fc4063b7b7fe87d1b.png), is stored in `rcond`. Additional workspace of size ![3 N](https://www.gnu.org/software/gsl/doc/html/_images/math/1b2b0ffdf792c9f396c727992622338aadd67f74.png) is required in `work`.



## Pivoted Cholesky Decomposition

A symmetric, positive definite square matrix ![A](https://www.gnu.org/software/gsl/doc/html/_images/math/3daf97c82cdab413e28d010c06770f224d64dc5a.png) has an alternate Cholesky decomposition into a product of a lower unit triangular matrix ![L](https://www.gnu.org/software/gsl/doc/html/_images/math/82359a7c56cc02d70a5ed0e98df631577f7d172b.png), a diagonal matrix ![D](https://www.gnu.org/software/gsl/doc/html/_images/math/0638e3c183dd8cfd1790c25e762f1a8c53ee5711.png) and ![L^T](https://www.gnu.org/software/gsl/doc/html/_images/math/36fd14d4f6a032c0d439613fa95d6e11c291bf09.png), given by ![L D L^T](https://www.gnu.org/software/gsl/doc/html/_images/math/3094d6ec30440e607d137156778e58e40fa1758f.png). This is equivalent to the Cholesky formulation discussed above, with the standard Cholesky lower triangular factor given by ![L D^{1 \over 2}](https://www.gnu.org/software/gsl/doc/html/_images/math/1e28ba5779262fef8c0268a70e538e0c335260d3.png). For ill-conditioned matrices, it can help to use a pivoting strategy to prevent the entries of ![D](https://www.gnu.org/software/gsl/doc/html/_images/math/0638e3c183dd8cfd1790c25e762f1a8c53ee5711.png) and ![L](https://www.gnu.org/software/gsl/doc/html/_images/math/82359a7c56cc02d70a5ed0e98df631577f7d172b.png) from growing too large, and also ensure ![D_1 \ge D_2 \ge \cdots \ge D_n > 0](https://www.gnu.org/software/gsl/doc/html/_images/math/47a9389bf30c7672dec34d927d110da32ee82d9c.png), where ![D_i](https://www.gnu.org/software/gsl/doc/html/_images/math/86a7c7eebe7ae8c4122594e145236e3a74d67ed9.png) are the diagonal entries of ![D](https://www.gnu.org/software/gsl/doc/html/_images/math/0638e3c183dd8cfd1790c25e762f1a8c53ee5711.png). The final decomposition is given by

![P A P^T = L D L^T](https://www.gnu.org/software/gsl/doc/html/_images/math/d4617f5028854027c9907bc6f12adbf8f4263cf7.png)

where ![P](https://www.gnu.org/software/gsl/doc/html/_images/math/d1a0a8c0064fe114c19c8fd71ffec7a98480f707.png) is a permutation matrix.

- int `gsl_linalg_pcholesky_decomp`([gsl_matrix](https://www.gnu.org/software/gsl/doc/html/vectors.html#c.gsl_matrix) * *A*, [gsl_permutation](https://www.gnu.org/software/gsl/doc/html/permutation.html#c.gsl_permutation) * *p*)

  This function factors the symmetric, positive-definite square matrix `A` into the Pivoted Cholesky decomposition ![P A P^T = L D L^T](https://www.gnu.org/software/gsl/doc/html/_images/math/320707d7dfe941864b1179369c884e959b14795a.png). On input, the values from the diagonal and lower-triangular part of the matrix `A` are used to construct the factorization. On output the diagonal of the input matrix `A` stores the diagonal elements of ![D](https://www.gnu.org/software/gsl/doc/html/_images/math/0638e3c183dd8cfd1790c25e762f1a8c53ee5711.png), and the lower triangular portion of `A`contains the matrix ![L](https://www.gnu.org/software/gsl/doc/html/_images/math/82359a7c56cc02d70a5ed0e98df631577f7d172b.png). Since ![L](https://www.gnu.org/software/gsl/doc/html/_images/math/82359a7c56cc02d70a5ed0e98df631577f7d172b.png) has ones on its diagonal these do not need to be explicitely stored. The upper triangular portion of `A` is unmodified. The permutation matrix ![P](https://www.gnu.org/software/gsl/doc/html/_images/math/d1a0a8c0064fe114c19c8fd71ffec7a98480f707.png) is stored in `p` on output.

- int `gsl_linalg_pcholesky_solve`(const [gsl_matrix](https://www.gnu.org/software/gsl/doc/html/vectors.html#c.gsl_matrix) * *LDLT*, const [gsl_permutation](https://www.gnu.org/software/gsl/doc/html/permutation.html#c.gsl_permutation) * *p*, const [gsl_vector](https://www.gnu.org/software/gsl/doc/html/vectors.html#c.gsl_vector) * *b*, [gsl_vector](https://www.gnu.org/software/gsl/doc/html/vectors.html#c.gsl_vector) * *x*)

  This function solves the system ![A x = b](https://www.gnu.org/software/gsl/doc/html/_images/math/bf1b0fe2f3273671dbd86939aa2ca8992c570e65.png) using the Pivoted Cholesky decomposition of ![A](https://www.gnu.org/software/gsl/doc/html/_images/math/3daf97c82cdab413e28d010c06770f224d64dc5a.png) held in the matrix `LDLT` and permutation `p` which must have been previously computed by [`gsl_linalg_pcholesky_decomp()`](https://www.gnu.org/software/gsl/doc/html/linalg.html#c.gsl_linalg_pcholesky_decomp).

- int `gsl_linalg_pcholesky_svx`(const [gsl_matrix](https://www.gnu.org/software/gsl/doc/html/vectors.html#c.gsl_matrix) * *LDLT*, const [gsl_permutation](https://www.gnu.org/software/gsl/doc/html/permutation.html#c.gsl_permutation) * *p*, [gsl_vector](https://www.gnu.org/software/gsl/doc/html/vectors.html#c.gsl_vector)* *x*)

  This function solves the system ![A x = b](https://www.gnu.org/software/gsl/doc/html/_images/math/bf1b0fe2f3273671dbd86939aa2ca8992c570e65.png) in-place using the Pivoted Cholesky decomposition of ![A](https://www.gnu.org/software/gsl/doc/html/_images/math/3daf97c82cdab413e28d010c06770f224d64dc5a.png) held in the matrix `LDLT` and permutation `p` which must have been previously computed by [`gsl_linalg_pcholesky_decomp()`](https://www.gnu.org/software/gsl/doc/html/linalg.html#c.gsl_linalg_pcholesky_decomp). On input, `x` contains the right hand side vector ![b](https://www.gnu.org/software/gsl/doc/html/_images/math/88082f6488798774b9cd95a59c1ab4c6cafe909f.png) which is replaced by the solution vector on output.

- int `gsl_linalg_pcholesky_decomp2`([gsl_matrix](https://www.gnu.org/software/gsl/doc/html/vectors.html#c.gsl_matrix) * *A*, [gsl_permutation](https://www.gnu.org/software/gsl/doc/html/permutation.html#c.gsl_permutation) * *p*, [gsl_vector](https://www.gnu.org/software/gsl/doc/html/vectors.html#c.gsl_vector) * *S*)

  This function computes the pivoted Cholesky factorization of the matrix ![S A S](https://www.gnu.org/software/gsl/doc/html/_images/math/c6b6a6dbdca6ac2a12b75bad4fbd13f57720c8d6.png), where the input matrix `A` is symmetric and positive definite, and the diagonal scaling matrix `S` is computed to reduce the condition number of `A` as much as possible. See [Cholesky Decomposition](https://www.gnu.org/software/gsl/doc/html/linalg.html#sec-cholesky-decomposition) for more information on the matrix `S`. The Pivoted Cholesky decomposition satisfies ![P S A S P^T = L D L^T](https://www.gnu.org/software/gsl/doc/html/_images/math/85004b30804e4e5945d274f68588079f276d636e.png). On input, the values from the diagonal and lower-triangular part of the matrix `A` are used to construct the factorization. On output the diagonal of the input matrix `A`stores the diagonal elements of ![D](https://www.gnu.org/software/gsl/doc/html/_images/math/0638e3c183dd8cfd1790c25e762f1a8c53ee5711.png), and the lower triangular portion of `A` contains the matrix ![L](https://www.gnu.org/software/gsl/doc/html/_images/math/82359a7c56cc02d70a5ed0e98df631577f7d172b.png). Since ![L](https://www.gnu.org/software/gsl/doc/html/_images/math/82359a7c56cc02d70a5ed0e98df631577f7d172b.png) has ones on its diagonal these do not need to be explicitely stored. The upper triangular portion of `A` is unmodified. The permutation matrix ![P](https://www.gnu.org/software/gsl/doc/html/_images/math/d1a0a8c0064fe114c19c8fd71ffec7a98480f707.png) is stored in `p` on output. The diagonal scaling transformation is stored in `S` on output.

- int `gsl_linalg_pcholesky_solve2`(const [gsl_matrix](https://www.gnu.org/software/gsl/doc/html/vectors.html#c.gsl_matrix) * *LDLT*, const [gsl_permutation](https://www.gnu.org/software/gsl/doc/html/permutation.html#c.gsl_permutation) * *p*, const [gsl_vector](https://www.gnu.org/software/gsl/doc/html/vectors.html#c.gsl_vector) * *S*, const [gsl_vector](https://www.gnu.org/software/gsl/doc/html/vectors.html#c.gsl_vector) * *b*, [gsl_vector](https://www.gnu.org/software/gsl/doc/html/vectors.html#c.gsl_vector) * *x*)

  This function solves the system ![(S A S) (S^{-1} x) = S b](https://www.gnu.org/software/gsl/doc/html/_images/math/0f3e2ccfdb50053406b8af953f5aac9f5e3dfa68.png) using the Pivoted Cholesky decomposition of ![S A S](https://www.gnu.org/software/gsl/doc/html/_images/math/c6b6a6dbdca6ac2a12b75bad4fbd13f57720c8d6.png) held in the matrix `LDLT`, permutation `p`, and vector `S`, which must have been previously computed by [`gsl_linalg_pcholesky_decomp2()`](https://www.gnu.org/software/gsl/doc/html/linalg.html#c.gsl_linalg_pcholesky_decomp2).

- int `gsl_linalg_pcholesky_svx2`(const [gsl_matrix](https://www.gnu.org/software/gsl/doc/html/vectors.html#c.gsl_matrix) * *LDLT*, const [gsl_permutation](https://www.gnu.org/software/gsl/doc/html/permutation.html#c.gsl_permutation) * *p*, const [gsl_vector](https://www.gnu.org/software/gsl/doc/html/vectors.html#c.gsl_vector) * *S*, [gsl_vector](https://www.gnu.org/software/gsl/doc/html/vectors.html#c.gsl_vector) * *x*)

  This function solves the system ![(S A S) (S^{-1} x) = S b](https://www.gnu.org/software/gsl/doc/html/_images/math/0f3e2ccfdb50053406b8af953f5aac9f5e3dfa68.png) in-place using the Pivoted Cholesky decomposition of ![S A S](https://www.gnu.org/software/gsl/doc/html/_images/math/c6b6a6dbdca6ac2a12b75bad4fbd13f57720c8d6.png) held in the matrix `LDLT`, permutation `p` and vector `S`, which must have been previously computed by [`gsl_linalg_pcholesky_decomp2()`](https://www.gnu.org/software/gsl/doc/html/linalg.html#c.gsl_linalg_pcholesky_decomp2). On input, `x` contains the right hand side vector ![b](https://www.gnu.org/software/gsl/doc/html/_images/math/88082f6488798774b9cd95a59c1ab4c6cafe909f.png) which is replaced by the solution vector on output.

- int `gsl_linalg_pcholesky_invert`(const [gsl_matrix](https://www.gnu.org/software/gsl/doc/html/vectors.html#c.gsl_matrix) * *LDLT*, const [gsl_permutation](https://www.gnu.org/software/gsl/doc/html/permutation.html#c.gsl_permutation) * *p*, [gsl_matrix](https://www.gnu.org/software/gsl/doc/html/vectors.html#c.gsl_matrix) * *Ainv*)

  This function computes the inverse of the matrix ![A](https://www.gnu.org/software/gsl/doc/html/_images/math/3daf97c82cdab413e28d010c06770f224d64dc5a.png), using the Pivoted Cholesky decomposition stored in `LDLT` and `p`. On output, the matrix `Ainv` contains ![A^{-1}](https://www.gnu.org/software/gsl/doc/html/_images/math/a6da7c3a6effc3b823cc440d3b450a82b5422f0d.png).

- int `gsl_linalg_pcholesky_rcond`(const [gsl_matrix](https://www.gnu.org/software/gsl/doc/html/vectors.html#c.gsl_matrix) * *LDLT*, const [gsl_permutation](https://www.gnu.org/software/gsl/doc/html/permutation.html#c.gsl_permutation) * *p*, double * *rcond*, [gsl_vector](https://www.gnu.org/software/gsl/doc/html/vectors.html#c.gsl_vector) * *work*)

  This function estimates the reciprocal condition number (using the 1-norm) of the symmetric positive definite matrix ![A](https://www.gnu.org/software/gsl/doc/html/_images/math/3daf97c82cdab413e28d010c06770f224d64dc5a.png), using its pivoted Cholesky decomposition provided in `LDLT`. The reciprocal condition number estimate, defined as ![1 / (||A||_1 \cdot ||A^{-1}||_1)](https://www.gnu.org/software/gsl/doc/html/_images/math/ad2d6cfcbfb6f15efb4fc25fc4063b7b7fe87d1b.png), is stored in `rcond`. Additional workspace of size ![3 N](https://www.gnu.org/software/gsl/doc/html/_images/math/1b2b0ffdf792c9f396c727992622338aadd67f74.png) is required in `work`.



## Modified Cholesky Decomposition

The modified Cholesky decomposition is suitable for solving systems ![A x = b](https://www.gnu.org/software/gsl/doc/html/_images/math/bf1b0fe2f3273671dbd86939aa2ca8992c570e65.png) where ![A](https://www.gnu.org/software/gsl/doc/html/_images/math/3daf97c82cdab413e28d010c06770f224d64dc5a.png) is a symmetric indefinite matrix. Such matrices arise in nonlinear optimization algorithms. The standard Cholesky decomposition requires a positive definite matrix and would fail in this case. Instead of resorting to a method like QR or SVD, which do not take into account the symmetry of the matrix, we can instead introduce a small perturbation to the matrix ![A](https://www.gnu.org/software/gsl/doc/html/_images/math/3daf97c82cdab413e28d010c06770f224d64dc5a.png) to make it positive definite, and then use a Cholesky decomposition on the perturbed matrix. The resulting decomposition satisfies

![P (A + E) P^T = L D L^T](https://www.gnu.org/software/gsl/doc/html/_images/math/f2ee7c3e93b4520b02237da70a933b6a91c108f9.png)

where ![P](https://www.gnu.org/software/gsl/doc/html/_images/math/d1a0a8c0064fe114c19c8fd71ffec7a98480f707.png) is a permutation matrix, ![E](https://www.gnu.org/software/gsl/doc/html/_images/math/03d5ff4c4bb576e884a9d0d9b69ca58d7117a184.png) is a diagonal perturbation matrix, ![L](https://www.gnu.org/software/gsl/doc/html/_images/math/82359a7c56cc02d70a5ed0e98df631577f7d172b.png) is unit lower triangular, and ![D](https://www.gnu.org/software/gsl/doc/html/_images/math/0638e3c183dd8cfd1790c25e762f1a8c53ee5711.png) is diagonal. If ![A](https://www.gnu.org/software/gsl/doc/html/_images/math/3daf97c82cdab413e28d010c06770f224d64dc5a.png) is sufficiently positive definite, then the perturbation matrix ![E](https://www.gnu.org/software/gsl/doc/html/_images/math/03d5ff4c4bb576e884a9d0d9b69ca58d7117a184.png) will be zero and this method is equivalent to the pivoted Cholesky algorithm. For indefinite matrices, the perturbation matrix ![E](https://www.gnu.org/software/gsl/doc/html/_images/math/03d5ff4c4bb576e884a9d0d9b69ca58d7117a184.png) is computed to ensure that ![A + E](https://www.gnu.org/software/gsl/doc/html/_images/math/47846a64dfecb9a65510de759cc683e6b287160a.png) is positive definite and well conditioned.

- int `gsl_linalg_mcholesky_decomp`([gsl_matrix](https://www.gnu.org/software/gsl/doc/html/vectors.html#c.gsl_matrix) * *A*, [gsl_permutation](https://www.gnu.org/software/gsl/doc/html/permutation.html#c.gsl_permutation) * *p*, [gsl_vector](https://www.gnu.org/software/gsl/doc/html/vectors.html#c.gsl_vector) * *E*)

  This function factors the symmetric, indefinite square matrix `A` into the Modified Cholesky decomposition ![P (A + E) P^T = L D L^T](https://www.gnu.org/software/gsl/doc/html/_images/math/e127a82487043c7ad0e68f842dcc1b2b3eb94e0b.png). On input, the values from the diagonal and lower-triangular part of the matrix `A` are used to construct the factorization. On output the diagonal of the input matrix `A` stores the diagonal elements of ![D](https://www.gnu.org/software/gsl/doc/html/_images/math/0638e3c183dd8cfd1790c25e762f1a8c53ee5711.png), and the lower triangular portion of `A`contains the matrix ![L](https://www.gnu.org/software/gsl/doc/html/_images/math/82359a7c56cc02d70a5ed0e98df631577f7d172b.png). Since ![L](https://www.gnu.org/software/gsl/doc/html/_images/math/82359a7c56cc02d70a5ed0e98df631577f7d172b.png) has ones on its diagonal these do not need to be explicitely stored. The upper triangular portion of `A` is unmodified. The permutation matrix ![P](https://www.gnu.org/software/gsl/doc/html/_images/math/d1a0a8c0064fe114c19c8fd71ffec7a98480f707.png) is stored in `p` on output. The diagonal perturbation matrix is stored in `E` on output. The parameter `E`may be set to NULL if it is not required.

- int `gsl_linalg_mcholesky_solve`(const [gsl_matrix](https://www.gnu.org/software/gsl/doc/html/vectors.html#c.gsl_matrix) * *LDLT*, const [gsl_permutation](https://www.gnu.org/software/gsl/doc/html/permutation.html#c.gsl_permutation) * *p*, const [gsl_vector](https://www.gnu.org/software/gsl/doc/html/vectors.html#c.gsl_vector) * *b*, [gsl_vector](https://www.gnu.org/software/gsl/doc/html/vectors.html#c.gsl_vector) * *x*)

  This function solves the perturbed system ![(A + E) x = b](https://www.gnu.org/software/gsl/doc/html/_images/math/0bd1bf7ba1926d76c0a499ffb17c7aefe7ff8710.png) using the Cholesky decomposition of ![A + E](https://www.gnu.org/software/gsl/doc/html/_images/math/47846a64dfecb9a65510de759cc683e6b287160a.png) held in the matrix `LDLT` and permutation `p` which must have been previously computed by [`gsl_linalg_mcholesky_decomp()`](https://www.gnu.org/software/gsl/doc/html/linalg.html#c.gsl_linalg_mcholesky_decomp).

- int `gsl_linalg_mcholesky_svx`(const [gsl_matrix](https://www.gnu.org/software/gsl/doc/html/vectors.html#c.gsl_matrix) * *LDLT*, const [gsl_permutation](https://www.gnu.org/software/gsl/doc/html/permutation.html#c.gsl_permutation) * *p*, [gsl_vector](https://www.gnu.org/software/gsl/doc/html/vectors.html#c.gsl_vector)* *x*)

  This function solves the perturbed system ![(A + E) x = b](https://www.gnu.org/software/gsl/doc/html/_images/math/0bd1bf7ba1926d76c0a499ffb17c7aefe7ff8710.png) in-place using the Cholesky decomposition of ![A + E](https://www.gnu.org/software/gsl/doc/html/_images/math/47846a64dfecb9a65510de759cc683e6b287160a.png) held in the matrix `LDLT` and permutation `p` which must have been previously computed by [`gsl_linalg_mcholesky_decomp()`](https://www.gnu.org/software/gsl/doc/html/linalg.html#c.gsl_linalg_mcholesky_decomp). On input, `x` contains the right hand side vector ![b](https://www.gnu.org/software/gsl/doc/html/_images/math/88082f6488798774b9cd95a59c1ab4c6cafe909f.png) which is replaced by the solution vector on output.

- int `gsl_linalg_mcholesky_rcond`(const [gsl_matrix](https://www.gnu.org/software/gsl/doc/html/vectors.html#c.gsl_matrix) * *LDLT*, const [gsl_permutation](https://www.gnu.org/software/gsl/doc/html/permutation.html#c.gsl_permutation) * *p*, double * *rcond*, [gsl_vector](https://www.gnu.org/software/gsl/doc/html/vectors.html#c.gsl_vector) * *work*)

  This function estimates the reciprocal condition number (using the 1-norm) of the perturbed matrix ![A + E](https://www.gnu.org/software/gsl/doc/html/_images/math/47846a64dfecb9a65510de759cc683e6b287160a.png), using its pivoted Cholesky decomposition provided in `LDLT`. The reciprocal condition number estimate, defined as ![1 / (||A + E||_1 \cdot ||(A + E)^{-1}||_1)](https://www.gnu.org/software/gsl/doc/html/_images/math/647987804cbbc05d094e05bf8892d2e62f18d2e7.png), is stored in `rcond`. Additional workspace of size ![3 N](https://www.gnu.org/software/gsl/doc/html/_images/math/1b2b0ffdf792c9f396c727992622338aadd67f74.png) is required in `work`.



## Tridiagonal Decomposition of Real Symmetric Matrices

A symmetric matrix ![A](https://www.gnu.org/software/gsl/doc/html/_images/math/3daf97c82cdab413e28d010c06770f224d64dc5a.png) can be factorized by similarity transformations into the form,

![A = Q T Q^T](https://www.gnu.org/software/gsl/doc/html/_images/math/385335428dd4a63b4a6dfa14833a985ab14971b3.png)

where ![Q](https://www.gnu.org/software/gsl/doc/html/_images/math/7fcf2b077fb63ef0b6eb5c350673266b94327847.png) is an orthogonal matrix and ![T](https://www.gnu.org/software/gsl/doc/html/_images/math/0dcfe86f94dfb3243931655ff3b097dd326e229f.png) is a symmetric tridiagonal matrix.

- int `gsl_linalg_symmtd_decomp`([gsl_matrix](https://www.gnu.org/software/gsl/doc/html/vectors.html#c.gsl_matrix) * *A*, [gsl_vector](https://www.gnu.org/software/gsl/doc/html/vectors.html#c.gsl_vector) * *tau*)

  This function factorizes the symmetric square matrix `A` into the symmetric tridiagonal decomposition ![Q T Q^T](https://www.gnu.org/software/gsl/doc/html/_images/math/00bee92f5dade5cc840b201629063fa99664bf10.png). On output the diagonal and subdiagonal part of the input matrix `A`contain the tridiagonal matrix ![T](https://www.gnu.org/software/gsl/doc/html/_images/math/0dcfe86f94dfb3243931655ff3b097dd326e229f.png). The remaining lower triangular part of the input matrix contains the Householder vectors which, together with the Householder coefficients `tau`, encode the orthogonal matrix ![Q](https://www.gnu.org/software/gsl/doc/html/_images/math/7fcf2b077fb63ef0b6eb5c350673266b94327847.png). This storage scheme is the same as used by LAPACK. The upper triangular part of `A` is not referenced.

- int `gsl_linalg_symmtd_unpack`(const [gsl_matrix](https://www.gnu.org/software/gsl/doc/html/vectors.html#c.gsl_matrix) * *A*, const [gsl_vector](https://www.gnu.org/software/gsl/doc/html/vectors.html#c.gsl_vector) * *tau*, [gsl_matrix](https://www.gnu.org/software/gsl/doc/html/vectors.html#c.gsl_matrix) * *Q*, [gsl_vector](https://www.gnu.org/software/gsl/doc/html/vectors.html#c.gsl_vector) * *diag*, [gsl_vector](https://www.gnu.org/software/gsl/doc/html/vectors.html#c.gsl_vector) * *subdiag*)

  This function unpacks the encoded symmetric tridiagonal decomposition (`A`, `tau`) obtained from [`gsl_linalg_symmtd_decomp()`](https://www.gnu.org/software/gsl/doc/html/linalg.html#c.gsl_linalg_symmtd_decomp) into the orthogonal matrix `Q`, the vector of diagonal elements `diag` and the vector of subdiagonal elements `subdiag`.

- int `gsl_linalg_symmtd_unpack_T`(const [gsl_matrix](https://www.gnu.org/software/gsl/doc/html/vectors.html#c.gsl_matrix) * *A*, [gsl_vector](https://www.gnu.org/software/gsl/doc/html/vectors.html#c.gsl_vector) * *diag*, [gsl_vector](https://www.gnu.org/software/gsl/doc/html/vectors.html#c.gsl_vector) * *subdiag*)

  This function unpacks the diagonal and subdiagonal of the encoded symmetric tridiagonal decomposition (`A`, `tau`) obtained from [`gsl_linalg_symmtd_decomp()`](https://www.gnu.org/software/gsl/doc/html/linalg.html#c.gsl_linalg_symmtd_decomp) into the vectors `diag`and `subdiag`.



## Tridiagonal Decomposition of Hermitian Matrices

A hermitian matrix ![A](https://www.gnu.org/software/gsl/doc/html/_images/math/3daf97c82cdab413e28d010c06770f224d64dc5a.png) can be factorized by similarity transformations into the form,

![A = U T U^T](https://www.gnu.org/software/gsl/doc/html/_images/math/41e4aaa86690a7afbdc69ef8eec98938fba3459b.png)

where ![U](https://www.gnu.org/software/gsl/doc/html/_images/math/3d0e8e66abfde1273e99e86e8a090a7073d842e3.png) is a unitary matrix and ![T](https://www.gnu.org/software/gsl/doc/html/_images/math/0dcfe86f94dfb3243931655ff3b097dd326e229f.png) is a real symmetric tridiagonal matrix.

- int `gsl_linalg_hermtd_decomp`(gsl_matrix_complex * *A*, gsl_vector_complex * *tau*)

  This function factorizes the hermitian matrix `A` into the symmetric tridiagonal decomposition ![U T U^T](https://www.gnu.org/software/gsl/doc/html/_images/math/fc3820568b8d8d4d9b70f8f62ca0cb620b7c9a7d.png). On output the real parts of the diagonal and subdiagonal part of the input matrix `A`contain the tridiagonal matrix ![T](https://www.gnu.org/software/gsl/doc/html/_images/math/0dcfe86f94dfb3243931655ff3b097dd326e229f.png). The remaining lower triangular part of the input matrix contains the Householder vectors which, together with the Householder coefficients `tau`, encode the unitary matrix ![U](https://www.gnu.org/software/gsl/doc/html/_images/math/3d0e8e66abfde1273e99e86e8a090a7073d842e3.png). This storage scheme is the same as used by LAPACK. The upper triangular part of `A` and imaginary parts of the diagonal are not referenced.

- int `gsl_linalg_hermtd_unpack`(const gsl_matrix_complex * *A*, const gsl_vector_complex * *tau*, gsl_matrix_complex * *U*, [gsl_vector](https://www.gnu.org/software/gsl/doc/html/vectors.html#c.gsl_vector) * *diag*, [gsl_vector](https://www.gnu.org/software/gsl/doc/html/vectors.html#c.gsl_vector) * *subdiag*)

  This function unpacks the encoded tridiagonal decomposition (`A`, `tau`) obtained from [`gsl_linalg_hermtd_decomp()`](https://www.gnu.org/software/gsl/doc/html/linalg.html#c.gsl_linalg_hermtd_decomp) into the unitary matrix `U`, the real vector of diagonal elements `diag` and the real vector of subdiagonal elements `subdiag`.

- int `gsl_linalg_hermtd_unpack_T`(const gsl_matrix_complex * *A*, [gsl_vector](https://www.gnu.org/software/gsl/doc/html/vectors.html#c.gsl_vector) * *diag*, [gsl_vector](https://www.gnu.org/software/gsl/doc/html/vectors.html#c.gsl_vector)* *subdiag*)

  This function unpacks the diagonal and subdiagonal of the encoded tridiagonal decomposition (`A`, `tau`) obtained from the [`gsl_linalg_hermtd_decomp()`](https://www.gnu.org/software/gsl/doc/html/linalg.html#c.gsl_linalg_hermtd_decomp) into the real vectors `diag` and `subdiag`.



## Hessenberg Decomposition of Real Matrices

A general real matrix ![A](https://www.gnu.org/software/gsl/doc/html/_images/math/3daf97c82cdab413e28d010c06770f224d64dc5a.png) can be decomposed by orthogonal similarity transformations into the form

![A = U H U^T](https://www.gnu.org/software/gsl/doc/html/_images/math/499f079d54111c2e0ecfb66248cc29bee766577c.png)

where ![U](https://www.gnu.org/software/gsl/doc/html/_images/math/3d0e8e66abfde1273e99e86e8a090a7073d842e3.png) is orthogonal and ![H](https://www.gnu.org/software/gsl/doc/html/_images/math/af029cc216351e86964a2146cb3a293d9336ddd0.png) is an upper Hessenberg matrix, meaning that it has zeros below the first subdiagonal. The Hessenberg reduction is the first step in the Schur decomposition for the nonsymmetric eigenvalue problem, but has applications in other areas as well.

- int `gsl_linalg_hessenberg_decomp`([gsl_matrix](https://www.gnu.org/software/gsl/doc/html/vectors.html#c.gsl_matrix) * *A*, [gsl_vector](https://www.gnu.org/software/gsl/doc/html/vectors.html#c.gsl_vector) * *tau*)

  This function computes the Hessenberg decomposition of the matrix `A` by applying the similarity transformation ![H = U^T A U](https://www.gnu.org/software/gsl/doc/html/_images/math/a1dc7a87c85e25a78a5f59b0152391e7bad2e252.png). On output, ![H](https://www.gnu.org/software/gsl/doc/html/_images/math/af029cc216351e86964a2146cb3a293d9336ddd0.png) is stored in the upper portion of `A`. The information required to construct the matrix ![U](https://www.gnu.org/software/gsl/doc/html/_images/math/3d0e8e66abfde1273e99e86e8a090a7073d842e3.png) is stored in the lower triangular portion of `A`. ![U](https://www.gnu.org/software/gsl/doc/html/_images/math/3d0e8e66abfde1273e99e86e8a090a7073d842e3.png)is a product of ![N - 2](https://www.gnu.org/software/gsl/doc/html/_images/math/4201c12bc78aa21b88d038a1053f8f8378aa2894.png) Householder matrices. The Householder vectors are stored in the lower portion of `A` (below the subdiagonal) and the Householder coefficients are stored in the vector `tau`. `tau` must be of length `N`.

- int `gsl_linalg_hessenberg_unpack`([gsl_matrix](https://www.gnu.org/software/gsl/doc/html/vectors.html#c.gsl_matrix) * *H*, [gsl_vector](https://www.gnu.org/software/gsl/doc/html/vectors.html#c.gsl_vector) * *tau*, [gsl_matrix](https://www.gnu.org/software/gsl/doc/html/vectors.html#c.gsl_matrix) * *U*)

  This function constructs the orthogonal matrix ![U](https://www.gnu.org/software/gsl/doc/html/_images/math/3d0e8e66abfde1273e99e86e8a090a7073d842e3.png) from the information stored in the Hessenberg matrix `H` along with the vector `tau`. `H` and `tau` are outputs from[`gsl_linalg_hessenberg_decomp()`](https://www.gnu.org/software/gsl/doc/html/linalg.html#c.gsl_linalg_hessenberg_decomp).

- int `gsl_linalg_hessenberg_unpack_accum`([gsl_matrix](https://www.gnu.org/software/gsl/doc/html/vectors.html#c.gsl_matrix) * *H*, [gsl_vector](https://www.gnu.org/software/gsl/doc/html/vectors.html#c.gsl_vector) * *tau*, [gsl_matrix](https://www.gnu.org/software/gsl/doc/html/vectors.html#c.gsl_matrix) * *V*)

  This function is similar to [`gsl_linalg_hessenberg_unpack()`](https://www.gnu.org/software/gsl/doc/html/linalg.html#c.gsl_linalg_hessenberg_unpack), except it accumulates the matrix `U`into `V`, so that ![V' = VU](https://www.gnu.org/software/gsl/doc/html/_images/math/339a53281590258730c58a3b584d5f6d5aa2bbf5.png). The matrix `V` must be initialized prior to calling this function. Setting `V` to the identity matrix provides the same result as [`gsl_linalg_hessenberg_unpack()`](https://www.gnu.org/software/gsl/doc/html/linalg.html#c.gsl_linalg_hessenberg_unpack). If `H` is order `N`, then `V` must have `N` columns but may have any number of rows.

- int `gsl_linalg_hessenberg_set_zero`([gsl_matrix](https://www.gnu.org/software/gsl/doc/html/vectors.html#c.gsl_matrix) * *H*)

  This function sets the lower triangular portion of `H`, below the subdiagonal, to zero. It is useful for clearing out the Householder vectors after calling [`gsl_linalg_hessenberg_decomp()`](https://www.gnu.org/software/gsl/doc/html/linalg.html#c.gsl_linalg_hessenberg_decomp).



## Hessenberg-Triangular Decomposition of Real Matrices

A general real matrix pair (![A](https://www.gnu.org/software/gsl/doc/html/_images/math/3daf97c82cdab413e28d010c06770f224d64dc5a.png), ![B](https://www.gnu.org/software/gsl/doc/html/_images/math/39ac5dd0978410e4f612dd00070a6227cb1ceb10.png)) can be decomposed by orthogonal similarity transformations into the form

![A &= U H V^T \\ B &= U R V^T](https://www.gnu.org/software/gsl/doc/html/_images/math/7a2b57d65fbf63a4707f7c32298f18f7c8d8b8ee.png)

where ![U](https://www.gnu.org/software/gsl/doc/html/_images/math/3d0e8e66abfde1273e99e86e8a090a7073d842e3.png) and ![V](https://www.gnu.org/software/gsl/doc/html/_images/math/632777989bf30753abbbf45a4a9b07ebf0fd36ad.png) are orthogonal, ![H](https://www.gnu.org/software/gsl/doc/html/_images/math/af029cc216351e86964a2146cb3a293d9336ddd0.png) is an upper Hessenberg matrix, and ![R](https://www.gnu.org/software/gsl/doc/html/_images/math/c9b67a086a723a5eec85a146053dd0a2d2378bac.png) is upper triangular. The Hessenberg-Triangular reduction is the first step in the generalized Schur decomposition for the generalized eigenvalue problem.

- int `gsl_linalg_hesstri_decomp`([gsl_matrix](https://www.gnu.org/software/gsl/doc/html/vectors.html#c.gsl_matrix) * *A*, [gsl_matrix](https://www.gnu.org/software/gsl/doc/html/vectors.html#c.gsl_matrix) * *B*, [gsl_matrix](https://www.gnu.org/software/gsl/doc/html/vectors.html#c.gsl_matrix) * *U*, [gsl_matrix](https://www.gnu.org/software/gsl/doc/html/vectors.html#c.gsl_matrix) * *V*, [gsl_vector](https://www.gnu.org/software/gsl/doc/html/vectors.html#c.gsl_vector) * *work*)

  This function computes the Hessenberg-Triangular decomposition of the matrix pair (`A`, `B`). On output, ![H](https://www.gnu.org/software/gsl/doc/html/_images/math/af029cc216351e86964a2146cb3a293d9336ddd0.png) is stored in `A`, and ![R](https://www.gnu.org/software/gsl/doc/html/_images/math/c9b67a086a723a5eec85a146053dd0a2d2378bac.png) is stored in `B`. If `U` and `V` are provided (they may be null), the similarity transformations are stored in them. Additional workspace of length ![N](https://www.gnu.org/software/gsl/doc/html/_images/math/99f8a3c97c19212ca3eee7d00dddc8d68b2d1c27.png) is needed in `work`.



## Bidiagonalization

A general matrix ![A](https://www.gnu.org/software/gsl/doc/html/_images/math/3daf97c82cdab413e28d010c06770f224d64dc5a.png) can be factorized by similarity transformations into the form,

![A = U B V^T](https://www.gnu.org/software/gsl/doc/html/_images/math/dcb3eb33d9407c672bd0c94ff9379a55f8922a72.png)

where ![U](https://www.gnu.org/software/gsl/doc/html/_images/math/3d0e8e66abfde1273e99e86e8a090a7073d842e3.png) and ![V](https://www.gnu.org/software/gsl/doc/html/_images/math/632777989bf30753abbbf45a4a9b07ebf0fd36ad.png) are orthogonal matrices and ![B](https://www.gnu.org/software/gsl/doc/html/_images/math/39ac5dd0978410e4f612dd00070a6227cb1ceb10.png) is a ![N](https://www.gnu.org/software/gsl/doc/html/_images/math/99f8a3c97c19212ca3eee7d00dddc8d68b2d1c27.png)-by-![N](https://www.gnu.org/software/gsl/doc/html/_images/math/99f8a3c97c19212ca3eee7d00dddc8d68b2d1c27.png) bidiagonal matrix with non-zero entries only on the diagonal and superdiagonal. The size of `U` is ![M](https://www.gnu.org/software/gsl/doc/html/_images/math/85ca40542b59eef42eec1488ce2784e307962a5f.png)-by-![N](https://www.gnu.org/software/gsl/doc/html/_images/math/99f8a3c97c19212ca3eee7d00dddc8d68b2d1c27.png) and the size of `V` is ![N](https://www.gnu.org/software/gsl/doc/html/_images/math/99f8a3c97c19212ca3eee7d00dddc8d68b2d1c27.png)-by-![N](https://www.gnu.org/software/gsl/doc/html/_images/math/99f8a3c97c19212ca3eee7d00dddc8d68b2d1c27.png).

- int `gsl_linalg_bidiag_decomp`([gsl_matrix](https://www.gnu.org/software/gsl/doc/html/vectors.html#c.gsl_matrix) * *A*, [gsl_vector](https://www.gnu.org/software/gsl/doc/html/vectors.html#c.gsl_vector) * *tau_U*, [gsl_vector](https://www.gnu.org/software/gsl/doc/html/vectors.html#c.gsl_vector) * *tau_V*)

  This function factorizes the ![M](https://www.gnu.org/software/gsl/doc/html/_images/math/85ca40542b59eef42eec1488ce2784e307962a5f.png)-by-![N](https://www.gnu.org/software/gsl/doc/html/_images/math/99f8a3c97c19212ca3eee7d00dddc8d68b2d1c27.png) matrix `A` into bidiagonal form ![U B V^T](https://www.gnu.org/software/gsl/doc/html/_images/math/4e2d498bfe1129c98e341cd15b70b9318ec13957.png). The diagonal and superdiagonal of the matrix ![B](https://www.gnu.org/software/gsl/doc/html/_images/math/39ac5dd0978410e4f612dd00070a6227cb1ceb10.png) are stored in the diagonal and superdiagonal of `A`. The orthogonal matrices ![U](https://www.gnu.org/software/gsl/doc/html/_images/math/3d0e8e66abfde1273e99e86e8a090a7073d842e3.png) and `V` are stored as compressed Householder vectors in the remaining elements of `A`. The Householder coefficients are stored in the vectors `tau_U` and `tau_V`. The length of `tau_U` must equal the number of elements in the diagonal of `A` and the length of `tau_V` should be one element shorter.

- int `gsl_linalg_bidiag_unpack`(const [gsl_matrix](https://www.gnu.org/software/gsl/doc/html/vectors.html#c.gsl_matrix) * *A*, const [gsl_vector](https://www.gnu.org/software/gsl/doc/html/vectors.html#c.gsl_vector) * *tau_U*, [gsl_matrix](https://www.gnu.org/software/gsl/doc/html/vectors.html#c.gsl_matrix) * *U*, const [gsl_vector](https://www.gnu.org/software/gsl/doc/html/vectors.html#c.gsl_vector) * *tau_V*, [gsl_matrix](https://www.gnu.org/software/gsl/doc/html/vectors.html#c.gsl_matrix) * *V*, [gsl_vector](https://www.gnu.org/software/gsl/doc/html/vectors.html#c.gsl_vector) * *diag*, [gsl_vector](https://www.gnu.org/software/gsl/doc/html/vectors.html#c.gsl_vector) * *superdiag*)

  This function unpacks the bidiagonal decomposition of `A` produced by[`gsl_linalg_bidiag_decomp()`](https://www.gnu.org/software/gsl/doc/html/linalg.html#c.gsl_linalg_bidiag_decomp), (`A`, `tau_U`, `tau_V`) into the separate orthogonal matrices `U`, `V` and the diagonal vector `diag` and superdiagonal `superdiag`. Note that `U` is stored as a compact ![M](https://www.gnu.org/software/gsl/doc/html/_images/math/85ca40542b59eef42eec1488ce2784e307962a5f.png)-by-![N](https://www.gnu.org/software/gsl/doc/html/_images/math/99f8a3c97c19212ca3eee7d00dddc8d68b2d1c27.png) orthogonal matrix satisfying ![U^T U = I](https://www.gnu.org/software/gsl/doc/html/_images/math/2fc4b0d0e951a7770948a5ef54eb18aa3f8d8105.png) for efficiency.

- int `gsl_linalg_bidiag_unpack2`([gsl_matrix](https://www.gnu.org/software/gsl/doc/html/vectors.html#c.gsl_matrix) * *A*, [gsl_vector](https://www.gnu.org/software/gsl/doc/html/vectors.html#c.gsl_vector) * *tau_U*, [gsl_vector](https://www.gnu.org/software/gsl/doc/html/vectors.html#c.gsl_vector) * *tau_V*, [gsl_matrix](https://www.gnu.org/software/gsl/doc/html/vectors.html#c.gsl_matrix) * *V*)

  This function unpacks the bidiagonal decomposition of `A` produced by[`gsl_linalg_bidiag_decomp()`](https://www.gnu.org/software/gsl/doc/html/linalg.html#c.gsl_linalg_bidiag_decomp), (`A`, `tau_U`, `tau_V`) into the separate orthogonal matrices `U`, `V` and the diagonal vector `diag` and superdiagonal `superdiag`. The matrix `U` is stored in-place in `A`.

- int `gsl_linalg_bidiag_unpack_B`(const [gsl_matrix](https://www.gnu.org/software/gsl/doc/html/vectors.html#c.gsl_matrix) * *A*, [gsl_vector](https://www.gnu.org/software/gsl/doc/html/vectors.html#c.gsl_vector) * *diag*, [gsl_vector](https://www.gnu.org/software/gsl/doc/html/vectors.html#c.gsl_vector) * *superdiag*)

  This function unpacks the diagonal and superdiagonal of the bidiagonal decomposition of `A`from [`gsl_linalg_bidiag_decomp()`](https://www.gnu.org/software/gsl/doc/html/linalg.html#c.gsl_linalg_bidiag_decomp), into the diagonal vector `diag` and superdiagonal vector `superdiag`.



## Givens Rotations

A Givens rotation is a rotation in the plane acting on two elements of a given vector. It can be represented in matrix form as

![G(i,j,\theta) = \left( \begin{matrix}   1 & \ldots & 0 & \ldots & 0 & \ldots & 0 \\   \vdots & \ddots & \vdots &  & \vdots & & \vdots \\   0 & \ldots & \cos{\theta} & \ldots & -\sin{\theta} & \ldots & 0 \\   \vdots &  & \vdots & \ddots & \vdots & & \vdots \\   0 & \ldots & \sin{\theta} & \ldots & \cos{\theta} & \ldots & 0 \\   \vdots &  & \vdots &  & \vdots & \ddots & \vdots \\   0 & \ldots & 0 & \ldots & 0 & \ldots & 1 \end{matrix} \right)](https://www.gnu.org/software/gsl/doc/html/_images/math/62a311e010ad219c3f5e2ee33e743a3404015429.png)

where the ![\cos{\theta}](https://www.gnu.org/software/gsl/doc/html/_images/math/626626ba8cf5a4db7da83c2fc8fac9aefe327a5a.png) and ![\sin{\theta}](https://www.gnu.org/software/gsl/doc/html/_images/math/51274bbc0a0d7d01096f5e591221b07c126e95bf.png) appear at the intersection of the ![i](https://www.gnu.org/software/gsl/doc/html/_images/math/1cc632900aae8b2837a1f383619c0ad753be7d29.png)-th and ![j](https://www.gnu.org/software/gsl/doc/html/_images/math/7bc2ad053edeecad39b867d60495ac2caa4555dc.png)-th rows and columns. When acting on a vector ![x](https://www.gnu.org/software/gsl/doc/html/_images/math/1c5715c8edd0bbc330baa0f0494f1f7951619754.png), ![G(i,j,\theta) x](https://www.gnu.org/software/gsl/doc/html/_images/math/3bfa1f2e87eda7012d0ce00e047fff4a660a3ef8.png) performs a rotation of the ![(i,j)](https://www.gnu.org/software/gsl/doc/html/_images/math/4b982cbeb3e7a1b695da41fa331848138a69542a.png) elements of ![x](https://www.gnu.org/software/gsl/doc/html/_images/math/1c5715c8edd0bbc330baa0f0494f1f7951619754.png). Givens rotations are typically used to introduce zeros in vectors, such as during the QR decomposition of a matrix. In this case, it is typically desired to find ![c](https://www.gnu.org/software/gsl/doc/html/_images/math/a2ff0b909e008fecae54f397fc3bd6ef4ae96b5c.png) and ![s](https://www.gnu.org/software/gsl/doc/html/_images/math/dd69d21b40087e08fc6827d2f3481f57462046e7.png) such that

![\left( \begin{matrix}   c & -s \\   s & c \end{matrix} \right) \left( \begin{matrix}   a \\   b \end{matrix} \right) = \left( \begin{matrix}   r \\   0 \end{matrix} \right)](https://www.gnu.org/software/gsl/doc/html/_images/math/2a2d445a6522fa3d6975d194039f445ee055d1d9.png)

with ![r = \sqrt{a^2 + b^2}](https://www.gnu.org/software/gsl/doc/html/_images/math/95081719b734ef36744006cc07cd51a89c5428c7.png).

- void `gsl_linalg_givens`(const double *a*, const double *b*, double * *c*, double * *s*)

  This function computes ![c = \cos{\theta}](https://www.gnu.org/software/gsl/doc/html/_images/math/ac648e26758a89a4dd9f65e6f7fc56b7c945dcd5.png) and ![s = \sin{\theta}](https://www.gnu.org/software/gsl/doc/html/_images/math/8133005be09f75dbc3ecf2adc3fbf3932491a712.png) so that the Givens matrix ![G(\theta)](https://www.gnu.org/software/gsl/doc/html/_images/math/b21d86d3d5c60257a460e78a2be0083f2d6c2f0a.png) acting on the vector ![(a,b)](https://www.gnu.org/software/gsl/doc/html/_images/math/2dc21aa88e97bfdf1a42cc3b66f9429479b9f536.png) produces ![(r, 0)](https://www.gnu.org/software/gsl/doc/html/_images/math/d0fffa3c570bf9e53ffe33801ac7a214e2d25648.png), with ![r = \sqrt{a^2 + b^2}](https://www.gnu.org/software/gsl/doc/html/_images/math/95081719b734ef36744006cc07cd51a89c5428c7.png).

- void `gsl_linalg_givens_gv`([gsl_vector](https://www.gnu.org/software/gsl/doc/html/vectors.html#c.gsl_vector) * *v*, const size_t *i*, const size_t *j*, const double *c*, const double *s*)

  This function applies the Givens rotation defined by ![c = \cos{\theta}](https://www.gnu.org/software/gsl/doc/html/_images/math/ac648e26758a89a4dd9f65e6f7fc56b7c945dcd5.png) and ![s = \sin{\theta}](https://www.gnu.org/software/gsl/doc/html/_images/math/8133005be09f75dbc3ecf2adc3fbf3932491a712.png) to the `i` and `j`elements of `v`. On output, ![(v(i),v(j)) \leftarrow G(\theta) (v(i),v(j))](https://www.gnu.org/software/gsl/doc/html/_images/math/f22152ce861db65fccd6fd843ecfa67da25c078f.png).



## Householder Transformations

A Householder transformation is a rank-1 modification of the identity matrix which can be used to zero out selected elements of a vector. A Householder matrix ![P](https://www.gnu.org/software/gsl/doc/html/_images/math/d1a0a8c0064fe114c19c8fd71ffec7a98480f707.png) takes the form,

![P = I - \tau v v^T](https://www.gnu.org/software/gsl/doc/html/_images/math/9f97ee262d1e22f81c97284d92534340d9b7a79c.png)

where ![v](https://www.gnu.org/software/gsl/doc/html/_images/math/a3ab041b283bfe2226babbccf0611362d5001323.png) is a vector (called the *Householder vector*) and ![\tau = 2/(v^T v)](https://www.gnu.org/software/gsl/doc/html/_images/math/3626d82889291b7190d367a441d03ea49d5e351c.png). The functions described in this section use the rank-1 structure of the Householder matrix to create and apply Householder transformations efficiently.

- double `gsl_linalg_householder_transform`([gsl_vector](https://www.gnu.org/software/gsl/doc/html/vectors.html#c.gsl_vector) * *w*)

- gsl_complex `gsl_linalg_complex_householder_transform`(gsl_vector_complex * *w*)

  This function prepares a Householder transformation ![P = I - \tau v v^T](https://www.gnu.org/software/gsl/doc/html/_images/math/8f4934c508dfc213b9580924141e8c286d1905e1.png) which can be used to zero all the elements of the input vector `w` except the first. On output the Householder vector `v` is stored in `w` and the scalar ![\tau](https://www.gnu.org/software/gsl/doc/html/_images/math/de9fbd420c3bcc6304590e768c308e743336d58d.png) is returned. The householder vector `v` is normalized so that `v[0] = 1`, however this 1 is not stored in the output vector. Instead, `w[0]` is set to the first element of the transformed vector, so that if ![u = P w](https://www.gnu.org/software/gsl/doc/html/_images/math/a7d83a5119a1fa1c14da4b02bbde3b4d13c33912.png), `w[0] = u[0]` on output and the remainder of ![u](https://www.gnu.org/software/gsl/doc/html/_images/math/977a622c61341d61de0eafdffff2943ab85ab92d.png) is zero.

- int `gsl_linalg_householder_hm`(double *tau*, const [gsl_vector](https://www.gnu.org/software/gsl/doc/html/vectors.html#c.gsl_vector) * *v*, [gsl_matrix](https://www.gnu.org/software/gsl/doc/html/vectors.html#c.gsl_matrix) * *A*)

- int `gsl_linalg_complex_householder_hm`(gsl_complex *tau*, const gsl_vector_complex * *v*, gsl_matrix_complex * *A*)

  This function applies the Householder matrix ![P](https://www.gnu.org/software/gsl/doc/html/_images/math/d1a0a8c0064fe114c19c8fd71ffec7a98480f707.png) defined by the scalar `tau` and the vector `v` to the left-hand side of the matrix `A`. On output the result ![P A](https://www.gnu.org/software/gsl/doc/html/_images/math/e36d37bbbc4f766a889f8b2a3de896e0b36bca44.png) is stored in `A`.

- int `gsl_linalg_householder_mh`(double *tau*, const [gsl_vector](https://www.gnu.org/software/gsl/doc/html/vectors.html#c.gsl_vector) * *v*, [gsl_matrix](https://www.gnu.org/software/gsl/doc/html/vectors.html#c.gsl_matrix) * *A*)

- int `gsl_linalg_complex_householder_mh`(gsl_complex *tau*, const gsl_vector_complex * *v*, gsl_matrix_complex * *A*)

  This function applies the Householder matrix ![P](https://www.gnu.org/software/gsl/doc/html/_images/math/d1a0a8c0064fe114c19c8fd71ffec7a98480f707.png) defined by the scalar `tau` and the vector `v` to the right-hand side of the matrix `A`. On output the result ![A P](https://www.gnu.org/software/gsl/doc/html/_images/math/677a652bb0ba6dfdd754330869340450aeaa2d7c.png) is stored in `A`.

- int `gsl_linalg_householder_hv`(double *tau*, const [gsl_vector](https://www.gnu.org/software/gsl/doc/html/vectors.html#c.gsl_vector) * *v*, [gsl_vector](https://www.gnu.org/software/gsl/doc/html/vectors.html#c.gsl_vector) * *w*)

- int `gsl_linalg_complex_householder_hv`(gsl_complex *tau*, const gsl_vector_complex * *v*, gsl_vector_complex * *w*)

  This function applies the Householder transformation ![P](https://www.gnu.org/software/gsl/doc/html/_images/math/d1a0a8c0064fe114c19c8fd71ffec7a98480f707.png) defined by the scalar `tau` and the vector `v` to the vector `w`. On output the result ![P w](https://www.gnu.org/software/gsl/doc/html/_images/math/8308df1687dfa9db16ec0bc25252e8e5764caf40.png) is stored in `w`.



## Householder solver for linear systems

- int `gsl_linalg_HH_solve`([gsl_matrix](https://www.gnu.org/software/gsl/doc/html/vectors.html#c.gsl_matrix) * *A*, const [gsl_vector](https://www.gnu.org/software/gsl/doc/html/vectors.html#c.gsl_vector) * *b*, [gsl_vector](https://www.gnu.org/software/gsl/doc/html/vectors.html#c.gsl_vector) * *x*)

  This function solves the system ![A x = b](https://www.gnu.org/software/gsl/doc/html/_images/math/bf1b0fe2f3273671dbd86939aa2ca8992c570e65.png) directly using Householder transformations. On output the solution is stored in `x` and `b` is not modified. The matrix `A` is destroyed by the Householder transformations.

- int `gsl_linalg_HH_svx`([gsl_matrix](https://www.gnu.org/software/gsl/doc/html/vectors.html#c.gsl_matrix) * *A*, [gsl_vector](https://www.gnu.org/software/gsl/doc/html/vectors.html#c.gsl_vector) * *x*)

  This function solves the system ![A x = b](https://www.gnu.org/software/gsl/doc/html/_images/math/bf1b0fe2f3273671dbd86939aa2ca8992c570e65.png) in-place using Householder transformations. On input `x` should contain the right-hand side ![b](https://www.gnu.org/software/gsl/doc/html/_images/math/88082f6488798774b9cd95a59c1ab4c6cafe909f.png), which is replaced by the solution on output. The matrix `A` is destroyed by the Householder transformations.



## Tridiagonal Systems

The functions described in this section efficiently solve symmetric, non-symmetric and cyclic tridiagonal systems with minimal storage. Note that the current implementations of these functions use a variant of Cholesky decomposition, so the tridiagonal matrix must be positive definite. For non-positive definite matrices, the functions return the error code `GSL_ESING`.

- int `gsl_linalg_solve_tridiag`(const [gsl_vector](https://www.gnu.org/software/gsl/doc/html/vectors.html#c.gsl_vector) * *diag*, const [gsl_vector](https://www.gnu.org/software/gsl/doc/html/vectors.html#c.gsl_vector) * *e*, const [gsl_vector](https://www.gnu.org/software/gsl/doc/html/vectors.html#c.gsl_vector) * *f*, const [gsl_vector](https://www.gnu.org/software/gsl/doc/html/vectors.html#c.gsl_vector) * *b*, [gsl_vector](https://www.gnu.org/software/gsl/doc/html/vectors.html#c.gsl_vector) * *x*)

  This function solves the general ![N](https://www.gnu.org/software/gsl/doc/html/_images/math/99f8a3c97c19212ca3eee7d00dddc8d68b2d1c27.png)-by-![N](https://www.gnu.org/software/gsl/doc/html/_images/math/99f8a3c97c19212ca3eee7d00dddc8d68b2d1c27.png) system ![A x = b](https://www.gnu.org/software/gsl/doc/html/_images/math/bf1b0fe2f3273671dbd86939aa2ca8992c570e65.png) where `A` is tridiagonal (![N \geq 2](https://www.gnu.org/software/gsl/doc/html/_images/math/9569337fe93593e8509bbe5929e48fe5085524ea.png)). The super-diagonal and sub-diagonal vectors `e` and `f` must be one element shorter than the diagonal vector `diag`. The form of `A` for the 4-by-4 case is shown below,![A = \left( \begin{matrix}   d_0&e_0&  0& 0\\   f_0&d_1&e_1& 0\\   0  &f_1&d_2&e_2\\   0  &0  &f_2&d_3 \end{matrix} \right)](https://www.gnu.org/software/gsl/doc/html/_images/math/37f068fef06b74534277c0d5b843b2e2ff1b0ad3.png)

- int `gsl_linalg_solve_symm_tridiag`(const [gsl_vector](https://www.gnu.org/software/gsl/doc/html/vectors.html#c.gsl_vector) * *diag*, const [gsl_vector](https://www.gnu.org/software/gsl/doc/html/vectors.html#c.gsl_vector) * *e*, const [gsl_vector](https://www.gnu.org/software/gsl/doc/html/vectors.html#c.gsl_vector) * *b*, [gsl_vector](https://www.gnu.org/software/gsl/doc/html/vectors.html#c.gsl_vector) * *x*)

  This function solves the general ![N](https://www.gnu.org/software/gsl/doc/html/_images/math/99f8a3c97c19212ca3eee7d00dddc8d68b2d1c27.png)-by-![N](https://www.gnu.org/software/gsl/doc/html/_images/math/99f8a3c97c19212ca3eee7d00dddc8d68b2d1c27.png) system ![A x = b](https://www.gnu.org/software/gsl/doc/html/_images/math/bf1b0fe2f3273671dbd86939aa2ca8992c570e65.png) where `A` is symmetric tridiagonal (![N \geq 2](https://www.gnu.org/software/gsl/doc/html/_images/math/9569337fe93593e8509bbe5929e48fe5085524ea.png)). The off-diagonal vector `e` must be one element shorter than the diagonal vector `diag`. The form of `A` for the 4-by-4 case is shown below,![A = \left( \begin{matrix}   d_0&e_0&  0& 0\\   e_0&d_1&e_1& 0\\   0  &e_1&d_2&e_2\\   0  &0  &e_2&d_3 \end{matrix} \right)](https://www.gnu.org/software/gsl/doc/html/_images/math/75f09e03ded45a7d209f84407f7a254a0b54f539.png)

- int `gsl_linalg_solve_cyc_tridiag`(const [gsl_vector](https://www.gnu.org/software/gsl/doc/html/vectors.html#c.gsl_vector) * *diag*, const [gsl_vector](https://www.gnu.org/software/gsl/doc/html/vectors.html#c.gsl_vector) * *e*, const [gsl_vector](https://www.gnu.org/software/gsl/doc/html/vectors.html#c.gsl_vector) * *f*, const [gsl_vector](https://www.gnu.org/software/gsl/doc/html/vectors.html#c.gsl_vector) * *b*, [gsl_vector](https://www.gnu.org/software/gsl/doc/html/vectors.html#c.gsl_vector) * *x*)

  This function solves the general ![N](https://www.gnu.org/software/gsl/doc/html/_images/math/99f8a3c97c19212ca3eee7d00dddc8d68b2d1c27.png)-by-![N](https://www.gnu.org/software/gsl/doc/html/_images/math/99f8a3c97c19212ca3eee7d00dddc8d68b2d1c27.png) system ![A x = b](https://www.gnu.org/software/gsl/doc/html/_images/math/bf1b0fe2f3273671dbd86939aa2ca8992c570e65.png) where `A` is cyclic tridiagonal (![N \geq 3](https://www.gnu.org/software/gsl/doc/html/_images/math/9022013752a3d81ace3c7de5697c9de906f12719.png)). The cyclic super-diagonal and sub-diagonal vectors `e` and `f` must have the same number of elements as the diagonal vector `diag`. The form of `A` for the 4-by-4 case is shown below,![A = \left( \begin{matrix}   d_0&e_0& 0 &f_3\\   f_0&d_1&e_1& 0 \\   0 &f_1&d_2&e_2\\   e_3& 0 &f_2&d_3 \end{matrix} \right)](https://www.gnu.org/software/gsl/doc/html/_images/math/efc1a18f5493abbe66ac41451c8d787a257ef965.png)

- int `gsl_linalg_solve_symm_cyc_tridiag`(const [gsl_vector](https://www.gnu.org/software/gsl/doc/html/vectors.html#c.gsl_vector) * *diag*, const [gsl_vector](https://www.gnu.org/software/gsl/doc/html/vectors.html#c.gsl_vector) * *e*, const [gsl_vector](https://www.gnu.org/software/gsl/doc/html/vectors.html#c.gsl_vector) * *b*, [gsl_vector](https://www.gnu.org/software/gsl/doc/html/vectors.html#c.gsl_vector) * *x*)

  This function solves the general ![N](https://www.gnu.org/software/gsl/doc/html/_images/math/99f8a3c97c19212ca3eee7d00dddc8d68b2d1c27.png)-by-![N](https://www.gnu.org/software/gsl/doc/html/_images/math/99f8a3c97c19212ca3eee7d00dddc8d68b2d1c27.png) system ![A x = b](https://www.gnu.org/software/gsl/doc/html/_images/math/bf1b0fe2f3273671dbd86939aa2ca8992c570e65.png) where `A` is symmetric cyclic tridiagonal (![N \geq 3](https://www.gnu.org/software/gsl/doc/html/_images/math/9022013752a3d81ace3c7de5697c9de906f12719.png)). The cyclic off-diagonal vector `e` must have the same number of elements as the diagonal vector `diag`. The form of `A` for the 4-by-4 case is shown below,![A = \left( \begin{matrix}   d_0&e_0& 0 &e_3\\   e_0&d_1&e_1& 0 \\   0 &e_1&d_2&e_2\\   e_3& 0 &e_2&d_3 \end{matrix} \right)](https://www.gnu.org/software/gsl/doc/html/_images/math/f07ca0c017bc4ca05e11503f19b909425214026b.png)



## Triangular Systems

- int `gsl_linalg_tri_upper_invert`([gsl_matrix](https://www.gnu.org/software/gsl/doc/html/vectors.html#c.gsl_matrix) * *T*)

- int `gsl_linalg_tri_lower_invert`([gsl_matrix](https://www.gnu.org/software/gsl/doc/html/vectors.html#c.gsl_matrix) * *T*)

- int `gsl_linalg_tri_upper_unit_invert`([gsl_matrix](https://www.gnu.org/software/gsl/doc/html/vectors.html#c.gsl_matrix) * *T*)

- int `gsl_linalg_tri_lower_unit_invert`([gsl_matrix](https://www.gnu.org/software/gsl/doc/html/vectors.html#c.gsl_matrix) * *T*)

  These functions calculate the in-place inverse of the triangular matrix `T`. When the `upper`prefix is specified, then the upper triangle of `T` is used, and when the `lower` prefix is specified, the lower triangle is used. If the `unit` prefix is specified, then the diagonal elements of the matrix `T` are taken as unity and are not referenced. Otherwise the diagonal elements are used in the inversion.

- int `gsl_linalg_tri_upper_rcond`(const [gsl_matrix](https://www.gnu.org/software/gsl/doc/html/vectors.html#c.gsl_matrix) * *T*, double * *rcond*, [gsl_vector](https://www.gnu.org/software/gsl/doc/html/vectors.html#c.gsl_vector) * *work*)

- int `gsl_linalg_tri_lower_rcond`(const [gsl_matrix](https://www.gnu.org/software/gsl/doc/html/vectors.html#c.gsl_matrix) * *T*, double * *rcond*, [gsl_vector](https://www.gnu.org/software/gsl/doc/html/vectors.html#c.gsl_vector) * *work*)

  These functions estimate the reciprocal condition number, in the 1-norm, of the upper or lower ![N](https://www.gnu.org/software/gsl/doc/html/_images/math/99f8a3c97c19212ca3eee7d00dddc8d68b2d1c27.png)-by-![N](https://www.gnu.org/software/gsl/doc/html/_images/math/99f8a3c97c19212ca3eee7d00dddc8d68b2d1c27.png) triangular matrix `T`. The reciprocal condition number is stored in `rcond` on output, and is defined by ![1 / (||T||_1 \cdot ||T^{-1}||_1)](https://www.gnu.org/software/gsl/doc/html/_images/math/a7e56b69a7e6b9f7a9c05171100e87ff8b072c9f.png). Additional workspace of size ![3 N](https://www.gnu.org/software/gsl/doc/html/_images/math/1b2b0ffdf792c9f396c727992622338aadd67f74.png) is required in `work`.





## Balancing

The process of balancing a matrix applies similarity transformations to make the rows and columns have comparable norms. This is useful, for example, to reduce roundoff errors in the solution of eigenvalue problems. Balancing a matrix ![A](https://www.gnu.org/software/gsl/doc/html/_images/math/3daf97c82cdab413e28d010c06770f224d64dc5a.png) consists of replacing ![A](https://www.gnu.org/software/gsl/doc/html/_images/math/3daf97c82cdab413e28d010c06770f224d64dc5a.png) with a similar matrix

![A' = D^{-1} A D](https://www.gnu.org/software/gsl/doc/html/_images/math/a59581fccae909bfed11c143b974f8668795c71e.png)

where ![D](https://www.gnu.org/software/gsl/doc/html/_images/math/0638e3c183dd8cfd1790c25e762f1a8c53ee5711.png) is a diagonal matrix whose entries are powers of the floating point radix.

- int `gsl_linalg_balance_matrix`([gsl_matrix](https://www.gnu.org/software/gsl/doc/html/vectors.html#c.gsl_matrix) * *A*, [gsl_vector](https://www.gnu.org/software/gsl/doc/html/vectors.html#c.gsl_vector) * *D*)

  This function replaces the matrix `A` with its balanced counterpart and stores the diagonal elements of the similarity transformation into the vector `D`.

## Examples

The following program solves the linear system ![A x = b](https://www.gnu.org/software/gsl/doc/html/_images/math/bf1b0fe2f3273671dbd86939aa2ca8992c570e65.png). The system to be solved is,

![\left( \begin{matrix}   0.18& 0.60& 0.57& 0.96\\   0.41& 0.24& 0.99& 0.58\\   0.14& 0.30& 0.97& 0.66\\   0.51& 0.13& 0.19& 0.85 \end{matrix} \right) \left( \begin{matrix}   x_0\\   x_1\\   x_2\\   x_3 \end{matrix} \right) = \left( \begin{matrix}   1.0\\   2.0\\   3.0\\   4.0 \end{matrix} \right)](https://www.gnu.org/software/gsl/doc/html/_images/math/f0d3c8b48a338ced25b57cb2057b6c07e1aa24a0.png)

and the solution is found using LU decomposition of the matrix ![A](https://www.gnu.org/software/gsl/doc/html/_images/math/3daf97c82cdab413e28d010c06770f224d64dc5a.png).

```
#include <stdio.h>
#include <gsl/gsl_linalg.h>

int
main (void)
{
  double a_data[] = { 0.18, 0.60, 0.57, 0.96,
                      0.41, 0.24, 0.99, 0.58,
                      0.14, 0.30, 0.97, 0.66,
                      0.51, 0.13, 0.19, 0.85 };

  double b_data[] = { 1.0, 2.0, 3.0, 4.0 };

  gsl_matrix_view m
    = gsl_matrix_view_array (a_data, 4, 4);

  gsl_vector_view b
    = gsl_vector_view_array (b_data, 4);

  gsl_vector *x = gsl_vector_alloc (4);

  int s;

  gsl_permutation * p = gsl_permutation_alloc (4);

  gsl_linalg_LU_decomp (&m.matrix, p, &s);

  gsl_linalg_LU_solve (&m.matrix, p, &b.vector, x);

  printf ("x = \n");
  gsl_vector_fprintf (stdout, x, "%g");

  gsl_permutation_free (p);
  gsl_vector_free (x);
  return 0;
}
```

Here is the output from the program,

```
x =
-4.05205
-12.6056
1.66091
8.69377
```

This can be verified by multiplying the solution ![x](https://www.gnu.org/software/gsl/doc/html/_images/math/1c5715c8edd0bbc330baa0f0494f1f7951619754.png) by the original matrix ![A](https://www.gnu.org/software/gsl/doc/html/_images/math/3daf97c82cdab413e28d010c06770f224d64dc5a.png) using GNU octave,

```
octave> A = [ 0.18, 0.60, 0.57, 0.96;
              0.41, 0.24, 0.99, 0.58;
              0.14, 0.30, 0.97, 0.66;
              0.51, 0.13, 0.19, 0.85 ];

octave> x = [ -4.05205; -12.6056; 1.66091; 8.69377];

octave> A * x
ans =
  1.0000
  2.0000
  3.0000
  4.0000
```

This reproduces the original right-hand side vector, ![b](https://www.gnu.org/software/gsl/doc/html/_images/math/88082f6488798774b9cd95a59c1ab4c6cafe909f.png), in accordance with the equation ![A x = b](https://www.gnu.org/software/gsl/doc/html/_images/math/bf1b0fe2f3273671dbd86939aa2ca8992c570e65.png).

## References and Further Reading

Further information on the algorithms described in this section can be found in the following book,

- G. H. Golub, C. F. Van Loan, “Matrix Computations” (3rd Ed, 1996), Johns Hopkins University Press, ISBN 0-8018-5414-8.

The LAPACK library is described in the following manual,

- *LAPACK Users’ Guide* (Third Edition, 1999), Published by SIAM, ISBN 0-89871-447-8

The LAPACK source code can be found at <http://www.netlib.org/lapack>, along with an online copy of the users guide.

The Modified Golub-Reinsch algorithm is described in the following paper,

- T.F. Chan, “An Improved Algorithm for Computing the Singular Value Decomposition”, ACM Transactions on Mathematical Software, 8 (1982), pp 72–83.

The Jacobi algorithm for singular value decomposition is described in the following papers,

- J.C. Nash, “A one-sided transformation method for the singular value decomposition and algebraic eigenproblem”, Computer Journal, Volume 18, Number 1 (1975), p 74–76
- J.C. Nash and S. Shlien “Simple algorithms for the partial singular value decomposition”, Computer Journal, Volume 30 (1987), p 268–275.
- J. Demmel, K. Veselic, “Jacobi’s Method is more accurate than QR”, Lapack Working Note 15 (LAWN-15), October 1989. Available from netlib, <http://www.netlib.org/lapack/> in the `lawns` or`lawnspdf` directories.

The algorithm for estimating a matrix condition number is described in the following paper,

- N. J. Higham, “FORTRAN codes for estimating the one-norm of a real or complex matrix, with applications to condition estimation”, ACM Trans. Math. Soft., vol. 14, no. 4, pp. 381-396, December 1988.
