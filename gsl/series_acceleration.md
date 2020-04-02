# Sparse Linear Algebra

This chapter describes functions for solving sparse linear systems of equations. The library provides linear algebra routines which operate directly on the [`gsl_spmatrix`](https://www.gnu.org/software/gsl/doc/html/spmatrix.html#c.gsl_spmatrix) and [`gsl_vector`](https://www.gnu.org/software/gsl/doc/html/vectors.html#c.gsl_vector) objects.

The functions described in this chapter are declared in the header file `gsl_splinalg.h`.



## Overview

This chapter is primarily concerned with the solution of the linear system

![A x = b](https://www.gnu.org/software/gsl/doc/html/_images/math/fadeea787fcb81a5f3b14b5e8f271a578720a4df.png)

where ![A](https://www.gnu.org/software/gsl/doc/html/_images/math/1c185d194f524d92bebad44e7b5fad9e4bded1d6.png) is a general square ![n](https://www.gnu.org/software/gsl/doc/html/_images/math/46ebcf5e3ac0161a5a3a0e7a494897e29c9bda2f.png)-by-![n](https://www.gnu.org/software/gsl/doc/html/_images/math/46ebcf5e3ac0161a5a3a0e7a494897e29c9bda2f.png) non-singular sparse matrix, ![x](https://www.gnu.org/software/gsl/doc/html/_images/math/b3da445c5b7bba3e36703d3763559ee8c510a111.png) is an unknown ![n](https://www.gnu.org/software/gsl/doc/html/_images/math/46ebcf5e3ac0161a5a3a0e7a494897e29c9bda2f.png)-by-![1](https://www.gnu.org/software/gsl/doc/html/_images/math/d1c44008a9848af77e068c9a70b57c975dc75ff6.png) vector, and ![b](https://www.gnu.org/software/gsl/doc/html/_images/math/97e8999ad6b5f2bbf1e37c90eee39c03aa01f0d4.png) is a given ![n](https://www.gnu.org/software/gsl/doc/html/_images/math/46ebcf5e3ac0161a5a3a0e7a494897e29c9bda2f.png)-by-1 right hand side vector. There exist many methods for solving such sparse linear systems, which broadly fall into either direct or iterative categories. Direct methods include LU and QR decompositions, while iterative methods start with an initial guess for the vector ![x](https://www.gnu.org/software/gsl/doc/html/_images/math/b3da445c5b7bba3e36703d3763559ee8c510a111.png) and update the guess through iteration until convergence. GSL does not currently provide any direct sparse solvers.



## Sparse Iterative Solvers

### Overview

Many practical iterative methods of solving large ![n](https://www.gnu.org/software/gsl/doc/html/_images/math/46ebcf5e3ac0161a5a3a0e7a494897e29c9bda2f.png)-by-![n](https://www.gnu.org/software/gsl/doc/html/_images/math/46ebcf5e3ac0161a5a3a0e7a494897e29c9bda2f.png) sparse linear systems involve projecting an approximate solution for `x` onto a subspace of ![{\bf R}^n](https://www.gnu.org/software/gsl/doc/html/_images/math/b7ea0d29b1b973641d4fb095334888e222c52c12.png). If we define a ![m](https://www.gnu.org/software/gsl/doc/html/_images/math/754fb0a482f2ae5e7ed35b7e606eddda70d561c3.png)-dimensional subspace ![{\cal K}](https://www.gnu.org/software/gsl/doc/html/_images/math/718868fc7d4be0fc4a7d76d7d64d55c2f9f0e614.png) as the subspace of approximations to the solution `x`, then ![m](https://www.gnu.org/software/gsl/doc/html/_images/math/754fb0a482f2ae5e7ed35b7e606eddda70d561c3.png) constraints must be imposed to determine the next approximation. These ![m](https://www.gnu.org/software/gsl/doc/html/_images/math/754fb0a482f2ae5e7ed35b7e606eddda70d561c3.png) constraints define another ![m](https://www.gnu.org/software/gsl/doc/html/_images/math/754fb0a482f2ae5e7ed35b7e606eddda70d561c3.png)-dimensional subspace denoted by ![{\cal L}](https://www.gnu.org/software/gsl/doc/html/_images/math/89518071f0cadb02a72c2e2f4c16d8ce417a832f.png). The subspace dimension ![m](https://www.gnu.org/software/gsl/doc/html/_images/math/754fb0a482f2ae5e7ed35b7e606eddda70d561c3.png) is typically chosen to be much smaller than ![n](https://www.gnu.org/software/gsl/doc/html/_images/math/46ebcf5e3ac0161a5a3a0e7a494897e29c9bda2f.png) in order to reduce the computational effort needed to generate the next approximate solution vector. The many iterative algorithms which exist differ mainly in their choice of ![{\cal K}](https://www.gnu.org/software/gsl/doc/html/_images/math/718868fc7d4be0fc4a7d76d7d64d55c2f9f0e614.png) and ![{\cal L}](https://www.gnu.org/software/gsl/doc/html/_images/math/89518071f0cadb02a72c2e2f4c16d8ce417a832f.png).

### Types of Sparse Iterative Solvers

The sparse linear algebra library provides the following types of iterative solvers:

- `gsl_splinalg_itersolve_type`

  `gsl_splinalg_itersolve_gmres`This specifies the Generalized Minimum Residual Method (GMRES). This is a projection method using ![{\cal K} = {\cal K}_m](https://www.gnu.org/software/gsl/doc/html/_images/math/a75552e2be0fcadaeb5798e19b01f007652ded2c.png) and ![{\cal L} = A {\cal K}_m](https://www.gnu.org/software/gsl/doc/html/_images/math/c98d9d2f8af15e6c45187961c4c8a04270405650.png) where ![{\cal K}_m](https://www.gnu.org/software/gsl/doc/html/_images/math/a204bf6d5d1a9e2d13b1824408b51c15022c4524.png) is the ![m](https://www.gnu.org/software/gsl/doc/html/_images/math/754fb0a482f2ae5e7ed35b7e606eddda70d561c3.png)-th Krylov subspace![{\cal K}_m = span \left\{ r_0, A r_0, ..., A^{m-1} r_0 \right\}](https://www.gnu.org/software/gsl/doc/html/_images/math/6b954859c8029ed9646dea742a69043ca95c9ba9.png)and ![r_0 = b - A x_0](https://www.gnu.org/software/gsl/doc/html/_images/math/49263b962b75da67181c6eed38e272eecd79a974.png) is the residual vector of the initial guess ![x_0](https://www.gnu.org/software/gsl/doc/html/_images/math/d9f93f759775cbe63c810ade453d4537fd7e73af.png). If ![m](https://www.gnu.org/software/gsl/doc/html/_images/math/754fb0a482f2ae5e7ed35b7e606eddda70d561c3.png) is set equal to ![n](https://www.gnu.org/software/gsl/doc/html/_images/math/46ebcf5e3ac0161a5a3a0e7a494897e29c9bda2f.png), then the Krylov subspace is ![{\bf R}^n](https://www.gnu.org/software/gsl/doc/html/_images/math/b7ea0d29b1b973641d4fb095334888e222c52c12.png) and GMRES will provide the exact solution `x`. However, the goal is for the method to arrive at a very good approximation to `x` using a much smaller subspace ![{\cal K}_m](https://www.gnu.org/software/gsl/doc/html/_images/math/a204bf6d5d1a9e2d13b1824408b51c15022c4524.png). By default, the GMRES method selects ![m = MIN(n,10)](https://www.gnu.org/software/gsl/doc/html/_images/math/b57fb0fc9352280436387bf8e709b7143121df48.png) but the user may specify a different value for ![m](https://www.gnu.org/software/gsl/doc/html/_images/math/754fb0a482f2ae5e7ed35b7e606eddda70d561c3.png). The GMRES storage requirements grow as ![O(n(m+1))](https://www.gnu.org/software/gsl/doc/html/_images/math/a90aa421cf36f7a521743b8c35b91a6dfbb5e30f.png) and the number of flops grow as ![O(4 m^2 n - 4 m^3 / 3)](https://www.gnu.org/software/gsl/doc/html/_images/math/a61260d38cd6eaeca40d843e6a858e94428bb0c8.png).In the below function [`gsl_splinalg_itersolve_iterate()`](https://www.gnu.org/software/gsl/doc/html/splinalg.html#c.gsl_splinalg_itersolve_iterate), one GMRES iteration is defined as projecting the approximate solution vector onto each Krylov subspace ![{\cal K}_1, ..., {\cal K}_m](https://www.gnu.org/software/gsl/doc/html/_images/math/38bec3894820f6f32fb06d5ecf449c2610896009.png), and so ![m](https://www.gnu.org/software/gsl/doc/html/_images/math/754fb0a482f2ae5e7ed35b7e606eddda70d561c3.png) can be kept smaller by “restarting” the method and calling [`gsl_splinalg_itersolve_iterate()`](https://www.gnu.org/software/gsl/doc/html/splinalg.html#c.gsl_splinalg_itersolve_iterate) multiple times, providing the updated approximation `x` to each new call. If the method is not adequately converging, the user may try increasing the parameter ![m](https://www.gnu.org/software/gsl/doc/html/_images/math/754fb0a482f2ae5e7ed35b7e606eddda70d561c3.png).GMRES is considered a robust general purpose iterative solver, however there are cases where the method stagnates if the matrix is not positive-definite and fails to reduce the residual until the very last projection onto the subspace ![{\cal K}_n = {\bf R}^n](https://www.gnu.org/software/gsl/doc/html/_images/math/e7f0c7d1ff9c1b6c26943ec1692880cc3799cac9.png). In these cases, preconditioning the linear system can help, but GSL does not currently provide any preconditioners.

### Iterating the Sparse Linear System

The following functions are provided to allocate storage for the sparse linear solvers and iterate the system to a solution.

- gsl_splinalg_itersolve * `gsl_splinalg_itersolve_alloc`(const [gsl_splinalg_itersolve_type](https://www.gnu.org/software/gsl/doc/html/splinalg.html#c.gsl_splinalg_itersolve_type) * *T*, const size_t *n*, const size_t *m*)

  This function allocates a workspace for the iterative solution of `n`-by-`n` sparse matrix systems. The iterative solver type is specified by `T`. The argument `m` specifies the size of the solution candidate subspace ![{\cal K}_m](https://www.gnu.org/software/gsl/doc/html/_images/math/a204bf6d5d1a9e2d13b1824408b51c15022c4524.png). The dimension `m` may be set to 0 in which case a reasonable default value is used.

- void `gsl_splinalg_itersolve_free`(gsl_splinalg_itersolve * *w*)

  This function frees the memory associated with the workspace `w`.

- const char * `gsl_splinalg_itersolve_name`(const gsl_splinalg_itersolve * *w*)

  This function returns a string pointer to the name of the solver.

- int `gsl_splinalg_itersolve_iterate`(const [gsl_spmatrix](https://www.gnu.org/software/gsl/doc/html/spmatrix.html#c.gsl_spmatrix) * *A*, const [gsl_vector](https://www.gnu.org/software/gsl/doc/html/vectors.html#c.gsl_vector) * *b*, const double *tol*, [gsl_vector](https://www.gnu.org/software/gsl/doc/html/vectors.html#c.gsl_vector) * *x*, gsl_splinalg_itersolve * *w*)

  This function performs one iteration of the iterative method for the sparse linear system specfied by the matrix `A`, right hand side vector `b` and solution vector `x`. On input, `x` must be set to an initial guess for the solution. On output, `x` is updated to give the current solution estimate. The parameter `tol` specifies the relative tolerance between the residual norm and norm of `b` in order to check for convergence. When the following condition is satisfied:![|| A x - b || \le tol \times || b ||](https://www.gnu.org/software/gsl/doc/html/_images/math/370fead79b5d8946a5688ca1c3d65de02d0b47f6.png)the method has converged, the function returns `GSL_SUCCESS` and the final solution is provided in `x`. Otherwise, the function returns `GSL_CONTINUE` to signal that more iterations are required. Here, ![|| \cdot ||](https://www.gnu.org/software/gsl/doc/html/_images/math/3ea77b20c54643fbaeb23e546f5f9f98f6670a80.png) represents the Euclidean norm. The input matrix `A` may be in triplet or compressed format.

- double `gsl_splinalg_itersolve_normr`(const gsl_splinalg_itersolve * *w*)

  This function returns the current residual norm ![||r|| = ||A x - b||](https://www.gnu.org/software/gsl/doc/html/_images/math/cc0e3e1b29dfe5628a6824f869eff10885848215.png), which is updated after each call to [`gsl_splinalg_itersolve_iterate()`](https://www.gnu.org/software/gsl/doc/html/splinalg.html#c.gsl_splinalg_itersolve_iterate).



## Examples

This example program demonstrates the sparse linear algebra routines on the solution of a simple 1D Poisson equation on ![[0,1]](https://www.gnu.org/software/gsl/doc/html/_images/math/2fc2df08bced2a54db898db367f3cd3befe03612.png):

![{d^2 u(x) \over dx^2} = f(x) = -\pi^2 \sin{(\pi x)}](https://www.gnu.org/software/gsl/doc/html/_images/math/905c01ea65ce6adfb6f42463f1cf6eb1572fcf91.png)

with boundary conditions ![u(0) = u(1) = 0](https://www.gnu.org/software/gsl/doc/html/_images/math/98f5a232b1d6ca2ac5d3bda70c61aa8200f68211.png). The analytic solution of this simple problem is ![u(x) = \sin{\pi x}](https://www.gnu.org/software/gsl/doc/html/_images/math/94c75244868290831933ca18379cfe2cdac083b3.png). We will solve this problem by finite differencing the left hand side to give

![{1 \over h^2} \left( u_{i+1} - 2 u_i + u_{i-1} \right) = f_i](https://www.gnu.org/software/gsl/doc/html/_images/math/4d7616c923a54de5b51d9f3485d2c21a43894384.png)

Defining a grid of ![N](https://www.gnu.org/software/gsl/doc/html/_images/math/79fd8768fbc05d51eb66ac61928673e6321a2d59.png) points, ![h = 1/(N-1)](https://www.gnu.org/software/gsl/doc/html/_images/math/a5e6f8a577c716283b46211d7c62fac77b8e1e99.png). In the finite difference equation above, ![u_0 = u_{N-1} = 0](https://www.gnu.org/software/gsl/doc/html/_images/math/3631e37bbac4b7ed60e2103ca730add186a758a9.png) are known from the boundary conditions, so we will only put the equations for ![i = 1, ..., N-2](https://www.gnu.org/software/gsl/doc/html/_images/math/85bb0bd2f0baaee2abe819a5f47c5beac9edb391.png) into the matrix system. The resulting ![(N-2) \times (N-2)](https://www.gnu.org/software/gsl/doc/html/_images/math/e690c9cc42c3d5e4f6fcea581fe5770556829660.png) matrix equation is

![{1 \over h^2} \left(   \begin{array}{cccccc}     -2 & 1 & 0 & 0 & \ldots & 0 \\     1 & -2 & 1 & 0 & \ldots & 0 \\     0 & 1 & -2 & 1 & \ldots & 0 \\     \vdots & \vdots & \ddots & \ddots & \ddots & \vdots \\     0 & \ldots & \ldots & 1 & -2 & 1 \\     0 & \ldots & \ldots & \ldots & 1 & -2   \end{array} \right) \left(   \begin{array}{c}     u_1 \\     u_2 \\     u_3 \\     \vdots \\     u_{N-3} \\     u_{N-2}   \end{array} \right) = \left(   \begin{array}{c}     f_1 \\     f_2 \\     f_3 \\     \vdots \\     f_{N-3} \\     f_{N-2}   \end{array} \right)](https://www.gnu.org/software/gsl/doc/html/_images/math/c265dae750da78fce4aacfd52fd8e6e44d6daea1.png)

An example program which constructs and solves this system is given below. The system is solved using the iterative GMRES solver. Here is the output of the program:

```
iter 0 residual = 4.297275996844e-11
Converged
```

showing that the method converged in a single iteration. The calculated solution is shown in [Fig. 40](https://www.gnu.org/software/gsl/doc/html/splinalg.html#fig-splinalg-poisson).

![_images/sparse_poisson.png](https://www.gnu.org/software/gsl/doc/html/_images/sparse_poisson.png)

Fig. 40 Solution of PDE

```
#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include <gsl/gsl_math.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_spmatrix.h>
#include <gsl/gsl_splinalg.h>

int
main()
{
  const size_t N = 100;                       /* number of grid points */
  const size_t n = N - 2;                     /* subtract 2 to exclude boundaries */
  const double h = 1.0 / (N - 1.0);           /* grid spacing */
  gsl_spmatrix *A = gsl_spmatrix_alloc(n ,n); /* triplet format */
  gsl_spmatrix *C;                            /* compressed format */
  gsl_vector *f = gsl_vector_alloc(n);        /* right hand side vector */
  gsl_vector *u = gsl_vector_alloc(n);        /* solution vector */
  size_t i;

  /* construct the sparse matrix for the finite difference equation */

  /* construct first row */
  gsl_spmatrix_set(A, 0, 0, -2.0);
  gsl_spmatrix_set(A, 0, 1, 1.0);

  /* construct rows [1:n-2] */
  for (i = 1; i < n - 1; ++i)
    {
      gsl_spmatrix_set(A, i, i + 1, 1.0);
      gsl_spmatrix_set(A, i, i, -2.0);
      gsl_spmatrix_set(A, i, i - 1, 1.0);
    }

  /* construct last row */
  gsl_spmatrix_set(A, n - 1, n - 1, -2.0);
  gsl_spmatrix_set(A, n - 1, n - 2, 1.0);

  /* scale by h^2 */
  gsl_spmatrix_scale(A, 1.0 / (h * h));

  /* construct right hand side vector */
  for (i = 0; i < n; ++i)
    {
      double xi = (i + 1) * h;
      double fi = -M_PI * M_PI * sin(M_PI * xi);
      gsl_vector_set(f, i, fi);
    }

  /* convert to compressed column format */
  C = gsl_spmatrix_ccs(A);

  /* now solve the system with the GMRES iterative solver */
  {
    const double tol = 1.0e-6;  /* solution relative tolerance */
    const size_t max_iter = 10; /* maximum iterations */
    const gsl_splinalg_itersolve_type *T = gsl_splinalg_itersolve_gmres;
    gsl_splinalg_itersolve *work =
      gsl_splinalg_itersolve_alloc(T, n, 0);
    size_t iter = 0;
    double residual;
    int status;

    /* initial guess u = 0 */
    gsl_vector_set_zero(u);

    /* solve the system A u = f */
    do
      {
        status = gsl_splinalg_itersolve_iterate(C, f, tol, u, work);

        /* print out residual norm ||A*u - f|| */
        residual = gsl_splinalg_itersolve_normr(work);
        fprintf(stderr, "iter %zu residual = %.12e\n", iter, residual);

        if (status == GSL_SUCCESS)
          fprintf(stderr, "Converged\n");
      }
    while (status == GSL_CONTINUE && ++iter < max_iter);

    /* output solution */
    for (i = 0; i < n; ++i)
      {
        double xi = (i + 1) * h;
        double u_exact = sin(M_PI * xi);
        double u_gsl = gsl_vector_get(u, i);

        printf("%f %.12e %.12e\n", xi, u_gsl, u_exact);
      }

    gsl_splinalg_itersolve_free(work);
  }

  gsl_spmatrix_free(A);
  gsl_spmatrix_free(C);
  gsl_vector_free(f);
  gsl_vector_free(u);

  return 0;
} /* main() */
```



## References and Further Reading

The implementation of the GMRES iterative solver closely follows the publications

- H. F. Walker, Implementation of the GMRES method using Householder transformations, SIAM J. Sci. Stat. Comput. 9(1), 1988.
- Y. Saad, Iterative methods for sparse linear systems, 2nd edition, SIAM, 2003.