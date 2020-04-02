# 多项式

This chapter describes functions for evaluating and solving polynomials. There are routines for finding real and complex roots of quadratic and cubic equations using analytic methods. An iterative polynomial solver is also available for finding the roots of general polynomials with real coefficients (of any order). The functions are declared in the header file `gsl_poly.h`.

本章描述用于多项式求值和求解的函数。含有用于通过解析方法寻找二次和三次方程的实根和复根的程序。对于寻找一个具有实系数(任何阶)的广义多项式的根也有一个迭代多项式求解器是可以使用的。函数被声明在头文件`gsl_poly.h`中。

## 多项式求值

The functions described here evaluate the polynomial

![P(x) = c[0] + c[1] x + c[2] x^2 + \dots + c[len-1] x^{len-1}](https://www.gnu.org/software/gsl/doc/html/_images/math/0cdb1361bfaba1c5ec1dc62b4c94e906455aa5e9.png)

using Horner’s method for stability. Inline versions of these functions are used when `HAVE_INLINE`is defined.

这里描述的函数对一下多项式进行求值

$P(x) = c[0]+c[1]x+c[2]x^2+...+c[len-1]x^{len-1}​$

- double `gsl_poly_eval`(const double *c[]*, const int *len*, const double *x*)

  This function evaluates a polynomial with real coefficients for the real variable `x`.

  这个函数对一个具有实系数的实变量`x`的多项式进行求值。

- gsl_complex `gsl_poly_complex_eval`(const double *c[]*, const int *len*, const gsl_complex *z*)

  This function evaluates a polynomial with real coefficients for the complex variable `z`.

  这个函数对一个具有实系数的复变量`z`的多项式是进行求值。

- gsl_complex `gsl_complex_poly_complex_eval`(const gsl_complex *c[]*, const int *len*, const gsl_complex *z*)

  This function evaluates a polynomial with complex coefficients for the complex variable `z`.

  这个函数返回一个具有复系数的复变量`z`的多项式进行求值。

- int `gsl_poly_eval_derivs`(const double *c[]*, const size_t *lenc*, const double *x*, double *res[]*, const size_t *lenres*)

  This function evaluates a polynomial and its derivatives storing the results in the array `res` of size `lenres`. The output array contains the values of ![d^k P(x)/d x^k](https://www.gnu.org/software/gsl/doc/html/_images/math/1600cca411950211f4183e112b84d963b9acf539.png) for the specified value of `x`starting with ![k = 0](https://www.gnu.org/software/gsl/doc/html/_images/math/7fe0629d90b771b38027bce3cadc1a36afe514bb.png).

  这个函数可以求一个多项式的值以及其导数，将其结果存储在具有尺寸`lenres`的数组`res`中。输出数组包含$$d^kP(x)/dx^k$$的对于特定的`x`的值的值，从$$k=0$$开始。



## 多项式的差商表示

The functions described here manipulate polynomials stored in Newton’s divided-difference representation. The use of divided-differences is described in Abramowitz & Stegun sections 25.1.4 and 25.2.26, and Burden and Faires, chapter 3, and discussed briefly below.

这里描述的函数操控以牛顿差商表示的多项式。差商的使用被描述在Abramowitz 和 Stegun的章节25.1.4 和 25.2.26以及Burden和Faires，第三章中，并会在下面得到简短描述。

Given a function ![f(x)](https://www.gnu.org/software/gsl/doc/html/_images/math/c7deb6ce5befe14135ebd23fa69801be7a796b15.png), an ![n](https://www.gnu.org/software/gsl/doc/html/_images/math/a24554f1502cf204e7aa24a6c064962c3504de48.png)th degree interpolating polynomial ![P_{n}(x)](https://www.gnu.org/software/gsl/doc/html/_images/math/3e2aa44dea8683f01996ed384aebd152865740c3.png) can be constructed which agrees with ![f](https://www.gnu.org/software/gsl/doc/html/_images/math/d947f897010fa03fc9854d6af00264296a50930c.png) at ![n+1](https://www.gnu.org/software/gsl/doc/html/_images/math/a3b5b5b94d164a7db425e53cf0d842a76e90d2da.png) distinct points ![x_0,x_1,...,x_{n}](https://www.gnu.org/software/gsl/doc/html/_images/math/a0aea8f09eb918ad6ecc6e3f35174842a37c1baa.png). This polynomial can be written in a form known as Newton’s divided-difference representation

![P_{n}(x) = f(x_0) + \sum_{k=1}^n [x_0,x_1,...,x_k] (x-x_0)(x-x_1) \cdots (x-x_{k-1})](https://www.gnu.org/software/gsl/doc/html/_images/math/044dd38880eb78dd29915b85f8d2a7a71a96a2d9.png)

给一个函数$$f(x)$$，一个$$n$$阶的插值多项式$$P_n(x)$$能被构建来在$$n+1$$个不同点$$x_0,x_1,...,x_n$$上与$$f$$相符。这个多项式能够以一个叫做牛顿差商表示的形式书写

$P_n(x) = f(x_0)+\displaystyle\sum_{k=1}^{n}{[x_0,x_1,...,x_k](x-x_0)(x-x_1)...(x-x_{k-1})}$

where the divided differences ![[x_0,x_1,...,x_k]](https://www.gnu.org/software/gsl/doc/html/_images/math/532aedc6dc9b02f5a48499a43080fe2d3a9c238c.png) are defined in section 25.1.4 of Abramowitz and Stegun. Additionally, it is possible to construct an interpolating polynomial of degree ![2n+1](https://www.gnu.org/software/gsl/doc/html/_images/math/c8f2db58f7120aae32ea289b45cc62c6a02130c5.png) which also matches the first derivatives of ![f](https://www.gnu.org/software/gsl/doc/html/_images/math/d947f897010fa03fc9854d6af00264296a50930c.png) at the points ![x_0,x_1,...,x_n](https://www.gnu.org/software/gsl/doc/html/_images/math/542e2ccb3ff98872f2b9730336305ae2bbe33126.png). This is called the Hermite interpolating polynomial and is defined as

![H_{2n+1}(x) = f(z_0) + \sum_{k=1}^{2n+1} [z_0,z_1,...,z_k] (x-z_0)(x-z_1) \cdots (x-z_{k-1})](https://www.gnu.org/software/gsl/doc/html/_images/math/59c62018c8806542b7b55d3503f53c18b4000966.png)

其中差商$$[x_0,x_1,...,x_k]$$被定义在Abramowitz和Stegun的25.1.4小节。此外，

where the elements of ![z = \{x_0,x_0,x_1,x_1,...,x_n,x_n\}](https://www.gnu.org/software/gsl/doc/html/_images/math/5de7f0979bcfbf5952f8abfdc09982a157961996.png) are defined by ![z_{2k} = z_{2k+1} = x_k](https://www.gnu.org/software/gsl/doc/html/_images/math/562b31e7bc6e44eba8aebee039c0e8866cfa5c2b.png). The divided-differences ![[z_0,z_1,...,z_k]](https://www.gnu.org/software/gsl/doc/html/_images/math/c18706b2d8130826f83c8bc4320be0ddbb083cfa.png) are discussed in Burden and Faires, section 3.4.

- int `gsl_poly_dd_init`(double *dd[]*, const double *xa[]*, const double *ya[]*, size_t *size*)

  This function computes a divided-difference representation of the interpolating polynomial for the points ![(x, y)](https://www.gnu.org/software/gsl/doc/html/_images/math/1e0f6e79f39b16046a35b11f9bf35a25cc49d5e1.png) stored in the arrays `xa` and `ya` of length `size`. On output the divided-differences of (`xa`, `ya`) are stored in the array `dd`, also of length `size`. Using the notation above, ![dd[k] = [x_0,x_1,...,x_k]](https://www.gnu.org/software/gsl/doc/html/_images/math/7dae05058ef54ed5653c682f30e6fe393bfa5536.png).

- double `gsl_poly_dd_eval`(const double *dd[]*, const double *xa[]*, const size_t *size*, const double *x*)

  This function evaluates the polynomial stored in divided-difference form in the arrays `dd` and `xa` of length `size` at the point `x`. An inline version of this function is used when `HAVE_INLINE` is defined.

- int `gsl_poly_dd_taylor`(double *c[]*, double *xp*, const double *dd[]*, const double *xa[]*, size_t *size*, double *w[]*)

  This function converts the divided-difference representation of a polynomial to a Taylor expansion. The divided-difference representation is supplied in the arrays `dd` and `xa` of length `size`. On output the Taylor coefficients of the polynomial expanded about the point `xp` are stored in the array `c` also of length `size`. A workspace of length `size` must be provided in the array `w`.

- int `gsl_poly_dd_hermite_init`(double *dd[]*, double *za[]*, const double *xa[]*, const double *ya[]*, const double *dya[]*, const size_t *size*)

  This function computes a divided-difference representation of the interpolating Hermite polynomial for the points ![(x,y)](https://www.gnu.org/software/gsl/doc/html/_images/math/7ecf269c4c3ea45e59922da2237627de837dbc5f.png) stored in the arrays `xa` and `ya` of length `size`. Hermite interpolation constructs polynomials which also match first derivatives ![dy/dx](https://www.gnu.org/software/gsl/doc/html/_images/math/c3a60a8e3f1873c7c661987ff5f6c9ad4a6b602e.png) which are provided in the array `dya` also of length `size`. The first derivatives can be incorported into the usual divided-difference algorithm by forming a new dataset ![z = \{x_0,x_0,x_1,x_1,...\}](https://www.gnu.org/software/gsl/doc/html/_images/math/a0999f04689e758dc75493150f951e07dc17d57e.png), which is stored in the array `za` of length 2*`size` on output. On output the divided-differences of the Hermite representation are stored in the array `dd`, also of length 2*`size`. Using the notation above, ![dd[k] = [z_0,z_1,...,z_k]](https://www.gnu.org/software/gsl/doc/html/_images/math/c718eaa9e3ec5ff290188f4c55f2e6090c49ab7d.png). The resulting Hermite polynomial can be evaluated by calling [`gsl_poly_dd_eval()`](https://www.gnu.org/software/gsl/doc/html/poly.html#c.gsl_poly_dd_eval) and using `za` for the input argument `xa`.



## Quadratic Equations

- int `gsl_poly_solve_quadratic`(double *a*, double *b*, double *c*, double * *x0*, double * *x1*)

  This function finds the real roots of the quadratic equation,![a x^2 + b x + c = 0](https://www.gnu.org/software/gsl/doc/html/_images/math/ab20e2443abdf7561a8025b512f9cedc11d7640d.png)The number of real roots (either zero, one or two) is returned, and their locations are stored in `x0` and `x1`. If no real roots are found then `x0` and `x1` are not modified. If one real root is found (i.e. if ![a=0](https://www.gnu.org/software/gsl/doc/html/_images/math/a63055f5b7d9ff061a6bb0f3b12347110e3fa06c.png)) then it is stored in `x0`. When two real roots are found they are stored in `x0` and `x1` in ascending order. The case of coincident roots is not considered special. For example ![(x-1)^2=0](https://www.gnu.org/software/gsl/doc/html/_images/math/7faf84d72ec62ea035312ec4e1d859d736bc566c.png) will have two roots, which happen to have exactly equal values.The number of roots found depends on the sign of the discriminant ![b^2 - 4 a c](https://www.gnu.org/software/gsl/doc/html/_images/math/2433939147739485b5629b7bcc077a1849169723.png). This will be subject to rounding and cancellation errors when computed in double precision, and will also be subject to errors if the coefficients of the polynomial are inexact. These errors may cause a discrete change in the number of roots. However, for polynomials with small integer coefficients the discriminant can always be computed exactly.

- int `gsl_poly_complex_solve_quadratic`(double *a*, double *b*, double *c*, gsl_complex * *z0*, gsl_complex * *z1*)

  This function finds the complex roots of the quadratic equation,![a z^2 + b z + c = 0](https://www.gnu.org/software/gsl/doc/html/_images/math/115ad1b884c5bcf9d404428618a1e6f043aa8a3b.png)The number of complex roots is returned (either one or two) and the locations of the roots are stored in `z0` and `z1`. The roots are returned in ascending order, sorted first by their real components and then by their imaginary components. If only one real root is found (i.e. if ![a=0](https://www.gnu.org/software/gsl/doc/html/_images/math/a63055f5b7d9ff061a6bb0f3b12347110e3fa06c.png)) then it is stored in `z0`.



## Cubic Equations

- int `gsl_poly_solve_cubic`(double *a*, double *b*, double *c*, double * *x0*, double * *x1*, double * *x2*)

  This function finds the real roots of the cubic equation,![x^3 + a x^2 + b x + c = 0](https://www.gnu.org/software/gsl/doc/html/_images/math/9908e1179b552b30892e342462b26ba7d60ffd35.png)with a leading coefficient of unity. The number of real roots (either one or three) is returned, and their locations are stored in `x0`, `x1` and `x2`. If one real root is found then only `x0` is modified. When three real roots are found they are stored in `x0`, `x1` and `x2` in ascending order. The case of coincident roots is not considered special. For example, the equation ![(x-1)^3=0](https://www.gnu.org/software/gsl/doc/html/_images/math/9ae8cde80ed4dc961a2b75f0add02ac69b91f718.png) will have three roots with exactly equal values. As in the quadratic case, finite precision may cause equal or closely-spaced real roots to move off the real axis into the complex plane, leading to a discrete change in the number of real roots.

- int `gsl_poly_complex_solve_cubic`(double *a*, double *b*, double *c*, gsl_complex * *z0*, gsl_complex * *z1*, gsl_complex * *z2*)

  This function finds the complex roots of the cubic equation,![z^3 + a z^2 + b z + c = 0](https://www.gnu.org/software/gsl/doc/html/_images/math/72a0cafd8c4600249253fcc99f9896ecb6db402e.png)The number of complex roots is returned (always three) and the locations of the roots are stored in `z0`, `z1` and `z2`. The roots are returned in ascending order, sorted first by their real components and then by their imaginary components.



## General Polynomial Equations

The roots of polynomial equations cannot be found analytically beyond the special cases of the quadratic, cubic and quartic equation. The algorithm described in this section uses an iterative method to find the approximate locations of roots of higher order polynomials.

- `gsl_poly_complex_workspace`

  This workspace contains parameters used for finding roots of general polynomials

- [gsl_poly_complex_workspace](https://www.gnu.org/software/gsl/doc/html/poly.html#c.gsl_poly_complex_workspace) * `gsl_poly_complex_workspace_alloc`(size_t *n*)

  This function allocates space for a [`gsl_poly_complex_workspace`](https://www.gnu.org/software/gsl/doc/html/poly.html#c.gsl_poly_complex_workspace) struct and a workspace suitable for solving a polynomial with `n` coefficients using the routine [`gsl_poly_complex_solve()`](https://www.gnu.org/software/gsl/doc/html/poly.html#c.gsl_poly_complex_solve).The function returns a pointer to the newly allocated [`gsl_poly_complex_workspace`](https://www.gnu.org/software/gsl/doc/html/poly.html#c.gsl_poly_complex_workspace) if no errors were detected, and a null pointer in the case of error.

- void `gsl_poly_complex_workspace_free`([gsl_poly_complex_workspace](https://www.gnu.org/software/gsl/doc/html/poly.html#c.gsl_poly_complex_workspace) * *w*)

  This function frees all the memory associated with the workspace `w`.

- int `gsl_poly_complex_solve`(const double * *a*, size_t *n*, [gsl_poly_complex_workspace](https://www.gnu.org/software/gsl/doc/html/poly.html#c.gsl_poly_complex_workspace) * *w*, gsl_complex_packed_ptr *z*)

  This function computes the roots of the general polynomial![P(x) = a_0 + a_1 x + a_2 x^2 + \cdots + a_{n-1} x^{n-1}](https://www.gnu.org/software/gsl/doc/html/_images/math/c5a7c66ab1442eb9b76202258cee42b7ca6305bc.png)using balanced-QR reduction of the companion matrix. The parameter `n` specifies the length of the coefficient array. The coefficient of the highest order term must be non-zero. The function requires a workspace `w` of the appropriate size. The ![n-1](https://www.gnu.org/software/gsl/doc/html/_images/math/0acbb742031c41a269215d223bd0d699a0cc0522.png) roots are returned in the packed complex array `z` of length ![2(n-1)](https://www.gnu.org/software/gsl/doc/html/_images/math/94ce50ed55100c028b83d1218a0728186764bf39.png), alternating real and imaginary parts.The function returns `GSL_SUCCESS` if all the roots are found. If the QR reduction does not converge, the error handler is invoked with an error code of `GSL_EFAILED`. Note that due to finite precision, roots of higher multiplicity are returned as a cluster of simple roots with reduced accuracy. The solution of polynomials with higher-order roots requires specialized algorithms that take the multiplicity structure into account (see e.g. Z. Zeng, Algorithm 835, ACM Transactions on Mathematical Software, Volume 30, Issue 2 (2004), pp 218–236).

## Examples

To demonstrate the use of the general polynomial solver we will take the polynomial ![P(x) = x^5 - 1](https://www.gnu.org/software/gsl/doc/html/_images/math/3139de50d2dcf290b00c02d58315f6cc1a4dd22f.png) which has these roots:

![1, e^{2\pi i / 5}, e^{4\pi i / 5}, e^{6\pi i / 5}, e^{8\pi i / 5}](https://www.gnu.org/software/gsl/doc/html/_images/math/bd1dfdc8b747dd009a50889dd4c04d7eecd2cbd0.png)

The following program will find these roots.

```
#include <stdio.h>
#include <gsl/gsl_poly.h>

int
main (void)
{
  int i;
  /* coefficients of P(x) =  -1 + x^5  */
  double a[6] = { -1, 0, 0, 0, 0, 1 };
  double z[10];

  gsl_poly_complex_workspace * w
      = gsl_poly_complex_workspace_alloc (6);

  gsl_poly_complex_solve (a, 6, w, z);

  gsl_poly_complex_workspace_free (w);

  for (i = 0; i < 5; i++)
    {
      printf ("z%d = %+.18f %+.18f\n",
              i, z[2*i], z[2*i+1]);
    }

  return 0;
}
```

The output of the program is

```
z0 = -0.809016994374947673 +0.587785252292473359
z1 = -0.809016994374947673 -0.587785252292473359
z2 = +0.309016994374947507 +0.951056516295152976
z3 = +0.309016994374947507 -0.951056516295152976
z4 = +0.999999999999999889 +0.000000000000000000
```

which agrees with the analytic result, ![z_n = \exp(2 \pi n i/5)](https://www.gnu.org/software/gsl/doc/html/_images/math/67c722396ef61b3f886770225579f23fa7da6799.png).

## References and Further Reading

The balanced-QR method and its error analysis are described in the following papers,

- R.S. Martin, G. Peters and J.H. Wilkinson, “The QR Algorithm for Real Hessenberg Matrices”, Numerische Mathematik, 14 (1970), 219–231.
- B.N. Parlett and C. Reinsch, “Balancing a Matrix for Calculation of Eigenvalues and Eigenvectors”, Numerische Mathematik, 13 (1969), 293–304.
- A. Edelman and H. Murakami, “Polynomial roots from companion matrix eigenvalues”, Mathematics of Computation, Vol.: 64, No.: 210 (1995), 763–776.

The formulas for divided differences are given in the following texts,

- Abramowitz and Stegun, Handbook of Mathematical Functions, Sections 25.1.4 and 25.2.26.
- R. L. Burden and J. D. Faires, Numerical Analysis, 9th edition, ISBN 0-538-73351-9, 2011.