# 特殊函数

This chapter describes the GSL special function library. The library includes routines for calculating the values of Airy functions, Bessel functions, Clausen functions, Coulomb wave functions, Coupling coefficients, the Dawson function, Debye functions, Dilogarithms, Elliptic integrals, Jacobi elliptic functions, Error functions, Exponential integrals, Fermi-Dirac functions, Gamma functions, Gegenbauer functions, Hermite polynomials and functions, Hypergeometric functions, Laguerre functions, Legendre functions and Spherical Harmonics, the Psi (Digamma) Function, Synchrotron functions, Transport functions, Trigonometric functions and Zeta functions. Each routine also computes an estimate of the numerical error in the calculated value of the function.

本章面搜狐GSL特殊函数库。本库包含计算艾里函数，贝塞尔函数，克劳森函数，库伦波函数，耦合系数，道森函数，德拜函数，二重对数，椭圆积分，雅可比椭圆函数，误差函数，指数积分，费米狄拉克函数，伽马函数，盖根堡函数，厄米特多项式和函数，超几何函数，拉盖尔函数，勒让德函数和球谐函数，Psi(双伽马)函数，同步加速函数，输运函数，三角函数和ζ函数的值的程序。每个程序也计算函数所计算的值的数值误差的一个估值。

The functions in this chapter are declared in individual header files, such as `gsl_sf_airy.h`, `gsl_sf_bessel.h`, etc. The complete set of header files can be included using the file `gsl_sf.h`.

本章描述的函数被声明在单独的头文件，比如`gsl_sf_airy.h`，`gsl_sf_bessel.h`等。头文件的完全集能通过文件`gsl_sf.h`进行包含。

## 用法

The special functions are available in two calling conventions, a *natural form* which returns the numerical value of the function and an *error-handling form* which returns an error code. The two types of function provide alternative ways of accessing the same underlying code.

特殊函数可通过两个调用约定进行获取，一个是*自然形式*其返回函数的数值，另一个是*误差处理形式*其返回一个误差代码。两个类型的函数提供了访问同样的底层代码的可选方式。

The *natural form* returns only the value of the function and can be used directly in mathematical expressions. For example, the following function call will compute the value of the Bessel function ![J_0(x)](https://www.gnu.org/software/gsl/doc/html/_images/math/436b1799db46b1ecb38cc2fde801e882814c7119.png):

```
double y = gsl_sf_bessel_J0 (x);
```

There is no way to access an error code or to estimate the error using this method. To allow access to this information the alternative error-handling form stores the value and error in a modifiable argument:

```
gsl_sf_result result;
int status = gsl_sf_bessel_J0_e (x, &result);
```

The error-handling functions have the suffix `_e`. The returned status value indicates error conditions such as overflow, underflow or loss of precision. If there are no errors the error-handling functions return `GSL_SUCCESS`.

## The gsl_sf_result struct

The error handling form of the special functions always calculate an error estimate along with the value of the result. Therefore, structures are provided for amalgamating a value and error estimate. These structures are declared in the header file `gsl_sf_result.h`.

The following struct contains value and error fields.

- `gsl_sf_result`

  `typedef struct {   double val;   double err; } gsl_sf_result; `The field `val` contains the value and the field `err` contains an estimate of the absolute error in the value.

In some cases, an overflow or underflow can be detected and handled by a function. In this case, it may be possible to return a scaling exponent as well as an error/value pair in order to save the result from exceeding the dynamic range of the built-in types. The following struct contains value and error fields as well as an exponent field such that the actual result is obtained as`result * 10^(e10)`.

- `gsl_sf_result_e10`

  `typedef struct {   double val;   double err;   int    e10; } gsl_sf_result_e10; `

## Modes

The goal of the library is to achieve double precision accuracy wherever possible. However the cost of evaluating some special functions to double precision can be significant, particularly where very high order terms are required. In these cases a `mode` argument, of type [`gsl_mode_t`](https://www.gnu.org/software/gsl/doc/html/specfunc.html#c.gsl_mode_t) allows the accuracy of the function to be reduced in order to improve performance. The following precision levels are available for the mode argument,

- `gsl_mode_t`

  `GSL_PREC_DOUBLE`Double-precision, a relative accuracy of approximately ![2 * 10^{-16}](https://www.gnu.org/software/gsl/doc/html/_images/math/21d63c4c48ddc2bde43e56c548478d05334402e5.png).`GSL_PREC_SINGLE`Single-precision, a relative accuracy of approximately ![10^{-7}](https://www.gnu.org/software/gsl/doc/html/_images/math/6c6f9131cb760e7552278fb3ce3d0adf0c44d6a4.png).`GSL_PREC_APPROX`Approximate values, a relative accuracy of approximately ![5 * 10^{-4}](https://www.gnu.org/software/gsl/doc/html/_images/math/b504f1f1d8160dd46e29db41cc61bdf5e8b3fce5.png).

The approximate mode provides the fastest evaluation at the lowest accuracy.

## Airy Functions and Derivatives





The Airy functions ![Ai(x)](https://www.gnu.org/software/gsl/doc/html/_images/math/57e20afe5e5473ee33592444a8240297320a1d79.png) and ![Bi(x)](https://www.gnu.org/software/gsl/doc/html/_images/math/0dc53af3809b4351f9a56b7e1fc932d0a3f4ed58.png) are defined by the integral representations,

![Ai(x) & = {1\over\pi} \int_0^\infty \cos(t^3/3 + xt ) \,dt \\ Bi(x) & = {1\over\pi} \int_0^\infty (e^{-t^3/3 + xt} + \sin(t^3/3 + xt)) \,dt](https://www.gnu.org/software/gsl/doc/html/_images/math/b3a74f498aa0fe283aed6cf51a15c6bb75a27e58.png)

For further information see Abramowitz & Stegun, Section 10.4. The Airy functions are defined in the header file `gsl_sf_airy.h`.

### Airy Functions

- double `gsl_sf_airy_Ai`(double *x*, [gsl_mode_t](https://www.gnu.org/software/gsl/doc/html/specfunc.html#c.gsl_mode_t) *mode*)

- int `gsl_sf_airy_Ai_e`(double *x*, [gsl_mode_t](https://www.gnu.org/software/gsl/doc/html/specfunc.html#c.gsl_mode_t) *mode*, [gsl_sf_result](https://www.gnu.org/software/gsl/doc/html/specfunc.html#c.gsl_sf_result) * *result*)

  These routines compute the Airy function ![Ai(x)](https://www.gnu.org/software/gsl/doc/html/_images/math/57e20afe5e5473ee33592444a8240297320a1d79.png) with an accuracy specified by [`mode`](https://www.gnu.org/software/gsl/doc/html/montecarlo.html#c.mode).

- double `gsl_sf_airy_Bi`(double *x*, [gsl_mode_t](https://www.gnu.org/software/gsl/doc/html/specfunc.html#c.gsl_mode_t) *mode*)

- int `gsl_sf_airy_Bi_e`(double *x*, [gsl_mode_t](https://www.gnu.org/software/gsl/doc/html/specfunc.html#c.gsl_mode_t) *mode*, [gsl_sf_result](https://www.gnu.org/software/gsl/doc/html/specfunc.html#c.gsl_sf_result) * *result*)

  These routines compute the Airy function ![Bi(x)](https://www.gnu.org/software/gsl/doc/html/_images/math/0dc53af3809b4351f9a56b7e1fc932d0a3f4ed58.png) with an accuracy specified by [`mode`](https://www.gnu.org/software/gsl/doc/html/montecarlo.html#c.mode).

- double `gsl_sf_airy_Ai_scaled`(double *x*, [gsl_mode_t](https://www.gnu.org/software/gsl/doc/html/specfunc.html#c.gsl_mode_t) *mode*)

- int `gsl_sf_airy_Ai_scaled_e`(double *x*, [gsl_mode_t](https://www.gnu.org/software/gsl/doc/html/specfunc.html#c.gsl_mode_t) *mode*, [gsl_sf_result](https://www.gnu.org/software/gsl/doc/html/specfunc.html#c.gsl_sf_result) * *result*)

  These routines compute a scaled version of the Airy function ![S_A(x) Ai(x)](https://www.gnu.org/software/gsl/doc/html/_images/math/1d31717e8cc910ec5dd5a06689ee631fdf5fbb8d.png). For ![x > 0](https://www.gnu.org/software/gsl/doc/html/_images/math/fe20eedd0aef40a657ad04835ff74fa8a73efb56.png) the scaling factor ![S_A(x)](https://www.gnu.org/software/gsl/doc/html/_images/math/7c054a5d1a71380023b979551f9ec0a46709710e.png) is ![\exp(+(2/3) x^{3/2})](https://www.gnu.org/software/gsl/doc/html/_images/math/c3bb32a20dab5bba222708b012c0304280f288b4.png), and is 1 for ![x < 0](https://www.gnu.org/software/gsl/doc/html/_images/math/93f2a3ce4b90fb07a7dc88c1bc60a6a43c93f45b.png).

- double `gsl_sf_airy_Bi_scaled`(double *x*, [gsl_mode_t](https://www.gnu.org/software/gsl/doc/html/specfunc.html#c.gsl_mode_t) *mode*)

- int `gsl_sf_airy_Bi_scaled_e`(double *x*, [gsl_mode_t](https://www.gnu.org/software/gsl/doc/html/specfunc.html#c.gsl_mode_t) *mode*, [gsl_sf_result](https://www.gnu.org/software/gsl/doc/html/specfunc.html#c.gsl_sf_result) * *result*)

  These routines compute a scaled version of the Airy function ![S_B(x) Bi(x)](https://www.gnu.org/software/gsl/doc/html/_images/math/870d1ce611612e9a2db950960617806b6b3b7a90.png). For ![x > 0](https://www.gnu.org/software/gsl/doc/html/_images/math/fe20eedd0aef40a657ad04835ff74fa8a73efb56.png) the scaling factor ![S_B(x)](https://www.gnu.org/software/gsl/doc/html/_images/math/506492848196327526b4ca1979ca187f118b8eab.png) is ![exp(-(2/3) x^{3/2})](https://www.gnu.org/software/gsl/doc/html/_images/math/839856b009b5682b75b50247bb88d65f99e2748e.png), and is 1 for ![x < 0](https://www.gnu.org/software/gsl/doc/html/_images/math/93f2a3ce4b90fb07a7dc88c1bc60a6a43c93f45b.png).

### Derivatives of Airy Functions

- double `gsl_sf_airy_Ai_deriv`(double *x*, [gsl_mode_t](https://www.gnu.org/software/gsl/doc/html/specfunc.html#c.gsl_mode_t) *mode*)

- int `gsl_sf_airy_Ai_deriv_e`(double *x*, [gsl_mode_t](https://www.gnu.org/software/gsl/doc/html/specfunc.html#c.gsl_mode_t) *mode*, [gsl_sf_result](https://www.gnu.org/software/gsl/doc/html/specfunc.html#c.gsl_sf_result) * *result*)

  These routines compute the Airy function derivative ![Ai'(x)](https://www.gnu.org/software/gsl/doc/html/_images/math/afbb77d674618893664e0e97712082a2f22d9689.png) with an accuracy specified by [`mode`](https://www.gnu.org/software/gsl/doc/html/montecarlo.html#c.mode).

- double `gsl_sf_airy_Bi_deriv`(double *x*, [gsl_mode_t](https://www.gnu.org/software/gsl/doc/html/specfunc.html#c.gsl_mode_t) *mode*)

- int `gsl_sf_airy_Bi_deriv_e`(double *x*, [gsl_mode_t](https://www.gnu.org/software/gsl/doc/html/specfunc.html#c.gsl_mode_t) *mode*, [gsl_sf_result](https://www.gnu.org/software/gsl/doc/html/specfunc.html#c.gsl_sf_result) * *result*)

  These routines compute the Airy function derivative ![Bi'(x)](https://www.gnu.org/software/gsl/doc/html/_images/math/d365d807fc54eaadc4a703f503cf251405f56ad0.png) with an accuracy specified by [`mode`](https://www.gnu.org/software/gsl/doc/html/montecarlo.html#c.mode).

- double `gsl_sf_airy_Ai_deriv_scaled`(double *x*, [gsl_mode_t](https://www.gnu.org/software/gsl/doc/html/specfunc.html#c.gsl_mode_t) *mode*)

- int `gsl_sf_airy_Ai_deriv_scaled_e`(double *x*, [gsl_mode_t](https://www.gnu.org/software/gsl/doc/html/specfunc.html#c.gsl_mode_t) *mode*, [gsl_sf_result](https://www.gnu.org/software/gsl/doc/html/specfunc.html#c.gsl_sf_result) * *result*)

  These routines compute the scaled Airy function derivative ![S_A(x) Ai'(x)](https://www.gnu.org/software/gsl/doc/html/_images/math/28d77b600eeafb0979f5a07613ee5ea1494f8784.png). For ![x > 0](https://www.gnu.org/software/gsl/doc/html/_images/math/fe20eedd0aef40a657ad04835ff74fa8a73efb56.png) the scaling factor ![S_A(x)](https://www.gnu.org/software/gsl/doc/html/_images/math/7c054a5d1a71380023b979551f9ec0a46709710e.png) is ![\exp(+(2/3) x^{3/2})](https://www.gnu.org/software/gsl/doc/html/_images/math/c3bb32a20dab5bba222708b012c0304280f288b4.png), and is 1 for ![x < 0](https://www.gnu.org/software/gsl/doc/html/_images/math/93f2a3ce4b90fb07a7dc88c1bc60a6a43c93f45b.png).

- double `gsl_sf_airy_Bi_deriv_scaled`(double *x*, [gsl_mode_t](https://www.gnu.org/software/gsl/doc/html/specfunc.html#c.gsl_mode_t) *mode*)

- int `gsl_sf_airy_Bi_deriv_scaled_e`(double *x*, [gsl_mode_t](https://www.gnu.org/software/gsl/doc/html/specfunc.html#c.gsl_mode_t) *mode*, [gsl_sf_result](https://www.gnu.org/software/gsl/doc/html/specfunc.html#c.gsl_sf_result) * *result*)

  These routines compute the scaled Airy function derivative ![S_B(x) Bi'(x)](https://www.gnu.org/software/gsl/doc/html/_images/math/c80fe8e3a3ffdc181cb90c1b61c9d62b4884dbb0.png). For ![x > 0](https://www.gnu.org/software/gsl/doc/html/_images/math/fe20eedd0aef40a657ad04835ff74fa8a73efb56.png) the scaling factor ![S_B(x)](https://www.gnu.org/software/gsl/doc/html/_images/math/506492848196327526b4ca1979ca187f118b8eab.png) is ![exp(-(2/3) x^{3/2})](https://www.gnu.org/software/gsl/doc/html/_images/math/839856b009b5682b75b50247bb88d65f99e2748e.png), and is 1 for ![x < 0](https://www.gnu.org/software/gsl/doc/html/_images/math/93f2a3ce4b90fb07a7dc88c1bc60a6a43c93f45b.png).

### Zeros of Airy Functions

- double `gsl_sf_airy_zero_Ai`(unsigned int *s*)

- int `gsl_sf_airy_zero_Ai_e`(unsigned int *s*, [gsl_sf_result](https://www.gnu.org/software/gsl/doc/html/specfunc.html#c.gsl_sf_result) * *result*)

  These routines compute the location of the `s`-th zero of the Airy function ![Ai(x)](https://www.gnu.org/software/gsl/doc/html/_images/math/57e20afe5e5473ee33592444a8240297320a1d79.png).

- double `gsl_sf_airy_zero_Bi`(unsigned int *s*)

- int `gsl_sf_airy_zero_Bi_e`(unsigned int *s*, [gsl_sf_result](https://www.gnu.org/software/gsl/doc/html/specfunc.html#c.gsl_sf_result) * *result*)

  These routines compute the location of the `s`-th zero of the Airy function ![Bi(x)](https://www.gnu.org/software/gsl/doc/html/_images/math/0dc53af3809b4351f9a56b7e1fc932d0a3f4ed58.png).

### Zeros of Derivatives of Airy Functions

- double `gsl_sf_airy_zero_Ai_deriv`(unsigned int *s*)

- int `gsl_sf_airy_zero_Ai_deriv_e`(unsigned int *s*, [gsl_sf_result](https://www.gnu.org/software/gsl/doc/html/specfunc.html#c.gsl_sf_result) * *result*)

  These routines compute the location of the `s`-th zero of the Airy function derivative ![Ai'(x)](https://www.gnu.org/software/gsl/doc/html/_images/math/afbb77d674618893664e0e97712082a2f22d9689.png).

- double `gsl_sf_airy_zero_Bi_deriv`(unsigned int *s*)

- int `gsl_sf_airy_zero_Bi_deriv_e`(unsigned int *s*, [gsl_sf_result](https://www.gnu.org/software/gsl/doc/html/specfunc.html#c.gsl_sf_result) * *result*)

  These routines compute the location of the `s`-th zero of the Airy function derivative ![Bi'(x)](https://www.gnu.org/software/gsl/doc/html/_images/math/d365d807fc54eaadc4a703f503cf251405f56ad0.png).

## Bessel Functions

The routines described in this section compute the Cylindrical Bessel functions ![J_n(x)](https://www.gnu.org/software/gsl/doc/html/_images/math/315c7555aa864ab7b76601b1714cd9c2e66883b3.png), ![Y_n(x)](https://www.gnu.org/software/gsl/doc/html/_images/math/d68f2b9626c69cf941c0281eb1296e9a7a55d6bd.png), Modified cylindrical Bessel functions ![I_n(x)](https://www.gnu.org/software/gsl/doc/html/_images/math/7e78a3c2669b08dafec3eee5c5339e16063277b3.png), ![K_n(x)](https://www.gnu.org/software/gsl/doc/html/_images/math/f44f0a148e38562c9d66c1d7ef29a270cb8a3e75.png), Spherical Bessel functions ![j_l(x)](https://www.gnu.org/software/gsl/doc/html/_images/math/2c58211abd5178158a7ccf83ba0d0ded869cb455.png), ![y_l(x)](https://www.gnu.org/software/gsl/doc/html/_images/math/a5237e6a66274662bef5fadebb08539203de721d.png), and Modified Spherical Bessel functions ![i_l(x)](https://www.gnu.org/software/gsl/doc/html/_images/math/66569b71033cd7295c2de3f85e2b52d48cb9c035.png), ![k_l(x)](https://www.gnu.org/software/gsl/doc/html/_images/math/0cae625b70bafb5ce956c6ddf2af444cdbdad4cd.png). For more information see Abramowitz & Stegun, Chapters 9 and 10. The Bessel functions are defined in the header file `gsl_sf_bessel.h`.

### Regular Cylindrical Bessel Functions







- double `gsl_sf_bessel_J0`(double *x*)

- int `gsl_sf_bessel_J0_e`(double *x*, [gsl_sf_result](https://www.gnu.org/software/gsl/doc/html/specfunc.html#c.gsl_sf_result) * *result*)

  These routines compute the regular cylindrical Bessel function of zeroth order, ![J_0(x)](https://www.gnu.org/software/gsl/doc/html/_images/math/436b1799db46b1ecb38cc2fde801e882814c7119.png).

- double `gsl_sf_bessel_J1`(double *x*)

- int `gsl_sf_bessel_J1_e`(double *x*, [gsl_sf_result](https://www.gnu.org/software/gsl/doc/html/specfunc.html#c.gsl_sf_result) * *result*)

  These routines compute the regular cylindrical Bessel function of first order, ![J_1(x)](https://www.gnu.org/software/gsl/doc/html/_images/math/cba66ff642120cf89928250899e1bf9a8dd258c2.png).

- double `gsl_sf_bessel_Jn`(int *n*, double *x*)

- int `gsl_sf_bessel_Jn_e`(int *n*, double *x*, [gsl_sf_result](https://www.gnu.org/software/gsl/doc/html/specfunc.html#c.gsl_sf_result) * *result*)

  These routines compute the regular cylindrical Bessel function of order `n`, ![J_n(x)](https://www.gnu.org/software/gsl/doc/html/_images/math/315c7555aa864ab7b76601b1714cd9c2e66883b3.png).

- int `gsl_sf_bessel_Jn_array`(int *nmin*, int *nmax*, double *x*, double *result_array[]*)

  This routine computes the values of the regular cylindrical Bessel functions ![J_n(x)](https://www.gnu.org/software/gsl/doc/html/_images/math/315c7555aa864ab7b76601b1714cd9c2e66883b3.png) for ![n](https://www.gnu.org/software/gsl/doc/html/_images/math/a24554f1502cf204e7aa24a6c064962c3504de48.png) from `nmin` to `nmax` inclusive, storing the results in the array `result_array`. The values are computed using recurrence relations for efficiency, and therefore may differ slightly from the exact values.

### Irregular Cylindrical Bessel Functions





- double `gsl_sf_bessel_Y0`(double *x*)

- int `gsl_sf_bessel_Y0_e`(double *x*, [gsl_sf_result](https://www.gnu.org/software/gsl/doc/html/specfunc.html#c.gsl_sf_result) * *result*)

  These routines compute the irregular cylindrical Bessel function of zeroth order, ![Y_0(x)](https://www.gnu.org/software/gsl/doc/html/_images/math/29ab2442fcaabad16038e2a3c2761c042d3b8f1e.png), for ![x>0](https://www.gnu.org/software/gsl/doc/html/_images/math/937b5840b2fb62e72aff4740c9d3c14db3b06720.png).

- double `gsl_sf_bessel_Y1`(double *x*)

- int `gsl_sf_bessel_Y1_e`(double *x*, [gsl_sf_result](https://www.gnu.org/software/gsl/doc/html/specfunc.html#c.gsl_sf_result) * *result*)

  These routines compute the irregular cylindrical Bessel function of first order, ![Y_1(x)](https://www.gnu.org/software/gsl/doc/html/_images/math/3c8a8846c5b5652a2699bf5e6fa106cf2be6ce15.png), for ![x>0](https://www.gnu.org/software/gsl/doc/html/_images/math/937b5840b2fb62e72aff4740c9d3c14db3b06720.png).

- double `gsl_sf_bessel_Yn`(int *n*, double *x*)

- int `gsl_sf_bessel_Yn_e`(int *n*, double *x*, [gsl_sf_result](https://www.gnu.org/software/gsl/doc/html/specfunc.html#c.gsl_sf_result) * *result*)

  These routines compute the irregular cylindrical Bessel function of order `n`, ![Y_n(x)](https://www.gnu.org/software/gsl/doc/html/_images/math/d68f2b9626c69cf941c0281eb1296e9a7a55d6bd.png), for ![x>0](https://www.gnu.org/software/gsl/doc/html/_images/math/937b5840b2fb62e72aff4740c9d3c14db3b06720.png).

- int `gsl_sf_bessel_Yn_array`(int *nmin*, int *nmax*, double *x*, double *result_array[]*)

  This routine computes the values of the irregular cylindrical Bessel functions ![Y_n(x)](https://www.gnu.org/software/gsl/doc/html/_images/math/d68f2b9626c69cf941c0281eb1296e9a7a55d6bd.png) for ![n](https://www.gnu.org/software/gsl/doc/html/_images/math/a24554f1502cf204e7aa24a6c064962c3504de48.png) from `nmin` to `nmax` inclusive, storing the results in the array `result_array`. The domain of the function is ![x>0](https://www.gnu.org/software/gsl/doc/html/_images/math/937b5840b2fb62e72aff4740c9d3c14db3b06720.png). The values are computed using recurrence relations for efficiency, and therefore may differ slightly from the exact values.

### Regular Modified Cylindrical Bessel Functions







- double `gsl_sf_bessel_I0`(double *x*)

- int `gsl_sf_bessel_I0_e`(double *x*, [gsl_sf_result](https://www.gnu.org/software/gsl/doc/html/specfunc.html#c.gsl_sf_result) * *result*)

  These routines compute the regular modified cylindrical Bessel function of zeroth order, ![I_0(x)](https://www.gnu.org/software/gsl/doc/html/_images/math/80b44c6d78fe7fb9a233c4887204aa666cd4fb2d.png).

- double `gsl_sf_bessel_I1`(double *x*)

- int `gsl_sf_bessel_I1_e`(double *x*, [gsl_sf_result](https://www.gnu.org/software/gsl/doc/html/specfunc.html#c.gsl_sf_result) * *result*)

  These routines compute the regular modified cylindrical Bessel function of first order, ![I_1(x)](https://www.gnu.org/software/gsl/doc/html/_images/math/64d3c5b2a2945edfdce0ce46efc26168c40c6a41.png).

- double `gsl_sf_bessel_In`(int *n*, double *x*)

- int `gsl_sf_bessel_In_e`(int *n*, double *x*, [gsl_sf_result](https://www.gnu.org/software/gsl/doc/html/specfunc.html#c.gsl_sf_result) * *result*)

  These routines compute the regular modified cylindrical Bessel function of order `n`, ![I_n(x)](https://www.gnu.org/software/gsl/doc/html/_images/math/7e78a3c2669b08dafec3eee5c5339e16063277b3.png).

- int `gsl_sf_bessel_In_array`(int *nmin*, int *nmax*, double *x*, double *result_array[]*)

  This routine computes the values of the regular modified cylindrical Bessel functions ![I_n(x)](https://www.gnu.org/software/gsl/doc/html/_images/math/7e78a3c2669b08dafec3eee5c5339e16063277b3.png) for ![n](https://www.gnu.org/software/gsl/doc/html/_images/math/a24554f1502cf204e7aa24a6c064962c3504de48.png) from `nmin` to `nmax` inclusive, storing the results in the array `result_array`. The start of the range `nmin` must be positive or zero. The values are computed using recurrence relations for efficiency, and therefore may differ slightly from the exact values.

- double `gsl_sf_bessel_I0_scaled`(double *x*)

- int `gsl_sf_bessel_I0_scaled_e`(double *x*, [gsl_sf_result](https://www.gnu.org/software/gsl/doc/html/specfunc.html#c.gsl_sf_result) * *result*)

  These routines compute the scaled regular modified cylindrical Bessel function of zeroth order ![\exp(-|x|) I_0(x)](https://www.gnu.org/software/gsl/doc/html/_images/math/8f7b507a6490a4ce5e08a8c9408e227f3c65a279.png).

- double `gsl_sf_bessel_I1_scaled`(double *x*)

- int `gsl_sf_bessel_I1_scaled_e`(double *x*, [gsl_sf_result](https://www.gnu.org/software/gsl/doc/html/specfunc.html#c.gsl_sf_result) * *result*)

  These routines compute the scaled regular modified cylindrical Bessel function of first order ![\exp(-|x|) I_1(x)](https://www.gnu.org/software/gsl/doc/html/_images/math/c66613bb121b21a819cef201621b6c9eb72f8308.png).

- double `gsl_sf_bessel_In_scaled`(int *n*, double *x*)

- int `gsl_sf_bessel_In_scaled_e`(int *n*, double *x*, [gsl_sf_result](https://www.gnu.org/software/gsl/doc/html/specfunc.html#c.gsl_sf_result) * *result*)

  These routines compute the scaled regular modified cylindrical Bessel function of order `n`, ![\exp(-|x|) I_n(x)](https://www.gnu.org/software/gsl/doc/html/_images/math/ddbf7f4c3316d6fed6d467ef03b6ee10a62a909c.png)

- int `gsl_sf_bessel_In_scaled_array`(int *nmin*, int *nmax*, double *x*, double *result_array[]*)

  This routine computes the values of the scaled regular cylindrical Bessel functions ![\exp(-|x|) I_n(x)](https://www.gnu.org/software/gsl/doc/html/_images/math/ddbf7f4c3316d6fed6d467ef03b6ee10a62a909c.png) for ![n](https://www.gnu.org/software/gsl/doc/html/_images/math/a24554f1502cf204e7aa24a6c064962c3504de48.png) from `nmin` to `nmax` inclusive, storing the results in the array`result_array`. The start of the range `nmin` must be positive or zero. The values are computed using recurrence relations for efficiency, and therefore may differ slightly from the exact values.

### Irregular Modified Cylindrical Bessel Functions





- double `gsl_sf_bessel_K0`(double *x*)

- int `gsl_sf_bessel_K0_e`(double *x*, [gsl_sf_result](https://www.gnu.org/software/gsl/doc/html/specfunc.html#c.gsl_sf_result) * *result*)

  These routines compute the irregular modified cylindrical Bessel function of zeroth order, ![K_0(x)](https://www.gnu.org/software/gsl/doc/html/_images/math/3a05961e29eabc493ba2a1f13fa6822f31314eb3.png), for ![x > 0](https://www.gnu.org/software/gsl/doc/html/_images/math/fe20eedd0aef40a657ad04835ff74fa8a73efb56.png).

- double `gsl_sf_bessel_K1`(double *x*)

- int `gsl_sf_bessel_K1_e`(double *x*, [gsl_sf_result](https://www.gnu.org/software/gsl/doc/html/specfunc.html#c.gsl_sf_result) * *result*)

  These routines compute the irregular modified cylindrical Bessel function of first order, ![K_1(x)](https://www.gnu.org/software/gsl/doc/html/_images/math/0effdda8be31ffbc39e330f72ce18536df2e5390.png), for ![x > 0](https://www.gnu.org/software/gsl/doc/html/_images/math/fe20eedd0aef40a657ad04835ff74fa8a73efb56.png).

- double `gsl_sf_bessel_Kn`(int *n*, double *x*)

- int `gsl_sf_bessel_Kn_e`(int *n*, double *x*, [gsl_sf_result](https://www.gnu.org/software/gsl/doc/html/specfunc.html#c.gsl_sf_result) * *result*)

  These routines compute the irregular modified cylindrical Bessel function of order `n`, ![K_n(x)](https://www.gnu.org/software/gsl/doc/html/_images/math/f44f0a148e38562c9d66c1d7ef29a270cb8a3e75.png), for ![x > 0](https://www.gnu.org/software/gsl/doc/html/_images/math/fe20eedd0aef40a657ad04835ff74fa8a73efb56.png).

- int `gsl_sf_bessel_Kn_array`(int *nmin*, int *nmax*, double *x*, double *result_array[]*)

  This routine computes the values of the irregular modified cylindrical Bessel functions ![K_n(x)](https://www.gnu.org/software/gsl/doc/html/_images/math/f44f0a148e38562c9d66c1d7ef29a270cb8a3e75.png)for ![n](https://www.gnu.org/software/gsl/doc/html/_images/math/a24554f1502cf204e7aa24a6c064962c3504de48.png) from `nmin` to `nmax` inclusive, storing the results in the array `result_array`. The start of the range `nmin` must be positive or zero. The domain of the function is ![x>0](https://www.gnu.org/software/gsl/doc/html/_images/math/937b5840b2fb62e72aff4740c9d3c14db3b06720.png). The values are computed using recurrence relations for efficiency, and therefore may differ slightly from the exact values.

- double `gsl_sf_bessel_K0_scaled`(double *x*)

- int `gsl_sf_bessel_K0_scaled_e`(double *x*, [gsl_sf_result](https://www.gnu.org/software/gsl/doc/html/specfunc.html#c.gsl_sf_result) * *result*)

  These routines compute the scaled irregular modified cylindrical Bessel function of zeroth order ![\exp(x) K_0(x)](https://www.gnu.org/software/gsl/doc/html/_images/math/b18990b88369bc9d4e1574b41859bf650e584d89.png) for ![x>0](https://www.gnu.org/software/gsl/doc/html/_images/math/937b5840b2fb62e72aff4740c9d3c14db3b06720.png).

- double `gsl_sf_bessel_K1_scaled`(double *x*)

- int `gsl_sf_bessel_K1_scaled_e`(double *x*, [gsl_sf_result](https://www.gnu.org/software/gsl/doc/html/specfunc.html#c.gsl_sf_result) * *result*)

  These routines compute the scaled irregular modified cylindrical Bessel function of first order ![\exp(x) K_1(x)](https://www.gnu.org/software/gsl/doc/html/_images/math/1f79b422556852bd5fe6722e2e48e99226ecfc63.png) for ![x>0](https://www.gnu.org/software/gsl/doc/html/_images/math/937b5840b2fb62e72aff4740c9d3c14db3b06720.png).

- double `gsl_sf_bessel_Kn_scaled`(int *n*, double *x*)

- int `gsl_sf_bessel_Kn_scaled_e`(int *n*, double *x*, [gsl_sf_result](https://www.gnu.org/software/gsl/doc/html/specfunc.html#c.gsl_sf_result) * *result*)

  These routines compute the scaled irregular modified cylindrical Bessel function of order `n`, ![\exp(x) K_n(x)](https://www.gnu.org/software/gsl/doc/html/_images/math/7af61114297df5e6da7339c5591aef5a647e3aab.png), for ![x>0](https://www.gnu.org/software/gsl/doc/html/_images/math/937b5840b2fb62e72aff4740c9d3c14db3b06720.png).

- int `gsl_sf_bessel_Kn_scaled_array`(int *nmin*, int *nmax*, double *x*, double *result_array[]*)

  This routine computes the values of the scaled irregular cylindrical Bessel functions ![\exp(x) K_n(x)](https://www.gnu.org/software/gsl/doc/html/_images/math/7af61114297df5e6da7339c5591aef5a647e3aab.png) for ![n](https://www.gnu.org/software/gsl/doc/html/_images/math/a24554f1502cf204e7aa24a6c064962c3504de48.png) from `nmin` to `nmax` inclusive, storing the results in the array `result_array`. The start of the range `nmin` must be positive or zero. The domain of the function is ![x>0](https://www.gnu.org/software/gsl/doc/html/_images/math/937b5840b2fb62e72aff4740c9d3c14db3b06720.png). The values are computed using recurrence relations for efficiency, and therefore may differ slightly from the exact values.

### Regular Spherical Bessel Functions







- double `gsl_sf_bessel_j0`(double *x*)

- int `gsl_sf_bessel_j0_e`(double *x*, [gsl_sf_result](https://www.gnu.org/software/gsl/doc/html/specfunc.html#c.gsl_sf_result) * *result*)

  These routines compute the regular spherical Bessel function of zeroth order, ![j_0(x) = \sin(x)/x](https://www.gnu.org/software/gsl/doc/html/_images/math/ca0a868e0489377a954208a35de310a397be44a1.png).

- double `gsl_sf_bessel_j1`(double *x*)

- int `gsl_sf_bessel_j1_e`(double *x*, [gsl_sf_result](https://www.gnu.org/software/gsl/doc/html/specfunc.html#c.gsl_sf_result) * *result*)

  These routines compute the regular spherical Bessel function of first order, ![j_1(x) = (\sin(x)/x - \cos(x))/x](https://www.gnu.org/software/gsl/doc/html/_images/math/7201882b61f67117d1b60f07971ff71a3bfdd2fb.png).

- double `gsl_sf_bessel_j2`(double *x*)

- int `gsl_sf_bessel_j2_e`(double *x*, [gsl_sf_result](https://www.gnu.org/software/gsl/doc/html/specfunc.html#c.gsl_sf_result) * *result*)

  These routines compute the regular spherical Bessel function of second order, ![j_2(x) = ((3/x^2 - 1)\sin(x) - 3\cos(x)/x)/x](https://www.gnu.org/software/gsl/doc/html/_images/math/1f1a973b7e1b1c1119fbd58cc0857fa95d9aab69.png).

- double `gsl_sf_bessel_jl`(int *l*, double *x*)

- int `gsl_sf_bessel_jl_e`(int *l*, double *x*, [gsl_sf_result](https://www.gnu.org/software/gsl/doc/html/specfunc.html#c.gsl_sf_result) * *result*)

  These routines compute the regular spherical Bessel function of order `l`, ![j_l(x)](https://www.gnu.org/software/gsl/doc/html/_images/math/2c58211abd5178158a7ccf83ba0d0ded869cb455.png), for ![l \geq 0](https://www.gnu.org/software/gsl/doc/html/_images/math/666da554674fafe8302f33271a318ec0e370cc44.png) and ![x \geq 0](https://www.gnu.org/software/gsl/doc/html/_images/math/c40948eefa2d293144a5313ae9075b8ced7e1798.png).

- int `gsl_sf_bessel_jl_array`(int *lmax*, double *x*, double *result_array[]*)

  This routine computes the values of the regular spherical Bessel functions ![j_l(x)](https://www.gnu.org/software/gsl/doc/html/_images/math/2c58211abd5178158a7ccf83ba0d0ded869cb455.png) for ![l](https://www.gnu.org/software/gsl/doc/html/_images/math/a7cc908d8b1de27d7caaa2eab1f302b760d3dbf0.png) from 0 to `lmax` inclusive for ![lmax \geq 0](https://www.gnu.org/software/gsl/doc/html/_images/math/78b5d4cf68bd74513d084afbf083c682c343551c.png) and ![x \geq 0](https://www.gnu.org/software/gsl/doc/html/_images/math/c40948eefa2d293144a5313ae9075b8ced7e1798.png), storing the results in the array `result_array`. The values are computed using recurrence relations for efficiency, and therefore may differ slightly from the exact values.

- int `gsl_sf_bessel_jl_steed_array`(int *lmax*, double *x*, double * *result_array*)

  This routine uses Steed’s method to compute the values of the regular spherical Bessel functions ![j_l(x)](https://www.gnu.org/software/gsl/doc/html/_images/math/2c58211abd5178158a7ccf83ba0d0ded869cb455.png) for ![l](https://www.gnu.org/software/gsl/doc/html/_images/math/a7cc908d8b1de27d7caaa2eab1f302b760d3dbf0.png) from 0 to `lmax` inclusive for ![lmax \geq 0](https://www.gnu.org/software/gsl/doc/html/_images/math/78b5d4cf68bd74513d084afbf083c682c343551c.png) and ![x \geq 0](https://www.gnu.org/software/gsl/doc/html/_images/math/c40948eefa2d293144a5313ae9075b8ced7e1798.png), storing the results in the array `result_array`. The Steed/Barnett algorithm is described in Comp. Phys. Comm. 21, 297 (1981). Steed’s method is more stable than the recurrence used in the other functions but is also slower.

### Irregular Spherical Bessel Functions





- double `gsl_sf_bessel_y0`(double *x*)

- int `gsl_sf_bessel_y0_e`(double *x*, [gsl_sf_result](https://www.gnu.org/software/gsl/doc/html/specfunc.html#c.gsl_sf_result) * *result*)

  These routines compute the irregular spherical Bessel function of zeroth order, ![y_0(x) = -\cos(x)/x](https://www.gnu.org/software/gsl/doc/html/_images/math/846001d6ff5db08a7c5f949eb83bfedac9ba731c.png).

- double `gsl_sf_bessel_y1`(double *x*)

- int `gsl_sf_bessel_y1_e`(double *x*, [gsl_sf_result](https://www.gnu.org/software/gsl/doc/html/specfunc.html#c.gsl_sf_result) * *result*)

  These routines compute the irregular spherical Bessel function of first order, ![y_1(x) = -(\cos(x)/x + \sin(x))/x](https://www.gnu.org/software/gsl/doc/html/_images/math/52ae64ac3a6a5bdb676321b6aa72300225296b4d.png).

- double `gsl_sf_bessel_y2`(double *x*)

- int `gsl_sf_bessel_y2_e`(double *x*, [gsl_sf_result](https://www.gnu.org/software/gsl/doc/html/specfunc.html#c.gsl_sf_result) * *result*)

  These routines compute the irregular spherical Bessel function of second order, ![y_2(x) = (-3/x^3 + 1/x)\cos(x) - (3/x^2)\sin(x)](https://www.gnu.org/software/gsl/doc/html/_images/math/a629140a63b9b566cdba2f5be55ebbcdce72a383.png).

- double `gsl_sf_bessel_yl`(int *l*, double *x*)

- int `gsl_sf_bessel_yl_e`(int *l*, double *x*, [gsl_sf_result](https://www.gnu.org/software/gsl/doc/html/specfunc.html#c.gsl_sf_result) * *result*)

  These routines compute the irregular spherical Bessel function of order `l`, ![y_l(x)](https://www.gnu.org/software/gsl/doc/html/_images/math/a5237e6a66274662bef5fadebb08539203de721d.png), for ![l \geq 0](https://www.gnu.org/software/gsl/doc/html/_images/math/666da554674fafe8302f33271a318ec0e370cc44.png).

- int `gsl_sf_bessel_yl_array`(int *lmax*, double *x*, double *result_array[]*)

  This routine computes the values of the irregular spherical Bessel functions ![y_l(x)](https://www.gnu.org/software/gsl/doc/html/_images/math/a5237e6a66274662bef5fadebb08539203de721d.png) for ![l](https://www.gnu.org/software/gsl/doc/html/_images/math/a7cc908d8b1de27d7caaa2eab1f302b760d3dbf0.png) from 0 to `lmax` inclusive for ![lmax \geq 0](https://www.gnu.org/software/gsl/doc/html/_images/math/78b5d4cf68bd74513d084afbf083c682c343551c.png), storing the results in the array `result_array`. The values are computed using recurrence relations for efficiency, and therefore may differ slightly from the exact values.

### Regular Modified Spherical Bessel Functions





The regular modified spherical Bessel functions ![i_l(x)](https://www.gnu.org/software/gsl/doc/html/_images/math/66569b71033cd7295c2de3f85e2b52d48cb9c035.png) are related to the modified Bessel functions of fractional order, ![i_l(x) = \sqrt{\pi/(2x)} I_{l+1/2}(x)](https://www.gnu.org/software/gsl/doc/html/_images/math/d33ed445eba81a7558f41f281e63ccf3a828c7a5.png)

- double `gsl_sf_bessel_i0_scaled`(double *x*)

- int `gsl_sf_bessel_i0_scaled_e`(double *x*, [gsl_sf_result](https://www.gnu.org/software/gsl/doc/html/specfunc.html#c.gsl_sf_result) * *result*)

  These routines compute the scaled regular modified spherical Bessel function of zeroth order, ![\exp(-|x|) i_0(x)](https://www.gnu.org/software/gsl/doc/html/_images/math/e9571a4347a75cdb2c0cf9b984fd806088b4f436.png).

- double `gsl_sf_bessel_i1_scaled`(double *x*)

- int `gsl_sf_bessel_i1_scaled_e`(double *x*, [gsl_sf_result](https://www.gnu.org/software/gsl/doc/html/specfunc.html#c.gsl_sf_result) * *result*)

  These routines compute the scaled regular modified spherical Bessel function of first order, ![\exp(-|x|) i_1(x)](https://www.gnu.org/software/gsl/doc/html/_images/math/2abfa2826dece283a5f67582a7d986853a75b4e9.png).

- double `gsl_sf_bessel_i2_scaled`(double *x*)

- int `gsl_sf_bessel_i2_scaled_e`(double *x*, [gsl_sf_result](https://www.gnu.org/software/gsl/doc/html/specfunc.html#c.gsl_sf_result) * *result*)

  These routines compute the scaled regular modified spherical Bessel function of second order, ![\exp(-|x|) i_2(x)](https://www.gnu.org/software/gsl/doc/html/_images/math/9bc65f80ad5e65004a0553cd25cbfcda10d6e28c.png)

- double `gsl_sf_bessel_il_scaled`(int *l*, double *x*)

- int `gsl_sf_bessel_il_scaled_e`(int *l*, double *x*, [gsl_sf_result](https://www.gnu.org/software/gsl/doc/html/specfunc.html#c.gsl_sf_result) * *result*)

  These routines compute the scaled regular modified spherical Bessel function of order `l`, ![\exp(-|x|) i_l(x)](https://www.gnu.org/software/gsl/doc/html/_images/math/45121f39bedbe0c3aed8414207448d10b0f4f9f7.png)

- int `gsl_sf_bessel_il_scaled_array`(int *lmax*, double *x*, double *result_array[]*)

  This routine computes the values of the scaled regular modified spherical Bessel functions ![\exp(-|x|) i_l(x)](https://www.gnu.org/software/gsl/doc/html/_images/math/45121f39bedbe0c3aed8414207448d10b0f4f9f7.png) for ![l](https://www.gnu.org/software/gsl/doc/html/_images/math/a7cc908d8b1de27d7caaa2eab1f302b760d3dbf0.png) from 0 to `lmax` inclusive for ![lmax \geq 0](https://www.gnu.org/software/gsl/doc/html/_images/math/78b5d4cf68bd74513d084afbf083c682c343551c.png), storing the results in the array `result_array`. The values are computed using recurrence relations for efficiency, and therefore may differ slightly from the exact values.

### Irregular Modified Spherical Bessel Functions



The irregular modified spherical Bessel functions ![k_l(x)](https://www.gnu.org/software/gsl/doc/html/_images/math/0cae625b70bafb5ce956c6ddf2af444cdbdad4cd.png) are related to the irregular modified Bessel functions of fractional order, ![k_l(x) = \sqrt{\pi/(2x)} K_{l+1/2}(x)](https://www.gnu.org/software/gsl/doc/html/_images/math/1e58ca3e97c9c954d69a653bcbad1ac774864000.png).

- double `gsl_sf_bessel_k0_scaled`(double *x*)

- int `gsl_sf_bessel_k0_scaled_e`(double *x*, [gsl_sf_result](https://www.gnu.org/software/gsl/doc/html/specfunc.html#c.gsl_sf_result) * *result*)

  These routines compute the scaled irregular modified spherical Bessel function of zeroth order, ![\exp(x) k_0(x)](https://www.gnu.org/software/gsl/doc/html/_images/math/4860fb894259fe8191bc4125f902e90784fd1893.png), for ![x>0](https://www.gnu.org/software/gsl/doc/html/_images/math/937b5840b2fb62e72aff4740c9d3c14db3b06720.png).

- double `gsl_sf_bessel_k1_scaled`(double *x*)

- int `gsl_sf_bessel_k1_scaled_e`(double *x*, [gsl_sf_result](https://www.gnu.org/software/gsl/doc/html/specfunc.html#c.gsl_sf_result) * *result*)

  These routines compute the scaled irregular modified spherical Bessel function of first order, ![\exp(x) k_1(x)](https://www.gnu.org/software/gsl/doc/html/_images/math/f5e0e7531d61064d39ae19c0d3e85301f65c7d41.png), for ![x>0](https://www.gnu.org/software/gsl/doc/html/_images/math/937b5840b2fb62e72aff4740c9d3c14db3b06720.png).

- double `gsl_sf_bessel_k2_scaled`(double *x*)

- int `gsl_sf_bessel_k2_scaled_e`(double *x*, [gsl_sf_result](https://www.gnu.org/software/gsl/doc/html/specfunc.html#c.gsl_sf_result) * *result*)

  These routines compute the scaled irregular modified spherical Bessel function of second order, ![\exp(x) k_2(x)](https://www.gnu.org/software/gsl/doc/html/_images/math/178a2a091cf856070eed95c1ce075872108ea25e.png), for ![x>0](https://www.gnu.org/software/gsl/doc/html/_images/math/937b5840b2fb62e72aff4740c9d3c14db3b06720.png).

- double `gsl_sf_bessel_kl_scaled`(int *l*, double *x*)

- int `gsl_sf_bessel_kl_scaled_e`(int *l*, double *x*, [gsl_sf_result](https://www.gnu.org/software/gsl/doc/html/specfunc.html#c.gsl_sf_result) * *result*)

  These routines compute the scaled irregular modified spherical Bessel function of order `l`, ![\exp(x) k_l(x)](https://www.gnu.org/software/gsl/doc/html/_images/math/0bead5300acd45a4770cca211b247ea4844fb684.png), for ![x>0](https://www.gnu.org/software/gsl/doc/html/_images/math/937b5840b2fb62e72aff4740c9d3c14db3b06720.png).

- int `gsl_sf_bessel_kl_scaled_array`(int *lmax*, double *x*, double *result_array[]*)

  This routine computes the values of the scaled irregular modified spherical Bessel functions ![\exp(x) k_l(x)](https://www.gnu.org/software/gsl/doc/html/_images/math/0bead5300acd45a4770cca211b247ea4844fb684.png) for ![l](https://www.gnu.org/software/gsl/doc/html/_images/math/a7cc908d8b1de27d7caaa2eab1f302b760d3dbf0.png) from 0 to `lmax` inclusive for ![lmax \geq 0](https://www.gnu.org/software/gsl/doc/html/_images/math/78b5d4cf68bd74513d084afbf083c682c343551c.png) and ![x>0](https://www.gnu.org/software/gsl/doc/html/_images/math/937b5840b2fb62e72aff4740c9d3c14db3b06720.png), storing the results in the array `result_array`. The values are computed using recurrence relations for efficiency, and therefore may differ slightly from the exact values.

### Regular Bessel Function—Fractional Order



- double `gsl_sf_bessel_Jnu`(double *nu*, double *x*)

- int `gsl_sf_bessel_Jnu_e`(double *nu*, double *x*, [gsl_sf_result](https://www.gnu.org/software/gsl/doc/html/specfunc.html#c.gsl_sf_result) * *result*)

  These routines compute the regular cylindrical Bessel function of fractional order ![\nu](https://www.gnu.org/software/gsl/doc/html/_images/math/b8314a673ae9c2de05235962b598194cee073946.png), ![J_\nu(x)](https://www.gnu.org/software/gsl/doc/html/_images/math/6e63fa9eb40ca493d0af71ffc8dcf179a841bfb1.png).

- int `gsl_sf_bessel_sequence_Jnu_e`(double *nu*, [gsl_mode_t](https://www.gnu.org/software/gsl/doc/html/specfunc.html#c.gsl_mode_t) *mode*, size_t *size*, double *v[]*)

  This function computes the regular cylindrical Bessel function of fractional order ![\nu](https://www.gnu.org/software/gsl/doc/html/_images/math/b8314a673ae9c2de05235962b598194cee073946.png), ![J_\nu(x)](https://www.gnu.org/software/gsl/doc/html/_images/math/6e63fa9eb40ca493d0af71ffc8dcf179a841bfb1.png), evaluated at a series of ![x](https://www.gnu.org/software/gsl/doc/html/_images/math/1c5715c8edd0bbc330baa0f0494f1f7951619754.png) values. The array `v` of length `size` contains the ![x](https://www.gnu.org/software/gsl/doc/html/_images/math/1c5715c8edd0bbc330baa0f0494f1f7951619754.png) values. They are assumed to be strictly ordered and positive. The array is over-written with the values of ![J_\nu(x_i)](https://www.gnu.org/software/gsl/doc/html/_images/math/2a0bf878a4b9a86a010a4e37586d1a369deaec57.png).

### Irregular Bessel Functions—Fractional Order

- double `gsl_sf_bessel_Ynu`(double *nu*, double *x*)

- int `gsl_sf_bessel_Ynu_e`(double *nu*, double *x*, [gsl_sf_result](https://www.gnu.org/software/gsl/doc/html/specfunc.html#c.gsl_sf_result) * *result*)

  These routines compute the irregular cylindrical Bessel function of fractional order ![\nu](https://www.gnu.org/software/gsl/doc/html/_images/math/b8314a673ae9c2de05235962b598194cee073946.png), ![Y_\nu(x)](https://www.gnu.org/software/gsl/doc/html/_images/math/730ee736db30f7c565689d357151f5c6ba59fa08.png).

### Regular Modified Bessel Functions—Fractional Order



- double `gsl_sf_bessel_Inu`(double *nu*, double *x*)

- int `gsl_sf_bessel_Inu_e`(double *nu*, double *x*, [gsl_sf_result](https://www.gnu.org/software/gsl/doc/html/specfunc.html#c.gsl_sf_result) * *result*)

  These routines compute the regular modified Bessel function of fractional order ![\nu](https://www.gnu.org/software/gsl/doc/html/_images/math/b8314a673ae9c2de05235962b598194cee073946.png), ![I_\nu(x)](https://www.gnu.org/software/gsl/doc/html/_images/math/cf219186493c9cfed1b5e761568ccb8e2cbf7b47.png) for ![x>0](https://www.gnu.org/software/gsl/doc/html/_images/math/937b5840b2fb62e72aff4740c9d3c14db3b06720.png), ![\nu>0](https://www.gnu.org/software/gsl/doc/html/_images/math/e3b5e1a0b5b2b6c20de392328ac0deb6fa07ee53.png).

- double `gsl_sf_bessel_Inu_scaled`(double *nu*, double *x*)

- int `gsl_sf_bessel_Inu_scaled_e`(double *nu*, double *x*, [gsl_sf_result](https://www.gnu.org/software/gsl/doc/html/specfunc.html#c.gsl_sf_result) * *result*)

  These routines compute the scaled regular modified Bessel function of fractional order ![\nu](https://www.gnu.org/software/gsl/doc/html/_images/math/b8314a673ae9c2de05235962b598194cee073946.png), ![\exp(-|x|)I_\nu(x)](https://www.gnu.org/software/gsl/doc/html/_images/math/d02e508fb115f16b9f5269d0a8cc1febed3ca457.png) for ![x>0](https://www.gnu.org/software/gsl/doc/html/_images/math/937b5840b2fb62e72aff4740c9d3c14db3b06720.png), ![\nu>0](https://www.gnu.org/software/gsl/doc/html/_images/math/e3b5e1a0b5b2b6c20de392328ac0deb6fa07ee53.png).

### Irregular Modified Bessel Functions—Fractional Order



- double `gsl_sf_bessel_Knu`(double *nu*, double *x*)

- int `gsl_sf_bessel_Knu_e`(double *nu*, double *x*, [gsl_sf_result](https://www.gnu.org/software/gsl/doc/html/specfunc.html#c.gsl_sf_result) * *result*)

  These routines compute the irregular modified Bessel function of fractional order ![\nu](https://www.gnu.org/software/gsl/doc/html/_images/math/b8314a673ae9c2de05235962b598194cee073946.png), ![K_\nu(x)](https://www.gnu.org/software/gsl/doc/html/_images/math/36f752d5190b0034a4393732afdf2365e8e25613.png) for ![x>0](https://www.gnu.org/software/gsl/doc/html/_images/math/937b5840b2fb62e72aff4740c9d3c14db3b06720.png), ![\nu>0](https://www.gnu.org/software/gsl/doc/html/_images/math/e3b5e1a0b5b2b6c20de392328ac0deb6fa07ee53.png).

- double `gsl_sf_bessel_lnKnu`(double *nu*, double *x*)

- int `gsl_sf_bessel_lnKnu_e`(double *nu*, double *x*, [gsl_sf_result](https://www.gnu.org/software/gsl/doc/html/specfunc.html#c.gsl_sf_result) * *result*)

  These routines compute the logarithm of the irregular modified Bessel function of fractional order ![\nu](https://www.gnu.org/software/gsl/doc/html/_images/math/b8314a673ae9c2de05235962b598194cee073946.png), ![\ln(K_\nu(x))](https://www.gnu.org/software/gsl/doc/html/_images/math/00ff48b1cb9fcc9f2183df2b46f62456a026ef26.png) for ![x>0](https://www.gnu.org/software/gsl/doc/html/_images/math/937b5840b2fb62e72aff4740c9d3c14db3b06720.png), ![\nu>0](https://www.gnu.org/software/gsl/doc/html/_images/math/e3b5e1a0b5b2b6c20de392328ac0deb6fa07ee53.png).

- double `gsl_sf_bessel_Knu_scaled`(double *nu*, double *x*)

- int `gsl_sf_bessel_Knu_scaled_e`(double *nu*, double *x*, [gsl_sf_result](https://www.gnu.org/software/gsl/doc/html/specfunc.html#c.gsl_sf_result) * *result*)

  These routines compute the scaled irregular modified Bessel function of fractional order ![\nu](https://www.gnu.org/software/gsl/doc/html/_images/math/b8314a673ae9c2de05235962b598194cee073946.png), ![\exp(+|x|) K_\nu(x)](https://www.gnu.org/software/gsl/doc/html/_images/math/cdccb6c6ad646f7bbfd8e5cee5d8297c0da77218.png) for ![x>0](https://www.gnu.org/software/gsl/doc/html/_images/math/937b5840b2fb62e72aff4740c9d3c14db3b06720.png), ![\nu>0](https://www.gnu.org/software/gsl/doc/html/_images/math/e3b5e1a0b5b2b6c20de392328ac0deb6fa07ee53.png).

### Zeros of Regular Bessel Functions



- double `gsl_sf_bessel_zero_J0`(unsigned int *s*)

- int `gsl_sf_bessel_zero_J0_e`(unsigned int *s*, [gsl_sf_result](https://www.gnu.org/software/gsl/doc/html/specfunc.html#c.gsl_sf_result) * *result*)

  These routines compute the location of the `s`-th positive zero of the Bessel function ![J_0(x)](https://www.gnu.org/software/gsl/doc/html/_images/math/436b1799db46b1ecb38cc2fde801e882814c7119.png).

- double `gsl_sf_bessel_zero_J1`(unsigned int *s*)

- int `gsl_sf_bessel_zero_J1_e`(unsigned int *s*, [gsl_sf_result](https://www.gnu.org/software/gsl/doc/html/specfunc.html#c.gsl_sf_result) * *result*)

  These routines compute the location of the `s`-th positive zero of the Bessel function ![J_1(x)](https://www.gnu.org/software/gsl/doc/html/_images/math/cba66ff642120cf89928250899e1bf9a8dd258c2.png).

- double `gsl_sf_bessel_zero_Jnu`(double *nu*, unsigned int *s*)

- int `gsl_sf_bessel_zero_Jnu_e`(double *nu*, unsigned int *s*, [gsl_sf_result](https://www.gnu.org/software/gsl/doc/html/specfunc.html#c.gsl_sf_result) * *result*)

  These routines compute the location of the `s`-th positive zero of the Bessel function ![J_\nu(x)](https://www.gnu.org/software/gsl/doc/html/_images/math/6e63fa9eb40ca493d0af71ffc8dcf179a841bfb1.png). The current implementation does not support negative values of `nu`.

## Clausen Functions

The Clausen function is defined by the following integral,

![Cl_2(x) = - \int_0^x dt \log{\left( 2 \sin{(t/2)} \right)}](https://www.gnu.org/software/gsl/doc/html/_images/math/ea081af5c019911c57b057319c5ae3878d62bd37.png)

It is related to the [dilogarithm](https://www.gnu.org/software/gsl/doc/html/specfunc.html#dilog-function) by ![Cl_2(\theta) = \Im Li_2(\exp(i\theta))](https://www.gnu.org/software/gsl/doc/html/_images/math/96bb623ec14a8fae798f6d9b93bff9fd8f43650b.png). The Clausen functions are declared in the header file `gsl_sf_clausen.h`.

- double `gsl_sf_clausen`(double *x*)

- int `gsl_sf_clausen_e`(double *x*, [gsl_sf_result](https://www.gnu.org/software/gsl/doc/html/specfunc.html#c.gsl_sf_result) * *result*)

  These routines compute the Clausen integral ![Cl_2(x)](https://www.gnu.org/software/gsl/doc/html/_images/math/20e1d4fe83fca8d3f250051432c1d2d6bbd58735.png).

## Coulomb Functions

The prototypes of the Coulomb functions are declared in the header file `gsl_sf_coulomb.h`. Both bound state and scattering solutions are available.

### Normalized Hydrogenic Bound States

- double `gsl_sf_hydrogenicR_1`(double *Z*, double *r*)

- int `gsl_sf_hydrogenicR_1_e`(double *Z*, double *r*, [gsl_sf_result](https://www.gnu.org/software/gsl/doc/html/specfunc.html#c.gsl_sf_result) * *result*)

  These routines compute the lowest-order normalized hydrogenic bound state radial wavefunction ![R_1 := 2Z \sqrt{Z} \exp(-Z r)](https://www.gnu.org/software/gsl/doc/html/_images/math/d8079535e5693847840e2bf2cae619b21b3b72ea.png).

- double `gsl_sf_hydrogenicR`(int *n*, int *l*, double *Z*, double *r*)

- int `gsl_sf_hydrogenicR_e`(int *n*, int *l*, double *Z*, double *r*, [gsl_sf_result](https://www.gnu.org/software/gsl/doc/html/specfunc.html#c.gsl_sf_result) * *result*)

  These routines compute the `n`-th normalized hydrogenic bound state radial wavefunction,![R_n := {2 Z^{3/2} \over n^2}  \left({2Z r \over n}\right)^l  \sqrt{(n-l-1)! \over (n+l)!} \exp(-Z r/n) L^{2l+1}_{n-l-1}(2Z r / n).](https://www.gnu.org/software/gsl/doc/html/_images/math/4d5f4ec277ea707a54dd8997d0cc51bfb6b423ac.png)where ![L^a_b(x)](https://www.gnu.org/software/gsl/doc/html/_images/math/1086ce668c54b26be30ee0ab5e5728d60bd52abf.png) is the [generalized Laguerre polynomial](https://www.gnu.org/software/gsl/doc/html/specfunc.html#laguerre-functions). The normalization is chosen such that the wavefunction ![\psi](https://www.gnu.org/software/gsl/doc/html/_images/math/ad4bb05f6787ca6228964e118165614ed5b10cb9.png) is given by ![\psi(n,l,r) = R_n Y_{lm}](https://www.gnu.org/software/gsl/doc/html/_images/math/6d0c94b9e67b235d88f9dde1f751a669bbafb7d2.png).

### Coulomb Wave Functions

The Coulomb wave functions ![F_L(\eta,x)](https://www.gnu.org/software/gsl/doc/html/_images/math/adeac6e6774635defab1eb89bf45e6d07bb65234.png), ![G_L(\eta,x)](https://www.gnu.org/software/gsl/doc/html/_images/math/d559c79d55cbb39dc34e82e01f77afe9b3a32f25.png) are described in Abramowitz & Stegun, Chapter 14. Because there can be a large dynamic range of values for these functions, overflows are handled gracefully. If an overflow occurs, `GSL_EOVRFLW` is signalled and exponent(s) are returned through the modifiable parameters `exp_F`, `exp_G`. The full solution can be reconstructed from the following relations,

![F_L(\eta,x) &= fc[k_L] * \exp(exp_F) \\ G_L(\eta,x) &= gc[k_L] * \exp(exp_G)](https://www.gnu.org/software/gsl/doc/html/_images/math/be53482f2bda2d071a10362f328d6aa8b0df026b.png)

![F_L'(\eta,x) &= fcp[k_L] * \exp(exp_F) \\ G_L'(\eta,x) &= gcp[k_L] * \exp(exp_G)](https://www.gnu.org/software/gsl/doc/html/_images/math/495197f0418465174a852fbffe4a1d605f4ebb45.png)

- int `gsl_sf_coulomb_wave_FG_e`(double *eta*, double *x*, double *L_F*, int *k*, [gsl_sf_result](https://www.gnu.org/software/gsl/doc/html/specfunc.html#c.gsl_sf_result) * *F*, [gsl_sf_result](https://www.gnu.org/software/gsl/doc/html/specfunc.html#c.gsl_sf_result)* *Fp*, [gsl_sf_result](https://www.gnu.org/software/gsl/doc/html/specfunc.html#c.gsl_sf_result) * *G*, [gsl_sf_result](https://www.gnu.org/software/gsl/doc/html/specfunc.html#c.gsl_sf_result) * *Gp*, double * *exp_F*, double * *exp_G*)

  This function computes the Coulomb wave functions ![F_L(\eta,x)](https://www.gnu.org/software/gsl/doc/html/_images/math/adeac6e6774635defab1eb89bf45e6d07bb65234.png), ![G_{L-k}(\eta,x)](https://www.gnu.org/software/gsl/doc/html/_images/math/44fe3793d138a48e6169882249752a9ceb5baa76.png) and their derivatives![F'_L(\eta,x)](https://www.gnu.org/software/gsl/doc/html/_images/math/54d1c2233c877cb3237716fbde1a2bbeb3d4516c.png), ![G'_{L-k}(\eta,x)](https://www.gnu.org/software/gsl/doc/html/_images/math/044ff984c1b94c26556dea11383bf5aa90a9b943.png) with respect to ![x](https://www.gnu.org/software/gsl/doc/html/_images/math/1c5715c8edd0bbc330baa0f0494f1f7951619754.png). The parameters are restricted to ![L, L-k > -1/2](https://www.gnu.org/software/gsl/doc/html/_images/math/388160d7e263d4228ebd46224cce898e3b03dab8.png), ![x > 0](https://www.gnu.org/software/gsl/doc/html/_images/math/fe20eedd0aef40a657ad04835ff74fa8a73efb56.png) and integer ![k](https://www.gnu.org/software/gsl/doc/html/_images/math/90645b5c2abc13d7de104fd87abeb2f19406e5f8.png). Note that ![L](https://www.gnu.org/software/gsl/doc/html/_images/math/82359a7c56cc02d70a5ed0e98df631577f7d172b.png) itself is not restricted to being an integer. The results are stored in the parameters F, G for the function values and `Fp`, `Gp` for the derivative values. If an overflow occurs, `GSL_EOVRFLW` is returned and scaling exponents are stored in the modifiable parameters `exp_F`, `exp_G`.

- int `gsl_sf_coulomb_wave_F_array`(double *L_min*, int *kmax*, double *eta*, double *x*, double *fc_array[]*, double * *F_exponent*)

  This function computes the Coulomb wave function ![F_L(\eta,x)](https://www.gnu.org/software/gsl/doc/html/_images/math/adeac6e6774635defab1eb89bf45e6d07bb65234.png) for ![L = Lmin \dots Lmin + kmax](https://www.gnu.org/software/gsl/doc/html/_images/math/d28fc5b0272531735e647d95c92772f8c6099d21.png), storing the results in `fc_array`. In the case of overflow the exponent is stored in `F_exponent`.

- int `gsl_sf_coulomb_wave_FG_array`(double *L_min*, int *kmax*, double *eta*, double *x*, double *fc_array[]*, double *gc_array[]*, double * *F_exponent*, double * *G_exponent*)

  This function computes the functions ![F_L(\eta,x)](https://www.gnu.org/software/gsl/doc/html/_images/math/adeac6e6774635defab1eb89bf45e6d07bb65234.png), ![G_L(\eta,x)](https://www.gnu.org/software/gsl/doc/html/_images/math/d559c79d55cbb39dc34e82e01f77afe9b3a32f25.png) for ![L = Lmin \dots Lmin + kmax](https://www.gnu.org/software/gsl/doc/html/_images/math/d28fc5b0272531735e647d95c92772f8c6099d21.png)storing the results in `fc_array` and `gc_array`. In the case of overflow the exponents are stored in `F_exponent` and `G_exponent`.

- int `gsl_sf_coulomb_wave_FGp_array`(double *L_min*, int *kmax*, double *eta*, double *x*, double *fc_array[]*, double *fcp_array[]*, double *gc_array[]*, double *gcp_array[]*, double * *F_exponent*, double * *G_exponent*)

  This function computes the functions ![F_L(\eta,x)](https://www.gnu.org/software/gsl/doc/html/_images/math/adeac6e6774635defab1eb89bf45e6d07bb65234.png), ![G_L(\eta,x)](https://www.gnu.org/software/gsl/doc/html/_images/math/d559c79d55cbb39dc34e82e01f77afe9b3a32f25.png) and their derivatives ![F'_L(\eta,x)](https://www.gnu.org/software/gsl/doc/html/_images/math/54d1c2233c877cb3237716fbde1a2bbeb3d4516c.png), ![G'_L(\eta,x)](https://www.gnu.org/software/gsl/doc/html/_images/math/0adf1322e387fe63ac588a6a472c573ed2b1a50d.png) for ![L = Lmin \dots Lmin + kmax](https://www.gnu.org/software/gsl/doc/html/_images/math/d28fc5b0272531735e647d95c92772f8c6099d21.png) storing the results in `fc_array`, `gc_array`, `fcp_array` and `gcp_array`. In the case of overflow the exponents are stored in `F_exponent` and `G_exponent`.

- int `gsl_sf_coulomb_wave_sphF_array`(double *L_min*, int *kmax*, double *eta*, double *x*, double *fc_array[]*, double *F_exponent[]*)

  This function computes the Coulomb wave function divided by the argument ![F_L(\eta, x)/x](https://www.gnu.org/software/gsl/doc/html/_images/math/b244e61a260cd48bc6853721bd72766daeea666b.png) for ![L = Lmin \dots Lmin + kmax](https://www.gnu.org/software/gsl/doc/html/_images/math/d28fc5b0272531735e647d95c92772f8c6099d21.png), storing the results in `fc_array`. In the case of overflow the exponent is stored in `F_exponent`. This function reduces to spherical Bessel functions in the limit ![\eta \to 0](https://www.gnu.org/software/gsl/doc/html/_images/math/b9f7f2f213e2788f0e51a5c8ecf7d7243c60c1d2.png).

### Coulomb Wave Function Normalization Constant

The Coulomb wave function normalization constant is defined in Abramowitz 14.1.7.

- int `gsl_sf_coulomb_CL_e`(double *L*, double *eta*, [gsl_sf_result](https://www.gnu.org/software/gsl/doc/html/specfunc.html#c.gsl_sf_result) * *result*)

  This function computes the Coulomb wave function normalization constant ![C_L(\eta)](https://www.gnu.org/software/gsl/doc/html/_images/math/bf26a926d6f1dcba66f35ef8c5164e8f001490cc.png) for ![L > -1](https://www.gnu.org/software/gsl/doc/html/_images/math/18d5ac22d92bcd6f076916f184e8fdebfc3ce837.png).

- int `gsl_sf_coulomb_CL_array`(double *Lmin*, int *kmax*, double *eta*, double *cl[]*)

  This function computes the Coulomb wave function normalization constant ![C_L(\eta)](https://www.gnu.org/software/gsl/doc/html/_images/math/bf26a926d6f1dcba66f35ef8c5164e8f001490cc.png) for ![L = Lmin \dots Lmin + kmax](https://www.gnu.org/software/gsl/doc/html/_images/math/d28fc5b0272531735e647d95c92772f8c6099d21.png), ![Lmin > -1](https://www.gnu.org/software/gsl/doc/html/_images/math/8a3b9b432920a2854871b772edc1caccca0d9cc1.png).

## Coupling Coefficients

The Wigner 3-j, 6-j and 9-j symbols give the coupling coefficients for combined angular momentum vectors. Since the arguments of the standard coupling coefficient functions are integer or half-integer, the arguments of the following functions are, by convention, integers equal to twice the actual spin value. For information on the 3-j coefficients see Abramowitz & Stegun, Section 27.9. The functions described in this section are declared in the header file `gsl_sf_coupling.h`.

### 3-j Symbols

- double `gsl_sf_coupling_3j`(int *two_ja*, int *two_jb*, int *two_jc*, int *two_ma*, int *two_mb*, int *two_mc*)

- int `gsl_sf_coupling_3j_e`(int *two_ja*, int *two_jb*, int *two_jc*, int *two_ma*, int *two_mb*, int *two_mc*, [gsl_sf_result](https://www.gnu.org/software/gsl/doc/html/specfunc.html#c.gsl_sf_result) * *result*)

  These routines compute the Wigner 3-j coefficient,![\left( \begin{array}{ccc}   ja & jb & jc \\   ma & mb & mc \end{array} \right)](https://www.gnu.org/software/gsl/doc/html/_images/math/18bc4dce0baee085ab6c6f0d0ab9fad0a5af26db.png)where the arguments are given in half-integer units, ![ja](https://www.gnu.org/software/gsl/doc/html/_images/math/22238bf5695828122fea788d666911edd1a57b78.png) = `two_ja`/2, ![ma](https://www.gnu.org/software/gsl/doc/html/_images/math/3c71c7e83c7e4476d3bcc472696ca1a640d69eb6.png) = `two_ma`/2, etc.

### 6-j Symbols

- double `gsl_sf_coupling_6j`(int *two_ja*, int *two_jb*, int *two_jc*, int *two_jd*, int *two_je*, int *two_jf*)

- int `gsl_sf_coupling_6j_e`(int *two_ja*, int *two_jb*, int *two_jc*, int *two_jd*, int *two_je*, int *two_jf*, [gsl_sf_result](https://www.gnu.org/software/gsl/doc/html/specfunc.html#c.gsl_sf_result) * *result*)

  These routines compute the Wigner 6-j coefficient,![\left\{ \begin{array}{ccc}   ja & jb & jc \\   jd & je & jf \end{array} \right\}](https://www.gnu.org/software/gsl/doc/html/_images/math/209d4823bf38bca06e321688192e90fc5659be30.png)where the arguments are given in half-integer units, ![ja](https://www.gnu.org/software/gsl/doc/html/_images/math/22238bf5695828122fea788d666911edd1a57b78.png) = `two_ja`/2, ![ma](https://www.gnu.org/software/gsl/doc/html/_images/math/3c71c7e83c7e4476d3bcc472696ca1a640d69eb6.png) = `two_ma`/2, etc.

### 9-j Symbols

- double `gsl_sf_coupling_9j`(int *two_ja*, int *two_jb*, int *two_jc*, int *two_jd*, int *two_je*, int *two_jf*, int *two_jg*, int *two_jh*, int *two_ji*)

- int `gsl_sf_coupling_9j_e`(int *two_ja*, int *two_jb*, int *two_jc*, int *two_jd*, int *two_je*, int *two_jf*, int *two_jg*, int *two_jh*, int *two_ji*, [gsl_sf_result](https://www.gnu.org/software/gsl/doc/html/specfunc.html#c.gsl_sf_result) * *result*)

  These routines compute the Wigner 9-j coefficient,![\left\{ \begin{array}{ccc}   ja & jb & jc \\   jd & je & jf \\   jg & jh & ji \end{array} \right\}](https://www.gnu.org/software/gsl/doc/html/_images/math/a1af4fac2d8f79e54ad4599239e2ec478b99db78.png)where the arguments are given in half-integer units, ![ja](https://www.gnu.org/software/gsl/doc/html/_images/math/22238bf5695828122fea788d666911edd1a57b78.png) = `two_ja`/2, ![ma](https://www.gnu.org/software/gsl/doc/html/_images/math/3c71c7e83c7e4476d3bcc472696ca1a640d69eb6.png) = `two_ma`/2, etc.

## Dawson Function

The Dawson integral is defined by

![\exp(-x^2) \int_0^x dt \exp(t^2)](https://www.gnu.org/software/gsl/doc/html/_images/math/47379b1d7c136e6d59f62efa04a3dc7d26aa98aa.png)

A table of Dawson’s integral can be found in Abramowitz & Stegun, Table 7.5. The Dawson functions are declared in the header file `gsl_sf_dawson.h`.

- double `gsl_sf_dawson`(double *x*)

- int `gsl_sf_dawson_e`(double *x*, [gsl_sf_result](https://www.gnu.org/software/gsl/doc/html/specfunc.html#c.gsl_sf_result) * *result*)

  These routines compute the value of Dawson’s integral for `x`.

## Debye Functions

The Debye functions ![D_n(x)](https://www.gnu.org/software/gsl/doc/html/_images/math/8d17c313f756651314139b9864b64d60bdd37c92.png) are defined by the following integral,

![D_n(x) = {n \over x^n} \int_0^x dt {t^n \over e^t - 1}](https://www.gnu.org/software/gsl/doc/html/_images/math/3f2c1690a8236258c09e5fbb2e4e9804a6b41524.png)

For further information see Abramowitz & Stegun, Section 27.1. The Debye functions are declared in the header file `gsl_sf_debye.h`.

- double `gsl_sf_debye_1`(double *x*)

- int `gsl_sf_debye_1_e`(double *x*, [gsl_sf_result](https://www.gnu.org/software/gsl/doc/html/specfunc.html#c.gsl_sf_result) * *result*)

  These routines compute the first-order Debye function ![D_1(x)](https://www.gnu.org/software/gsl/doc/html/_images/math/627b94ff9c89d2e8e0595a1d52d8e2059354b040.png).

- double `gsl_sf_debye_2`(double *x*)

- int `gsl_sf_debye_2_e`(double *x*, [gsl_sf_result](https://www.gnu.org/software/gsl/doc/html/specfunc.html#c.gsl_sf_result) * *result*)

  These routines compute the second-order Debye function ![D_2(x)](https://www.gnu.org/software/gsl/doc/html/_images/math/d45b3ec09efcf32484682e4377381fb3ccb131ef.png).

- double `gsl_sf_debye_3`(double *x*)

- int `gsl_sf_debye_3_e`(double *x*, [gsl_sf_result](https://www.gnu.org/software/gsl/doc/html/specfunc.html#c.gsl_sf_result) * *result*)

  These routines compute the third-order Debye function ![D_3(x)](https://www.gnu.org/software/gsl/doc/html/_images/math/f10629c69d05704caf730e8692709d25e5579e2c.png).

- double `gsl_sf_debye_4`(double *x*)

- int `gsl_sf_debye_4_e`(double *x*, [gsl_sf_result](https://www.gnu.org/software/gsl/doc/html/specfunc.html#c.gsl_sf_result) * *result*)

  These routines compute the fourth-order Debye function ![D_4(x)](https://www.gnu.org/software/gsl/doc/html/_images/math/a47512fee325cf5facf4e82d21f873d59c6088c0.png).

- double `gsl_sf_debye_5`(double *x*)

- int `gsl_sf_debye_5_e`(double *x*, [gsl_sf_result](https://www.gnu.org/software/gsl/doc/html/specfunc.html#c.gsl_sf_result) * *result*)

  These routines compute the fifth-order Debye function ![D_5(x)](https://www.gnu.org/software/gsl/doc/html/_images/math/2db8cf35d4445ae412682c477f94d54beeb6a602.png).

- double `gsl_sf_debye_6`(double *x*)

- int `gsl_sf_debye_6_e`(double *x*, [gsl_sf_result](https://www.gnu.org/software/gsl/doc/html/specfunc.html#c.gsl_sf_result) * *result*)

  These routines compute the sixth-order Debye function ![D_6(x)](https://www.gnu.org/software/gsl/doc/html/_images/math/e50319e165456d724d175768b5eb79cc59b9322f.png).



## Dilogarithm

The dilogarithm is defined as

![Li_2(z) = - \int_0^z ds {\log{(1-s)} \over s}](https://www.gnu.org/software/gsl/doc/html/_images/math/4f3aa8c780bd0ef55e8ac1fa0ea11e111bb9420c.png)

The functions described in this section are declared in the header file `gsl_sf_dilog.h`.

### Real Argument

- double `gsl_sf_dilog`(double *x*)

- int `gsl_sf_dilog_e`(double *x*, [gsl_sf_result](https://www.gnu.org/software/gsl/doc/html/specfunc.html#c.gsl_sf_result) * *result*)

  These routines compute the dilogarithm for a real argument. In Lewin’s notation this is ![Li_2(x)](https://www.gnu.org/software/gsl/doc/html/_images/math/7c92315ba8351a9386b396fc26b873ba7572b337.png), the real part of the dilogarithm of a real ![x](https://www.gnu.org/software/gsl/doc/html/_images/math/1c5715c8edd0bbc330baa0f0494f1f7951619754.png). It is defined by the integral representation![Li_2(x) = - \Re \int_0^x ds \log(1-s) / s](https://www.gnu.org/software/gsl/doc/html/_images/math/404ed5be516a369ecec542d9b067032e7027a0b1.png)Note that ![\Im(Li_2(x)) = 0](https://www.gnu.org/software/gsl/doc/html/_images/math/41dd2a9ac37083f26f7dea867ff41c0cd0dbfb58.png) for ![x \le 1](https://www.gnu.org/software/gsl/doc/html/_images/math/af51b2cc76aa1f1673ba23686d949d7e09d1bbd7.png), and ![-\pi\log(x)](https://www.gnu.org/software/gsl/doc/html/_images/math/7ba5684b0fb2e04dcc545e1de0f0e3cdcb3f548b.png) for ![x > 1](https://www.gnu.org/software/gsl/doc/html/_images/math/d6a25338ca0f2a6fb29e9110ca0b76d103f6463f.png).Note that Abramowitz & Stegun refer to the Spence integral ![S(x) = Li_2(1 - x)](https://www.gnu.org/software/gsl/doc/html/_images/math/aeef5c1759286da6105ed88574fa26252ab7f190.png) as the dilogarithm rather than ![Li_2(x)](https://www.gnu.org/software/gsl/doc/html/_images/math/7c92315ba8351a9386b396fc26b873ba7572b337.png).

### Complex Argument

- int `gsl_sf_complex_dilog_e`(double *r*, double *theta*, [gsl_sf_result](https://www.gnu.org/software/gsl/doc/html/specfunc.html#c.gsl_sf_result) * *result_re*, [gsl_sf_result](https://www.gnu.org/software/gsl/doc/html/specfunc.html#c.gsl_sf_result) * *result_im*)

  This function computes the full complex-valued dilogarithm for the complex argument ![z = r \exp(i \theta)](https://www.gnu.org/software/gsl/doc/html/_images/math/c7d264113099aa36d113b3ab51c972dff2517edc.png). The real and imaginary parts of the result are returned in `result_re`, `result_im`.

## Elementary Operations



The following functions allow for the propagation of errors when combining quantities by multiplication. The functions are declared in the header file `gsl_sf_elementary.h`.

- double `gsl_sf_multiply`(double *x*, double *y*)

- int `gsl_sf_multiply_e`(double *x*, double *y*, [gsl_sf_result](https://www.gnu.org/software/gsl/doc/html/specfunc.html#c.gsl_sf_result) * *result*)

  This function multiplies `x` and `y` storing the product and its associated error in `result`.

- int `gsl_sf_multiply_err_e`(double *x*, double *dx*, double *y*, double *dy*, [gsl_sf_result](https://www.gnu.org/software/gsl/doc/html/specfunc.html#c.gsl_sf_result) * *result*)

  This function multiplies `x` and `y` with associated absolute errors `dx` and `dy`. The product ![xy \pm xy \sqrt{(dx/x)^2 +(dy/y)^2}](https://www.gnu.org/software/gsl/doc/html/_images/math/7c90e361853d0636106c4426ad1c429fa8d285b0.png) is stored in `result`.

## Elliptic Integrals

The functions described in this section are declared in the header file `gsl_sf_ellint.h`. Further information about the elliptic integrals can be found in Abramowitz & Stegun, Chapter 17.

### Definition of Legendre Forms

The Legendre forms of elliptic integrals ![F(\phi,k)](https://www.gnu.org/software/gsl/doc/html/_images/math/5e447e3b2cba8230362e65c49302a2e312186a13.png), ![E(\phi,k)](https://www.gnu.org/software/gsl/doc/html/_images/math/ee8224eb0d55eb0dab4e8d270ca1f764b2b9a677.png) and ![\Pi(\phi,k,n)](https://www.gnu.org/software/gsl/doc/html/_images/math/4684364912954c9c8aac26813ae3f2ee297e093c.png) are defined by,

![F(\phi,k)   &= \int_0^\phi dt {1 \over \sqrt{(1 - k^2 \sin^2(t))}} \\ E(\phi,k)   &= \int_0^\phi dt   \sqrt{(1 - k^2 \sin^2(t))} \\ \Pi(\phi,k,n) &= \int_0^\phi dt {1 \over (1 + n \sin^2(t)) \sqrt{1 - k^2 \sin^2(t)}}](https://www.gnu.org/software/gsl/doc/html/_images/math/0ace4878e091aabca13a6c9a40b9965523dfdfc9.png)

The complete Legendre forms are denoted by ![K(k) = F(\pi/2, k)](https://www.gnu.org/software/gsl/doc/html/_images/math/f24d101353e6f6713f5952691fe92e140d0c9710.png) and ![E(k) = E(\pi/2, k)](https://www.gnu.org/software/gsl/doc/html/_images/math/4afa2bb3949f2684570c026cf54939659097c007.png).

The notation used here is based on Carlson, “Numerische Mathematik” 33 (1979) 1 and differs slightly from that used by Abramowitz & Stegun, where the functions are given in terms of the parameter ![m = k^2](https://www.gnu.org/software/gsl/doc/html/_images/math/9d5667426cec24427e500eb97a3096029ad992e9.png) and ![n](https://www.gnu.org/software/gsl/doc/html/_images/math/a24554f1502cf204e7aa24a6c064962c3504de48.png) is replaced by ![-n](https://www.gnu.org/software/gsl/doc/html/_images/math/5f9cb80424cf824b016af9b551419eae94ca974f.png).

### Definition of Carlson Forms

The Carlson symmetric forms of elliptical integrals ![RC(x,y)](https://www.gnu.org/software/gsl/doc/html/_images/math/1ba238bde832a3013a8245cd4968e3bdf8f31bfd.png), ![RD(x,y,z)](https://www.gnu.org/software/gsl/doc/html/_images/math/d01cff478d89a62023a0c37832307a9b13d28aab.png), ![RF(x,y,z)](https://www.gnu.org/software/gsl/doc/html/_images/math/73a62113b3068835188237682b0fc2a7a0961159.png) and ![RJ(x,y,z,p)](https://www.gnu.org/software/gsl/doc/html/_images/math/2f3ddc98d528fff74d7369ecd655c003cb54a8be.png) are defined by,

![RC(x,y)   &= 1/2 \int_0^\infty dt (t+x)^{-1/2} (t+y)^{-1} \\ RD(x,y,z) &= 3/2 \int_0^\infty dt (t+x)^{-1/2} (t+y)^{-1/2} (t+z)^{-3/2} \\ RF(x,y,z) &= 1/2 \int_0^\infty dt (t+x)^{-1/2} (t+y)^{-1/2} (t+z)^{-1/2} \\ RJ(x,y,z,p) &= 3/2 \int_0^\infty dt (t+x)^{-1/2} (t+y)^{-1/2} (t+z)^{-1/2} (t+p)^{-1}](https://www.gnu.org/software/gsl/doc/html/_images/math/01a22a5de9b3c268af02eeddd274db257d26795f.png)

### Legendre Form of Complete Elliptic Integrals

- double `gsl_sf_ellint_Kcomp`(double *k*, [gsl_mode_t](https://www.gnu.org/software/gsl/doc/html/specfunc.html#c.gsl_mode_t) *mode*)

- int `gsl_sf_ellint_Kcomp_e`(double *k*, [gsl_mode_t](https://www.gnu.org/software/gsl/doc/html/specfunc.html#c.gsl_mode_t) *mode*, [gsl_sf_result](https://www.gnu.org/software/gsl/doc/html/specfunc.html#c.gsl_sf_result) * *result*)

  These routines compute the complete elliptic integral ![K(k)](https://www.gnu.org/software/gsl/doc/html/_images/math/6a10339baae0e1e1a76fd45517d6eea3e4081a84.png) to the accuracy specified by the mode variable [`mode`](https://www.gnu.org/software/gsl/doc/html/montecarlo.html#c.mode). Note that Abramowitz & Stegun define this function in terms of the parameter ![m = k^2](https://www.gnu.org/software/gsl/doc/html/_images/math/9d5667426cec24427e500eb97a3096029ad992e9.png).

- double `gsl_sf_ellint_Ecomp`(double *k*, [gsl_mode_t](https://www.gnu.org/software/gsl/doc/html/specfunc.html#c.gsl_mode_t) *mode*)

- int `gsl_sf_ellint_Ecomp_e`(double *k*, [gsl_mode_t](https://www.gnu.org/software/gsl/doc/html/specfunc.html#c.gsl_mode_t) *mode*, [gsl_sf_result](https://www.gnu.org/software/gsl/doc/html/specfunc.html#c.gsl_sf_result) * *result*)

  These routines compute the complete elliptic integral ![E(k)](https://www.gnu.org/software/gsl/doc/html/_images/math/24ac77f7ca35cffadcd52e4c4594cc6d53531a7b.png) to the accuracy specified by the mode variable [`mode`](https://www.gnu.org/software/gsl/doc/html/montecarlo.html#c.mode). Note that Abramowitz & Stegun define this function in terms of the parameter ![m = k^2](https://www.gnu.org/software/gsl/doc/html/_images/math/9d5667426cec24427e500eb97a3096029ad992e9.png).

- double `gsl_sf_ellint_Pcomp`(double *k*, double *n*, [gsl_mode_t](https://www.gnu.org/software/gsl/doc/html/specfunc.html#c.gsl_mode_t) *mode*)

- int `gsl_sf_ellint_Pcomp_e`(double *k*, double *n*, [gsl_mode_t](https://www.gnu.org/software/gsl/doc/html/specfunc.html#c.gsl_mode_t) *mode*, [gsl_sf_result](https://www.gnu.org/software/gsl/doc/html/specfunc.html#c.gsl_sf_result) * *result*)

  These routines compute the complete elliptic integral ![\Pi(k,n)](https://www.gnu.org/software/gsl/doc/html/_images/math/92934b6dbee225af6069de22c2d1a3ba900252f3.png) to the accuracy specified by the mode variable [`mode`](https://www.gnu.org/software/gsl/doc/html/montecarlo.html#c.mode). Note that Abramowitz & Stegun define this function in terms of the parameters ![m = k^2](https://www.gnu.org/software/gsl/doc/html/_images/math/9d5667426cec24427e500eb97a3096029ad992e9.png) and ![\sin^2(\alpha) = k^2](https://www.gnu.org/software/gsl/doc/html/_images/math/ee65d5a743ecc0a7eb0abb32bc1ee23fe9dc1e50.png), with the change of sign ![n \to -n](https://www.gnu.org/software/gsl/doc/html/_images/math/bf297db77a13d0778e7247ba8a239a883dc1cd2d.png).

### Legendre Form of Incomplete Elliptic Integrals

- double `gsl_sf_ellint_F`(double *phi*, double *k*, [gsl_mode_t](https://www.gnu.org/software/gsl/doc/html/specfunc.html#c.gsl_mode_t) *mode*)

- int `gsl_sf_ellint_F_e`(double *phi*, double *k*, [gsl_mode_t](https://www.gnu.org/software/gsl/doc/html/specfunc.html#c.gsl_mode_t) *mode*, [gsl_sf_result](https://www.gnu.org/software/gsl/doc/html/specfunc.html#c.gsl_sf_result) * *result*)

  These routines compute the incomplete elliptic integral ![F(\phi,k)](https://www.gnu.org/software/gsl/doc/html/_images/math/5e447e3b2cba8230362e65c49302a2e312186a13.png) to the accuracy specified by the mode variable [`mode`](https://www.gnu.org/software/gsl/doc/html/montecarlo.html#c.mode). Note that Abramowitz & Stegun define this function in terms of the parameter ![m = k^2](https://www.gnu.org/software/gsl/doc/html/_images/math/9d5667426cec24427e500eb97a3096029ad992e9.png).

- double `gsl_sf_ellint_E`(double *phi*, double *k*, [gsl_mode_t](https://www.gnu.org/software/gsl/doc/html/specfunc.html#c.gsl_mode_t) *mode*)

- int `gsl_sf_ellint_E_e`(double *phi*, double *k*, [gsl_mode_t](https://www.gnu.org/software/gsl/doc/html/specfunc.html#c.gsl_mode_t) *mode*, [gsl_sf_result](https://www.gnu.org/software/gsl/doc/html/specfunc.html#c.gsl_sf_result) * *result*)

  These routines compute the incomplete elliptic integral ![E(\phi,k)](https://www.gnu.org/software/gsl/doc/html/_images/math/ee8224eb0d55eb0dab4e8d270ca1f764b2b9a677.png) to the accuracy specified by the mode variable [`mode`](https://www.gnu.org/software/gsl/doc/html/montecarlo.html#c.mode). Note that Abramowitz & Stegun define this function in terms of the parameter ![m = k^2](https://www.gnu.org/software/gsl/doc/html/_images/math/9d5667426cec24427e500eb97a3096029ad992e9.png).

- double `gsl_sf_ellint_P`(double *phi*, double *k*, double *n*, [gsl_mode_t](https://www.gnu.org/software/gsl/doc/html/specfunc.html#c.gsl_mode_t) *mode*)

- int `gsl_sf_ellint_P_e`(double *phi*, double *k*, double *n*, [gsl_mode_t](https://www.gnu.org/software/gsl/doc/html/specfunc.html#c.gsl_mode_t) *mode*, [gsl_sf_result](https://www.gnu.org/software/gsl/doc/html/specfunc.html#c.gsl_sf_result) * *result*)

  These routines compute the incomplete elliptic integral ![\Pi(\phi,k,n)](https://www.gnu.org/software/gsl/doc/html/_images/math/4684364912954c9c8aac26813ae3f2ee297e093c.png) to the accuracy specified by the mode variable [`mode`](https://www.gnu.org/software/gsl/doc/html/montecarlo.html#c.mode). Note that Abramowitz & Stegun define this function in terms of the parameters ![m = k^2](https://www.gnu.org/software/gsl/doc/html/_images/math/9d5667426cec24427e500eb97a3096029ad992e9.png) and ![\sin^2(\alpha) = k^2](https://www.gnu.org/software/gsl/doc/html/_images/math/ee65d5a743ecc0a7eb0abb32bc1ee23fe9dc1e50.png), with the change of sign ![n \to -n](https://www.gnu.org/software/gsl/doc/html/_images/math/bf297db77a13d0778e7247ba8a239a883dc1cd2d.png).

- double `gsl_sf_ellint_D`(double *phi*, double *k*, [gsl_mode_t](https://www.gnu.org/software/gsl/doc/html/specfunc.html#c.gsl_mode_t) *mode*)

- int `gsl_sf_ellint_D_e`(double *phi*, double *k*, [gsl_mode_t](https://www.gnu.org/software/gsl/doc/html/specfunc.html#c.gsl_mode_t) *mode*, [gsl_sf_result](https://www.gnu.org/software/gsl/doc/html/specfunc.html#c.gsl_sf_result) * *result*)

  These functions compute the incomplete elliptic integral ![D(\phi,k)](https://www.gnu.org/software/gsl/doc/html/_images/math/e8025c7151ab59ebe5480ab9887c33149de789ea.png) which is defined through the Carlson form ![RD(x,y,z)](https://www.gnu.org/software/gsl/doc/html/_images/math/d01cff478d89a62023a0c37832307a9b13d28aab.png) by the following relation,![D(\phi,k) = {1 \over 3} (\sin \phi)^3 RD (1-\sin^2(\phi), 1-k^2 \sin^2(\phi), 1)](https://www.gnu.org/software/gsl/doc/html/_images/math/409f3750cb4a632f778b2a2d423afd74d788e5ab.png)

### Carlson Forms

- double `gsl_sf_ellint_RC`(double *x*, double *y*, [gsl_mode_t](https://www.gnu.org/software/gsl/doc/html/specfunc.html#c.gsl_mode_t) *mode*)

- int `gsl_sf_ellint_RC_e`(double *x*, double *y*, [gsl_mode_t](https://www.gnu.org/software/gsl/doc/html/specfunc.html#c.gsl_mode_t) *mode*, [gsl_sf_result](https://www.gnu.org/software/gsl/doc/html/specfunc.html#c.gsl_sf_result) * *result*)

  These routines compute the incomplete elliptic integral ![RC(x,y)](https://www.gnu.org/software/gsl/doc/html/_images/math/1ba238bde832a3013a8245cd4968e3bdf8f31bfd.png) to the accuracy specified by the mode variable [`mode`](https://www.gnu.org/software/gsl/doc/html/montecarlo.html#c.mode).

- double `gsl_sf_ellint_RD`(double *x*, double *y*, double *z*, [gsl_mode_t](https://www.gnu.org/software/gsl/doc/html/specfunc.html#c.gsl_mode_t) *mode*)

- int `gsl_sf_ellint_RD_e`(double *x*, double *y*, double *z*, [gsl_mode_t](https://www.gnu.org/software/gsl/doc/html/specfunc.html#c.gsl_mode_t) *mode*, [gsl_sf_result](https://www.gnu.org/software/gsl/doc/html/specfunc.html#c.gsl_sf_result) * *result*)

  These routines compute the incomplete elliptic integral ![RD(x,y,z)](https://www.gnu.org/software/gsl/doc/html/_images/math/d01cff478d89a62023a0c37832307a9b13d28aab.png) to the accuracy specified by the mode variable [`mode`](https://www.gnu.org/software/gsl/doc/html/montecarlo.html#c.mode).

- double `gsl_sf_ellint_RF`(double *x*, double *y*, double *z*, [gsl_mode_t](https://www.gnu.org/software/gsl/doc/html/specfunc.html#c.gsl_mode_t) *mode*)

- int `gsl_sf_ellint_RF_e`(double *x*, double *y*, double *z*, [gsl_mode_t](https://www.gnu.org/software/gsl/doc/html/specfunc.html#c.gsl_mode_t) *mode*, [gsl_sf_result](https://www.gnu.org/software/gsl/doc/html/specfunc.html#c.gsl_sf_result) * *result*)

  These routines compute the incomplete elliptic integral ![RF(x,y,z)](https://www.gnu.org/software/gsl/doc/html/_images/math/73a62113b3068835188237682b0fc2a7a0961159.png) to the accuracy specified by the mode variable [`mode`](https://www.gnu.org/software/gsl/doc/html/montecarlo.html#c.mode).

- double `gsl_sf_ellint_RJ`(double *x*, double *y*, double *z*, double *p*, [gsl_mode_t](https://www.gnu.org/software/gsl/doc/html/specfunc.html#c.gsl_mode_t) *mode*)

- int `gsl_sf_ellint_RJ_e`(double *x*, double *y*, double *z*, double *p*, [gsl_mode_t](https://www.gnu.org/software/gsl/doc/html/specfunc.html#c.gsl_mode_t) *mode*, [gsl_sf_result](https://www.gnu.org/software/gsl/doc/html/specfunc.html#c.gsl_sf_result) * *result*)

  These routines compute the incomplete elliptic integral ![RJ(x,y,z,p)](https://www.gnu.org/software/gsl/doc/html/_images/math/2f3ddc98d528fff74d7369ecd655c003cb54a8be.png) to the accuracy specified by the mode variable [`mode`](https://www.gnu.org/software/gsl/doc/html/montecarlo.html#c.mode).

## Elliptic Functions (Jacobi)

The Jacobian Elliptic functions are defined in Abramowitz & Stegun, Chapter 16. The functions are declared in the header file `gsl_sf_elljac.h`.

- int `gsl_sf_elljac_e`(double *u*, double *m*, double * *sn*, double * *cn*, double * *dn*)

  This function computes the Jacobian elliptic functions ![sn(u|m)](https://www.gnu.org/software/gsl/doc/html/_images/math/3a5977e328983c8403b7c9e43bf776e7ba9eb9aa.png), ![cn(u|m)](https://www.gnu.org/software/gsl/doc/html/_images/math/e8a4c1140eb5666609298a3fdc32db0ea752c8db.png), ![dn(u|m)](https://www.gnu.org/software/gsl/doc/html/_images/math/e5153f7d14c4ca4f0fa59a3111f15437b892541f.png) by descending Landen transformations.

## Error Functions

The error function is described in Abramowitz & Stegun, Chapter 7. The functions in this section are declared in the header file `gsl_sf_erf.h`.

### Error Function

- double `gsl_sf_erf`(double *x*)

- int `gsl_sf_erf_e`(double *x*, [gsl_sf_result](https://www.gnu.org/software/gsl/doc/html/specfunc.html#c.gsl_sf_result) * *result*)

  These routines compute the error function ![\erf(x)](https://www.gnu.org/software/gsl/doc/html/_images/math/430aa8dd66cdd8ad745750e56c516ca8172c26b7.png), where ![\erf(x) = (2/\sqrt{\pi}) \int_0^x dt \exp(-t^2)](https://www.gnu.org/software/gsl/doc/html/_images/math/ef6a7bb90e4a83c7b49c050a8f07a41af16f4b17.png).

### Complementary Error Function

- double `gsl_sf_erfc`(double *x*)

- int `gsl_sf_erfc_e`(double *x*, [gsl_sf_result](https://www.gnu.org/software/gsl/doc/html/specfunc.html#c.gsl_sf_result) * *result*)

  These routines compute the complementary error function ![\erfc(x) = 1 - \erf(x) = (2/\sqrt{\pi}) \int_x^\infty \exp(-t^2)](https://www.gnu.org/software/gsl/doc/html/_images/math/3641d3f15a784faf7101a6612c2fc7bf7396bae3.png)

### Log Complementary Error Function

- double `gsl_sf_log_erfc`(double *x*)

- int `gsl_sf_log_erfc_e`(double *x*, [gsl_sf_result](https://www.gnu.org/software/gsl/doc/html/specfunc.html#c.gsl_sf_result) * *result*)

  These routines compute the logarithm of the complementary error function ![\log(\erfc(x))](https://www.gnu.org/software/gsl/doc/html/_images/math/fc23483520afb7e29540d9259d5a683dddd75f57.png).

### Probability functions

The probability functions for the Normal or Gaussian distribution are described in Abramowitz & Stegun, Section 26.2.

- double `gsl_sf_erf_Z`(double *x*)

- int `gsl_sf_erf_Z_e`(double *x*, [gsl_sf_result](https://www.gnu.org/software/gsl/doc/html/specfunc.html#c.gsl_sf_result) * *result*)

  These routines compute the Gaussian probability density function ![Z(x) = (1/\sqrt{2\pi}) \exp(-x^2/2)](https://www.gnu.org/software/gsl/doc/html/_images/math/a33fb125184c3ed19fda5137da0fc22fcb972722.png)

- double `gsl_sf_erf_Q`(double *x*)

- int `gsl_sf_erf_Q_e`(double *x*, [gsl_sf_result](https://www.gnu.org/software/gsl/doc/html/specfunc.html#c.gsl_sf_result) * *result*)

  These routines compute the upper tail of the Gaussian probability function ![Q(x) = (1/\sqrt{2\pi}) \int_x^\infty dt \exp(-t^2/2)](https://www.gnu.org/software/gsl/doc/html/_images/math/f4d52fa3bd54ae5355a89280f88c7d419c58adc8.png)

The *hazard function* for the normal distribution, also known as the inverse Mills’ ratio, is defined as,

![h(x) = {Z(x) \over Q(x)} = \sqrt{2 \over \pi} {\exp(-x^2 / 2) \over \erfc(x/\sqrt 2)}](https://www.gnu.org/software/gsl/doc/html/_images/math/c35b7e236351da5af4836c38d701b180b3c86efb.png)

It decreases rapidly as ![x](https://www.gnu.org/software/gsl/doc/html/_images/math/1c5715c8edd0bbc330baa0f0494f1f7951619754.png) approaches ![-\infty](https://www.gnu.org/software/gsl/doc/html/_images/math/342dc8da9165373fac8f4b9024eea6d6e76c53d9.png) and asymptotes to ![h(x) \sim x](https://www.gnu.org/software/gsl/doc/html/_images/math/9b89525fb4b8a42cdd02fe3072f3a5fc2b0e49a6.png) as ![x](https://www.gnu.org/software/gsl/doc/html/_images/math/1c5715c8edd0bbc330baa0f0494f1f7951619754.png) approaches ![+\infty](https://www.gnu.org/software/gsl/doc/html/_images/math/6be41aa884c6c25ee576afe5476eea358adc2931.png).

- double `gsl_sf_hazard`(double *x*)

- int `gsl_sf_hazard_e`(double *x*, [gsl_sf_result](https://www.gnu.org/software/gsl/doc/html/specfunc.html#c.gsl_sf_result) * *result*)

  These routines compute the hazard function for the normal distribution.

## Exponential Functions

The functions described in this section are declared in the header file `gsl_sf_exp.h`.

### Exponential Function

- double `gsl_sf_exp`(double *x*)

- int `gsl_sf_exp_e`(double *x*, [gsl_sf_result](https://www.gnu.org/software/gsl/doc/html/specfunc.html#c.gsl_sf_result) * *result*)

  These routines provide an exponential function ![\exp(x)](https://www.gnu.org/software/gsl/doc/html/_images/math/bccd4c1b6a65e230c3bdd85f36186dd29b47c1fe.png) using GSL semantics and error checking.

- int `gsl_sf_exp_e10_e`(double *x*, [gsl_sf_result_e10](https://www.gnu.org/software/gsl/doc/html/specfunc.html#c.gsl_sf_result_e10) * *result*)

  This function computes the exponential ![\exp(x)](https://www.gnu.org/software/gsl/doc/html/_images/math/bccd4c1b6a65e230c3bdd85f36186dd29b47c1fe.png) using the [`gsl_sf_result_e10`](https://www.gnu.org/software/gsl/doc/html/specfunc.html#c.gsl_sf_result_e10) type to return a result with extended range. This function may be useful if the value of ![\exp(x)](https://www.gnu.org/software/gsl/doc/html/_images/math/bccd4c1b6a65e230c3bdd85f36186dd29b47c1fe.png) would overflow the numeric range of `double`.

- double `gsl_sf_exp_mult`(double *x*, double *y*)

- int `gsl_sf_exp_mult_e`(double *x*, double *y*, [gsl_sf_result](https://www.gnu.org/software/gsl/doc/html/specfunc.html#c.gsl_sf_result) * *result*)

  These routines exponentiate `x` and multiply by the factor `y` to return the product ![y \exp(x)](https://www.gnu.org/software/gsl/doc/html/_images/math/e09333a5a2a24858c7d4a4623708d8cbc34fae03.png).

- int `gsl_sf_exp_mult_e10_e`(const double *x*, const double *y*, [gsl_sf_result_e10](https://www.gnu.org/software/gsl/doc/html/specfunc.html#c.gsl_sf_result_e10) * *result*)

  This function computes the product ![y \exp(x)](https://www.gnu.org/software/gsl/doc/html/_images/math/e09333a5a2a24858c7d4a4623708d8cbc34fae03.png) using the [`gsl_sf_result_e10`](https://www.gnu.org/software/gsl/doc/html/specfunc.html#c.gsl_sf_result_e10) type to return a result with extended numeric range.

### Relative Exponential Functions

- double `gsl_sf_expm1`(double *x*)

- int `gsl_sf_expm1_e`(double *x*, [gsl_sf_result](https://www.gnu.org/software/gsl/doc/html/specfunc.html#c.gsl_sf_result) * *result*)

  These routines compute the quantity ![\exp(x)-1](https://www.gnu.org/software/gsl/doc/html/_images/math/40d15d7dfac571e38c04c5f83ac3cb15c43ecff5.png) using an algorithm that is accurate for small ![x](https://www.gnu.org/software/gsl/doc/html/_images/math/1c5715c8edd0bbc330baa0f0494f1f7951619754.png).

- double `gsl_sf_exprel`(double *x*)

- int `gsl_sf_exprel_e`(double *x*, [gsl_sf_result](https://www.gnu.org/software/gsl/doc/html/specfunc.html#c.gsl_sf_result) * *result*)

  These routines compute the quantity ![(\exp(x)-1)/x](https://www.gnu.org/software/gsl/doc/html/_images/math/9202e39663c549b9a8af08539bf365ede7f03cbd.png) using an algorithm that is accurate for small `x`. For small `x` the algorithm is based on the expansion ![(\exp(x)-1)/x = 1 + x/2 + x^2/(2*3) + x^3/(2*3*4) + \dots](https://www.gnu.org/software/gsl/doc/html/_images/math/2e0466cb1f5f07bddab6e9525b630c5ae0409f14.png).

- double `gsl_sf_exprel_2`(double *x*)

- int `gsl_sf_exprel_2_e`(double *x*, [gsl_sf_result](https://www.gnu.org/software/gsl/doc/html/specfunc.html#c.gsl_sf_result) * *result*)

  These routines compute the quantity ![2(\exp(x)-1-x)/x^2](https://www.gnu.org/software/gsl/doc/html/_images/math/e531b57d5a1b11f8eaa7933759f8c8e187d7c3a3.png) using an algorithm that is accurate for small `x`. For small `x` the algorithm is based on the expansion ![2(\exp(x)-1-x)/x^2 = 1 + x/3 + x^2/(3*4) + x^3/(3*4*5) + \dots](https://www.gnu.org/software/gsl/doc/html/_images/math/0729d090c90338b5a0667a84a49c5f36412dc1df.png).

- double `gsl_sf_exprel_n`(int *n*, double *x*)

- int `gsl_sf_exprel_n_e`(int *n*, double *x*, [gsl_sf_result](https://www.gnu.org/software/gsl/doc/html/specfunc.html#c.gsl_sf_result) * *result*)

  These routines compute the ![N](https://www.gnu.org/software/gsl/doc/html/_images/math/99f8a3c97c19212ca3eee7d00dddc8d68b2d1c27.png)-relative exponential, which is the `n`-th generalization of the functions [`gsl_sf_exprel()`](https://www.gnu.org/software/gsl/doc/html/specfunc.html#c.gsl_sf_exprel) and [`gsl_sf_exprel_2()`](https://www.gnu.org/software/gsl/doc/html/specfunc.html#c.gsl_sf_exprel_2). The ![N](https://www.gnu.org/software/gsl/doc/html/_images/math/99f8a3c97c19212ca3eee7d00dddc8d68b2d1c27.png)-relative exponential is given by,![\hbox{exprel}_N(x)             &= N!/x^N \left(\exp(x) - \sum_{k=0}^{N-1} x^k/k!\right)\cr             &= 1 + x/(N+1) + x^2/((N+1)(N+2)) + \dots\cr             &= {}_1F_1(1,1+N,x)\cr](https://www.gnu.org/software/gsl/doc/html/_images/math/77bffc0ba109094f73abd82439982d16c57638e9.png)

### Exponentiation With Error Estimate

- int `gsl_sf_exp_err_e`(double *x*, double *dx*, [gsl_sf_result](https://www.gnu.org/software/gsl/doc/html/specfunc.html#c.gsl_sf_result) * *result*)

  This function exponentiates `x` with an associated absolute error `dx`.

- int `gsl_sf_exp_err_e10_e`(double *x*, double *dx*, [gsl_sf_result_e10](https://www.gnu.org/software/gsl/doc/html/specfunc.html#c.gsl_sf_result_e10) * *result*)

  This function exponentiates a quantity `x` with an associated absolute error `dx` using the [`gsl_sf_result_e10`](https://www.gnu.org/software/gsl/doc/html/specfunc.html#c.gsl_sf_result_e10) type to return a result with extended range.

- int `gsl_sf_exp_mult_err_e`(double *x*, double *dx*, double *y*, double *dy*, [gsl_sf_result](https://www.gnu.org/software/gsl/doc/html/specfunc.html#c.gsl_sf_result) * *result*)

  This routine computes the product ![y \exp(x)](https://www.gnu.org/software/gsl/doc/html/_images/math/e09333a5a2a24858c7d4a4623708d8cbc34fae03.png) for the quantities `x`, `y` with associated absolute errors `dx`, `dy`.

- int `gsl_sf_exp_mult_err_e10_e`(double *x*, double *dx*, double *y*, double *dy*, [gsl_sf_result_e10](https://www.gnu.org/software/gsl/doc/html/specfunc.html#c.gsl_sf_result_e10) * *result*)

  This routine computes the product ![y \exp(x)](https://www.gnu.org/software/gsl/doc/html/_images/math/e09333a5a2a24858c7d4a4623708d8cbc34fae03.png) for the quantities `x`, `y` with associated absolute errors `dx`, `dy` using the [`gsl_sf_result_e10`](https://www.gnu.org/software/gsl/doc/html/specfunc.html#c.gsl_sf_result_e10) type to return a result with extended range.

## Exponential Integrals

Information on the exponential integrals can be found in Abramowitz & Stegun, Chapter 5. These functions are declared in the header file `gsl_sf_expint.h`.

### Exponential Integral



- double `gsl_sf_expint_E1`(double *x*)

- int `gsl_sf_expint_E1_e`(double *x*, [gsl_sf_result](https://www.gnu.org/software/gsl/doc/html/specfunc.html#c.gsl_sf_result) * *result*)

  These routines compute the exponential integral ![E_1(x)](https://www.gnu.org/software/gsl/doc/html/_images/math/2b6fa2c7d905c3c0c5165a390108fc52e8bd0489.png),![E_1(x) := \Re \int_1^\infty dt \exp(-xt)/t.](https://www.gnu.org/software/gsl/doc/html/_images/math/3f765852dc0cd637c04f6235eb0c53b4a72b5662.png)

- double `gsl_sf_expint_E2`(double *x*)

- int `gsl_sf_expint_E2_e`(double *x*, [gsl_sf_result](https://www.gnu.org/software/gsl/doc/html/specfunc.html#c.gsl_sf_result) * *result*)

  These routines compute the second-order exponential integral ![E_2(x)](https://www.gnu.org/software/gsl/doc/html/_images/math/555735f76a0655c8254096f25cab2da776a8d4e8.png),![E_2(x) := \Re \int_1^\infty dt \exp(-xt)/t^2](https://www.gnu.org/software/gsl/doc/html/_images/math/451511996f2130559525b79b311d1dc11627513b.png)

- double `gsl_sf_expint_En`(int *n*, double *x*)

- int `gsl_sf_expint_En_e`(int *n*, double *x*, [gsl_sf_result](https://www.gnu.org/software/gsl/doc/html/specfunc.html#c.gsl_sf_result) * *result*)

  These routines compute the exponential integral ![E_n(x)](https://www.gnu.org/software/gsl/doc/html/_images/math/73f7e9821936d6faa492edfb8aa70ef9bd58545c.png) of order `n`,![E_n(x) := \Re \int_1^\infty dt \exp(-xt)/t^n.](https://www.gnu.org/software/gsl/doc/html/_images/math/95e0a233346bb489ff1c6a6f21039b4dd01bedb0.png)

### Ei(x)

- double `gsl_sf_expint_Ei`(double *x*)

- int `gsl_sf_expint_Ei_e`(double *x*, [gsl_sf_result](https://www.gnu.org/software/gsl/doc/html/specfunc.html#c.gsl_sf_result) * *result*)

  These routines compute the exponential integral ![Ei(x)](https://www.gnu.org/software/gsl/doc/html/_images/math/2e62669695a3b5dd0585f90d46873493dd0cfa35.png),![\hbox{Ei}(x) = - PV \left( \int_{-x}^\infty dt \exp(-t)/t \right)](https://www.gnu.org/software/gsl/doc/html/_images/math/9fb544911cfdb1ba9e7bf5fabcb4e86746333a59.png)where ![PV](https://www.gnu.org/software/gsl/doc/html/_images/math/ed5a0eb45496e8fab0012c5e22ed5a26caa85318.png) denotes the principal value of the integral.

### Hyperbolic Integrals



- double `gsl_sf_Shi`(double *x*)

- int `gsl_sf_Shi_e`(double *x*, [gsl_sf_result](https://www.gnu.org/software/gsl/doc/html/specfunc.html#c.gsl_sf_result) * *result*)

  These routines compute the integral![\hbox{Shi}(x) = \int_0^x dt \sinh(t)/t](https://www.gnu.org/software/gsl/doc/html/_images/math/1b2adce6b024da29ad0009770ab3dc9967942b60.png)

- double `gsl_sf_Chi`(double *x*)

- int `gsl_sf_Chi_e`(double *x*, [gsl_sf_result](https://www.gnu.org/software/gsl/doc/html/specfunc.html#c.gsl_sf_result) * *result*)

  These routines compute the integral![\hbox{Chi}(x) := \Re \left[ \gamma_E + \log(x) + \int_0^x dt (\cosh(t)-1)/t \right]](https://www.gnu.org/software/gsl/doc/html/_images/math/8eea6c8cf7de8a706e7db9245903ea2f5edd2853.png)where ![\gamma_E](https://www.gnu.org/software/gsl/doc/html/_images/math/ec64dfe5360b14f3f96df7880960225416b68d06.png) is the Euler constant (available as the macro `M_EULER`).

### Ei_3(x)

- double `gsl_sf_expint_3`(double *x*)

- int `gsl_sf_expint_3_e`(double *x*, [gsl_sf_result](https://www.gnu.org/software/gsl/doc/html/specfunc.html#c.gsl_sf_result) * *result*)

  These routines compute the third-order exponential integral![{\rm Ei}_3(x) = \int_0^x dt \exp(-t^3)](https://www.gnu.org/software/gsl/doc/html/_images/math/e6d7154ac46c0cfef5f44a17390f244ed26289a5.png)for ![x \ge 0](https://www.gnu.org/software/gsl/doc/html/_images/math/fe351bceb81b6a8053c3711170955c252c9b2fc6.png).

### Trigonometric Integrals



- double `gsl_sf_Si`(const double *x*)

- int `gsl_sf_Si_e`(double *x*, [gsl_sf_result](https://www.gnu.org/software/gsl/doc/html/specfunc.html#c.gsl_sf_result) * *result*)

  These routines compute the Sine integral![\hbox{Si}(x) = \int_0^x dt \sin(t)/t](https://www.gnu.org/software/gsl/doc/html/_images/math/68457761c34349a932250bbf233c33640bddbbc3.png)

- double `gsl_sf_Ci`(const double *x*)

- int `gsl_sf_Ci_e`(double *x*, [gsl_sf_result](https://www.gnu.org/software/gsl/doc/html/specfunc.html#c.gsl_sf_result) * *result*)

  These routines compute the Cosine integral![\hbox{Ci}(x) = -\int_x^\infty dt \cos(t)/t](https://www.gnu.org/software/gsl/doc/html/_images/math/f565b1e777ea586ff1e00a976121ff14333f31d6.png)for ![x > 0](https://www.gnu.org/software/gsl/doc/html/_images/math/fe20eedd0aef40a657ad04835ff74fa8a73efb56.png)

### Arctangent Integral



- double `gsl_sf_atanint`(double *x*)

- int `gsl_sf_atanint_e`(double *x*, [gsl_sf_result](https://www.gnu.org/software/gsl/doc/html/specfunc.html#c.gsl_sf_result) * *result*)

  These routines compute the Arctangent integral, which is defined as![\hbox{AtanInt}(x) = \int_0^x dt \arctan(t)/t](https://www.gnu.org/software/gsl/doc/html/_images/math/8e7b7e5107decee158f23da6244d79723046b363.png)

## Fermi-Dirac Function

The functions described in this section are declared in the header file `gsl_sf_fermi_dirac.h`.

### Complete Fermi-Dirac Integrals

The complete Fermi-Dirac integral ![F_j(x)](https://www.gnu.org/software/gsl/doc/html/_images/math/a55973e8cd63e64ac2ef6ca0914105b73df180da.png) is given by,

![F_j(x) := {1\over\Gamma(j+1)} \int_0^\infty dt {t^j  \over (\exp(t-x) + 1)}](https://www.gnu.org/software/gsl/doc/html/_images/math/45f67ae999cf45c454c3351c7259b9698be6552e.png)

Note that the Fermi-Dirac integral is sometimes defined without the normalisation factor in other texts.

- double `gsl_sf_fermi_dirac_m1`(double *x*)

- int `gsl_sf_fermi_dirac_m1_e`(double *x*, [gsl_sf_result](https://www.gnu.org/software/gsl/doc/html/specfunc.html#c.gsl_sf_result) * *result*)

  These routines compute the complete Fermi-Dirac integral with an index of ![-1](https://www.gnu.org/software/gsl/doc/html/_images/math/e26aecc9561355dd0e6fbf2742d0ff0d6b58c864.png). This integral is given by ![F_{-1}(x) = e^x / (1 + e^x)](https://www.gnu.org/software/gsl/doc/html/_images/math/f04ade6cb2a3e543c47529cfa7e679f0bfc317c0.png).

- double `gsl_sf_fermi_dirac_0`(double *x*)

- int `gsl_sf_fermi_dirac_0_e`(double *x*, [gsl_sf_result](https://www.gnu.org/software/gsl/doc/html/specfunc.html#c.gsl_sf_result) * *result*)

  These routines compute the complete Fermi-Dirac integral with an index of ![0](https://www.gnu.org/software/gsl/doc/html/_images/math/3b9bed1b0ccefe4f2e5e9be7ff710ff847c91a5c.png). This integral is given by ![F_0(x) = \ln(1 + e^x)](https://www.gnu.org/software/gsl/doc/html/_images/math/1c0ed5f8dd437859af14f0bc740a344e23e910c4.png).

- double `gsl_sf_fermi_dirac_1`(double *x*)

- int `gsl_sf_fermi_dirac_1_e`(double *x*, [gsl_sf_result](https://www.gnu.org/software/gsl/doc/html/specfunc.html#c.gsl_sf_result) * *result*)

  These routines compute the complete Fermi-Dirac integral with an index of ![1](https://www.gnu.org/software/gsl/doc/html/_images/math/585581da152cce359927e7b9eeadef5ca7c3a42d.png), ![F_1(x) = \int_0^\infty dt (t /(\exp(t-x)+1))](https://www.gnu.org/software/gsl/doc/html/_images/math/de47889ee88ee3c7b0fa59c3994d3874ecc83f9f.png).

- double `gsl_sf_fermi_dirac_2`(double *x*)

- int `gsl_sf_fermi_dirac_2_e`(double *x*, [gsl_sf_result](https://www.gnu.org/software/gsl/doc/html/specfunc.html#c.gsl_sf_result) * *result*)

  These routines compute the complete Fermi-Dirac integral with an index of ![2](https://www.gnu.org/software/gsl/doc/html/_images/math/29cce9750a86cad81782e7abce5a772e0d0ea393.png), ![F_2(x) = (1/2) \int_0^\infty dt (t^2 /(\exp(t-x)+1))](https://www.gnu.org/software/gsl/doc/html/_images/math/c600259f5e95bcb5e843cbfc1896a7e05fbb7414.png).

- double `gsl_sf_fermi_dirac_int`(int *j*, double *x*)

- int `gsl_sf_fermi_dirac_int_e`(int *j*, double *x*, [gsl_sf_result](https://www.gnu.org/software/gsl/doc/html/specfunc.html#c.gsl_sf_result) * *result*)

  These routines compute the complete Fermi-Dirac integral with an integer index of ![j](https://www.gnu.org/software/gsl/doc/html/_images/math/7bc2ad053edeecad39b867d60495ac2caa4555dc.png), ![F_j(x) = (1/\Gamma(j+1)) \int_0^\infty dt (t^j /(\exp(t-x)+1))](https://www.gnu.org/software/gsl/doc/html/_images/math/98a7d033e3dfa46968392f9f4a0fcfe33c8576fc.png).

- double `gsl_sf_fermi_dirac_mhalf`(double *x*)

- int `gsl_sf_fermi_dirac_mhalf_e`(double *x*, [gsl_sf_result](https://www.gnu.org/software/gsl/doc/html/specfunc.html#c.gsl_sf_result) * *result*)

  These routines compute the complete Fermi-Dirac integral ![F_{-1/2}(x)](https://www.gnu.org/software/gsl/doc/html/_images/math/f6db8629c4cb8d04b3829b5462b916db76048611.png).

- double `gsl_sf_fermi_dirac_half`(double *x*)

- int `gsl_sf_fermi_dirac_half_e`(double *x*, [gsl_sf_result](https://www.gnu.org/software/gsl/doc/html/specfunc.html#c.gsl_sf_result) * *result*)

  These routines compute the complete Fermi-Dirac integral ![F_{1/2}(x)](https://www.gnu.org/software/gsl/doc/html/_images/math/8677216c776cd49eca6367c6bf50b2472a9714f1.png).

- double `gsl_sf_fermi_dirac_3half`(double *x*)

- int `gsl_sf_fermi_dirac_3half_e`(double *x*, [gsl_sf_result](https://www.gnu.org/software/gsl/doc/html/specfunc.html#c.gsl_sf_result) * *result*)

  These routines compute the complete Fermi-Dirac integral ![F_{3/2}(x)](https://www.gnu.org/software/gsl/doc/html/_images/math/2347c2663b821f1e6b3fbab9a0940e5150c388a6.png).

### Incomplete Fermi-Dirac Integrals

The incomplete Fermi-Dirac integral ![F_j(x,b)](https://www.gnu.org/software/gsl/doc/html/_images/math/437354e92a3a544a46e3b79d9dff093a4c96f22a.png) is given by,

![F_j(x,b) := {1\over\Gamma(j+1)} \int_b^\infty dt {t^j  \over (\exp(t-x) + 1)}](https://www.gnu.org/software/gsl/doc/html/_images/math/b97e5f87df02d8f1f16b15a341139bc476a166c3.png)

- double `gsl_sf_fermi_dirac_inc_0`(double *x*, double *b*)

- int `gsl_sf_fermi_dirac_inc_0_e`(double *x*, double *b*, [gsl_sf_result](https://www.gnu.org/software/gsl/doc/html/specfunc.html#c.gsl_sf_result) * *result*)

  These routines compute the incomplete Fermi-Dirac integral with an index of zero, ![F_0(x,b) = \ln(1 + e^{b-x}) - (b-x)](https://www.gnu.org/software/gsl/doc/html/_images/math/9273ccbc165533418dcc588106affba8884cb9ee.png)

## Gamma and Beta Functions

The following routines compute the gamma and beta functions in their full and incomplete forms, as well as various kinds of factorials. The functions described in this section are declared in the header file `gsl_sf_gamma.h`.

### Gamma Functions

The Gamma function is defined by the following integral,

![\Gamma(x) = \int_0^{\infty} dt t^{x-1} \exp(-t)](https://www.gnu.org/software/gsl/doc/html/_images/math/c9889dd894a5626246a18a33e6f5547a835bbed2.png)

It is related to the factorial function by ![\Gamma(n) = (n-1)!](https://www.gnu.org/software/gsl/doc/html/_images/math/647035237ce454f6a955ec874559b2b9c3758c27.png) for positive integer ![n](https://www.gnu.org/software/gsl/doc/html/_images/math/a24554f1502cf204e7aa24a6c064962c3504de48.png). Further information on the Gamma function can be found in Abramowitz & Stegun, Chapter 6.

- double `gsl_sf_gamma`(double *x*)

- int `gsl_sf_gamma_e`(double *x*, [gsl_sf_result](https://www.gnu.org/software/gsl/doc/html/specfunc.html#c.gsl_sf_result) * *result*)

  These routines compute the Gamma function ![\Gamma(x)](https://www.gnu.org/software/gsl/doc/html/_images/math/96679e2506d89c9a4317f56fda713283cfb96964.png), subject to ![x](https://www.gnu.org/software/gsl/doc/html/_images/math/1c5715c8edd0bbc330baa0f0494f1f7951619754.png) not being a negative integer or zero. The function is computed using the real Lanczos method. The maximum value of ![x](https://www.gnu.org/software/gsl/doc/html/_images/math/1c5715c8edd0bbc330baa0f0494f1f7951619754.png) such that ![\Gamma(x)](https://www.gnu.org/software/gsl/doc/html/_images/math/96679e2506d89c9a4317f56fda713283cfb96964.png) is not considered an overflow is given by the macro `GSL_SF_GAMMA_XMAX` and is 171.0.



- double `gsl_sf_lngamma`(double *x*)

- int `gsl_sf_lngamma_e`(double *x*, [gsl_sf_result](https://www.gnu.org/software/gsl/doc/html/specfunc.html#c.gsl_sf_result) * *result*)

  These routines compute the logarithm of the Gamma function, ![\log(\Gamma(x))](https://www.gnu.org/software/gsl/doc/html/_images/math/e2618f080a4b2fe5fe4ab73c2e4ab80aacfc2ab4.png), subject to ![x](https://www.gnu.org/software/gsl/doc/html/_images/math/1c5715c8edd0bbc330baa0f0494f1f7951619754.png) not being a negative integer or zero. For ![x < 0](https://www.gnu.org/software/gsl/doc/html/_images/math/93f2a3ce4b90fb07a7dc88c1bc60a6a43c93f45b.png) the real part of ![\log(\Gamma(x))](https://www.gnu.org/software/gsl/doc/html/_images/math/e2618f080a4b2fe5fe4ab73c2e4ab80aacfc2ab4.png) is returned, which is equivalent to ![\log(|\Gamma(x)|)](https://www.gnu.org/software/gsl/doc/html/_images/math/2de54976b4acc44372738df1426c6da95bfd6ba1.png). The function is computed using the real Lanczos method.

- int `gsl_sf_lngamma_sgn_e`(double *x*, [gsl_sf_result](https://www.gnu.org/software/gsl/doc/html/specfunc.html#c.gsl_sf_result) * *result_lg*, double * *sgn*)

  This routine computes the sign of the gamma function and the logarithm of its magnitude, subject to ![x](https://www.gnu.org/software/gsl/doc/html/_images/math/1c5715c8edd0bbc330baa0f0494f1f7951619754.png) not being a negative integer or zero. The function is computed using the real Lanczos method. The value of the gamma function and its error can be reconstructed using the relation ![\Gamma(x) = sgn * \exp(result\_lg)](https://www.gnu.org/software/gsl/doc/html/_images/math/90d47125bfbc4056eea56291fbb7791e0f8aea0b.png), taking into account the two components of `result_lg`.



- double `gsl_sf_gammastar`(double *x*)

- int `gsl_sf_gammastar_e`(double *x*, [gsl_sf_result](https://www.gnu.org/software/gsl/doc/html/specfunc.html#c.gsl_sf_result) * *result*)

  These routines compute the regulated Gamma Function ![\Gamma^*(x)](https://www.gnu.org/software/gsl/doc/html/_images/math/99db1b5b85a51ce48be941959f6242715590964f.png) for ![x > 0](https://www.gnu.org/software/gsl/doc/html/_images/math/fe20eedd0aef40a657ad04835ff74fa8a73efb56.png). The regulated gamma function is given by,![\Gamma^*(x) &= \Gamma(x)/(\sqrt{2\pi} x^{(x-1/2)} \exp(-x))\cr             &= \left(1 + {1 \over 12x} + ...\right) \quad\hbox{for~} x\to \infty\cr](https://www.gnu.org/software/gsl/doc/html/_images/math/2b84f8cc3777c46ca3c37e335e692ceed7a9c807.png)and is a useful suggestion of Temme.



- double `gsl_sf_gammainv`(double *x*)

- int `gsl_sf_gammainv_e`(double *x*, [gsl_sf_result](https://www.gnu.org/software/gsl/doc/html/specfunc.html#c.gsl_sf_result) * *result*)

  These routines compute the reciprocal of the gamma function, ![1/\Gamma(x)](https://www.gnu.org/software/gsl/doc/html/_images/math/77ddc26223bd6fe31ecd896d4c302633eb03bdf8.png) using the real Lanczos method.



- int `gsl_sf_lngamma_complex_e`(double *zr*, double *zi*, [gsl_sf_result](https://www.gnu.org/software/gsl/doc/html/specfunc.html#c.gsl_sf_result) * *lnr*, [gsl_sf_result](https://www.gnu.org/software/gsl/doc/html/specfunc.html#c.gsl_sf_result) * *arg*)

  This routine computes ![\log(\Gamma(z))](https://www.gnu.org/software/gsl/doc/html/_images/math/f9f05c861560fd9b23b208c0ab4b05b6c8538510.png) for complex ![z = z_r + i z_i](https://www.gnu.org/software/gsl/doc/html/_images/math/2bd1785b0de8324fc35b8e3a6810476232e4ccf1.png) and ![z](https://www.gnu.org/software/gsl/doc/html/_images/math/4c9e578e25d40c5b0c5cad91a1e2f34e58b5be5e.png) not a negative integer or zero, using the complex Lanczos method. The returned parameters are ![lnr = \log|\Gamma(z)|](https://www.gnu.org/software/gsl/doc/html/_images/math/e91d15034237b4471c65157be18b5e50c59c35a7.png) and ![arg = \arg(\Gamma(z))](https://www.gnu.org/software/gsl/doc/html/_images/math/c8b3cc2d1f7dcc1fdb6e1c64c16527f98e488c7e.png) in ![(-\pi,\pi]](https://www.gnu.org/software/gsl/doc/html/_images/math/2c4cbfed45ff565dab518bf9e6dcdb1558490157.png). Note that the phase part (`arg`) is not well-determined when ![|z|](https://www.gnu.org/software/gsl/doc/html/_images/math/e21e85b701e7cd55df38706041140e741b3ab735.png) is very large, due to inevitable roundoff in restricting to ![(-\pi,\pi]](https://www.gnu.org/software/gsl/doc/html/_images/math/2c4cbfed45ff565dab518bf9e6dcdb1558490157.png). This will result in a `GSL_ELOSS`error when it occurs. The absolute value part (`lnr`), however, never suffers from loss of precision.

### Factorials

Although factorials can be computed from the Gamma function, using the relation ![n! = \Gamma(n+1)](https://www.gnu.org/software/gsl/doc/html/_images/math/5ee202ecd4ec1ae2d79850d8c339ff28308779ee.png)for non-negative integer ![n](https://www.gnu.org/software/gsl/doc/html/_images/math/a24554f1502cf204e7aa24a6c064962c3504de48.png), it is usually more efficient to call the functions in this section, particularly for small values of ![n](https://www.gnu.org/software/gsl/doc/html/_images/math/a24554f1502cf204e7aa24a6c064962c3504de48.png), whose factorial values are maintained in hardcoded tables.



- double `gsl_sf_fact`(unsigned int *n*)

- int `gsl_sf_fact_e`(unsigned int *n*, [gsl_sf_result](https://www.gnu.org/software/gsl/doc/html/specfunc.html#c.gsl_sf_result) * *result*)

  These routines compute the factorial ![n!](https://www.gnu.org/software/gsl/doc/html/_images/math/1ad3a44f5e6828b120d23481384a0c7683f0c666.png). The factorial is related to the Gamma function by ![n! = \Gamma(n+1)](https://www.gnu.org/software/gsl/doc/html/_images/math/5ee202ecd4ec1ae2d79850d8c339ff28308779ee.png). The maximum value of ![n](https://www.gnu.org/software/gsl/doc/html/_images/math/a24554f1502cf204e7aa24a6c064962c3504de48.png) such that ![n!](https://www.gnu.org/software/gsl/doc/html/_images/math/1ad3a44f5e6828b120d23481384a0c7683f0c666.png) is not considered an overflow is given by the macro `GSL_SF_FACT_NMAX` and is 170.



- double `gsl_sf_doublefact`(unsigned int *n*)

- int `gsl_sf_doublefact_e`(unsigned int *n*, [gsl_sf_result](https://www.gnu.org/software/gsl/doc/html/specfunc.html#c.gsl_sf_result) * *result*)

  These routines compute the double factorial ![n!! = n(n-2)(n-4) \dots](https://www.gnu.org/software/gsl/doc/html/_images/math/af7ab3e818aec89cf10336d929423394a53e04d8.png). The maximum value of ![n](https://www.gnu.org/software/gsl/doc/html/_images/math/a24554f1502cf204e7aa24a6c064962c3504de48.png) such that ![n!!](https://www.gnu.org/software/gsl/doc/html/_images/math/a8b0d69aab7775ea873ba5e6e25d97023b56f120.png) is not considered an overflow is given by the macro `GSL_SF_DOUBLEFACT_NMAX` and is 297.



- double `gsl_sf_lnfact`(unsigned int *n*)

- int `gsl_sf_lnfact_e`(unsigned int *n*, [gsl_sf_result](https://www.gnu.org/software/gsl/doc/html/specfunc.html#c.gsl_sf_result) * *result*)

  These routines compute the logarithm of the factorial of `n`, ![\log(n!)](https://www.gnu.org/software/gsl/doc/html/_images/math/605920e333e5f1f2ba0af33f856b51fe6d96962c.png). The algorithm is faster than computing ![\ln(\Gamma(n+1))](https://www.gnu.org/software/gsl/doc/html/_images/math/5ae62107df9fddf3ead9f40a3ce5ec6b81544d84.png) via [`gsl_sf_lngamma()`](https://www.gnu.org/software/gsl/doc/html/specfunc.html#c.gsl_sf_lngamma) for ![n < 170](https://www.gnu.org/software/gsl/doc/html/_images/math/409d30a4d79d3fadc040ada4ccd59d2d46f22d0c.png), but defers for larger `n`.



- double `gsl_sf_lndoublefact`(unsigned int *n*)

- int `gsl_sf_lndoublefact_e`(unsigned int *n*, [gsl_sf_result](https://www.gnu.org/software/gsl/doc/html/specfunc.html#c.gsl_sf_result) * *result*)

  These routines compute the logarithm of the double factorial of `n`, ![\log(n!!)](https://www.gnu.org/software/gsl/doc/html/_images/math/bcc47b5cf535160c43b338df4551b9064cdd3360.png).



- double `gsl_sf_choose`(unsigned int *n*, unsigned int *m*)

- int `gsl_sf_choose_e`(unsigned int *n*, unsigned int *m*, [gsl_sf_result](https://www.gnu.org/software/gsl/doc/html/specfunc.html#c.gsl_sf_result) * *result*)

  These routines compute the combinatorial factor `n choose m` ![= n!/(m!(n-m)!)](https://www.gnu.org/software/gsl/doc/html/_images/math/62470c4f451736a5c7c10beef395298e37bd00ab.png)



- double `gsl_sf_lnchoose`(unsigned int *n*, unsigned int *m*)

- int `gsl_sf_lnchoose_e`(unsigned int *n*, unsigned int *m*, [gsl_sf_result](https://www.gnu.org/software/gsl/doc/html/specfunc.html#c.gsl_sf_result) * *result*)

  These routines compute the logarithm of `n choose m`. This is equivalent to the sum ![\log(n!) - \log(m!) - \log((n-m)!)](https://www.gnu.org/software/gsl/doc/html/_images/math/1d40ee54f05bda4ab03812a7898d7de5300d8716.png).



- double `gsl_sf_taylorcoeff`(int *n*, double *x*)

- int `gsl_sf_taylorcoeff_e`(int *n*, double *x*, [gsl_sf_result](https://www.gnu.org/software/gsl/doc/html/specfunc.html#c.gsl_sf_result) * *result*)

  These routines compute the Taylor coefficient ![x^n / n!](https://www.gnu.org/software/gsl/doc/html/_images/math/afd2ff50c5c58890fc8876d1aea4589a562279d6.png) for ![x \ge 0](https://www.gnu.org/software/gsl/doc/html/_images/math/fe351bceb81b6a8053c3711170955c252c9b2fc6.png), ![n \ge 0](https://www.gnu.org/software/gsl/doc/html/_images/math/a4a2e96a00e9bf71dda87e79f4f76f7ad0bca7cf.png)



### Pochhammer Symbol



- double `gsl_sf_poch`(double *a*, double *x*)

- int `gsl_sf_poch_e`(double *a*, double *x*, [gsl_sf_result](https://www.gnu.org/software/gsl/doc/html/specfunc.html#c.gsl_sf_result) * *result*)

  These routines compute the Pochhammer symbol ![(a)_x = \Gamma(a + x)/\Gamma(a)](https://www.gnu.org/software/gsl/doc/html/_images/math/43533fcbddab6b7d1878fc90bbd2a40f601ebcea.png). The Pochhammer symbol is also known as the Apell symbol and sometimes written as ![(a,x)](https://www.gnu.org/software/gsl/doc/html/_images/math/617841550dd061524387c96dc7bec7067912881b.png). When ![a](https://www.gnu.org/software/gsl/doc/html/_images/math/122c953b5f811e9257bc7845f38503ebe1dff7de.png) and ![a + x](https://www.gnu.org/software/gsl/doc/html/_images/math/e41ca8134ada3ce7685069b1aa4e4003f4fdd6e8.png)are negative integers or zero, the limiting value of the ratio is returned.



- double `gsl_sf_lnpoch`(double *a*, double *x*)

- int `gsl_sf_lnpoch_e`(double *a*, double *x*, [gsl_sf_result](https://www.gnu.org/software/gsl/doc/html/specfunc.html#c.gsl_sf_result) * *result*)

  These routines compute the logarithm of the Pochhammer symbol, ![\log((a)_x) = \log(\Gamma(a + x)/\Gamma(a))](https://www.gnu.org/software/gsl/doc/html/_images/math/1fec710d8b7e6ce89eeba0c58c432f0d92cdf90c.png).

- int `gsl_sf_lnpoch_sgn_e`(double *a*, double *x*, [gsl_sf_result](https://www.gnu.org/software/gsl/doc/html/specfunc.html#c.gsl_sf_result) * *result*, double * *sgn*)

  These routines compute the sign of the Pochhammer symbol and the logarithm of its magnitude. The computed parameters are ![result = \log(|(a)_x|)](https://www.gnu.org/software/gsl/doc/html/_images/math/64c19b08278e7038f064fc19f68b8b582bf8b1de.png) with a corresponding error term, and ![sgn = \sgn((a)_x)](https://www.gnu.org/software/gsl/doc/html/_images/math/7560d861dc9e49f0114e813cae7781baba5fbe42.png) where ![(a)_x = \Gamma(a + x)/\Gamma(a)](https://www.gnu.org/software/gsl/doc/html/_images/math/43533fcbddab6b7d1878fc90bbd2a40f601ebcea.png).



- double `gsl_sf_pochrel`(double *a*, double *x*)

- int `gsl_sf_pochrel_e`(double *a*, double *x*, [gsl_sf_result](https://www.gnu.org/software/gsl/doc/html/specfunc.html#c.gsl_sf_result) * *result*)

  These routines compute the relative Pochhammer symbol ![((a)_x - 1)/x](https://www.gnu.org/software/gsl/doc/html/_images/math/80aa3b6fc74427de7771e08b1f81b194f3c70bfa.png) where ![(a)_x = \Gamma(a + x)/\Gamma(a)](https://www.gnu.org/software/gsl/doc/html/_images/math/43533fcbddab6b7d1878fc90bbd2a40f601ebcea.png).

### Incomplete Gamma Functions



- double `gsl_sf_gamma_inc`(double *a*, double *x*)

- int `gsl_sf_gamma_inc_e`(double *a*, double *x*, [gsl_sf_result](https://www.gnu.org/software/gsl/doc/html/specfunc.html#c.gsl_sf_result) * *result*)

  These functions compute the unnormalized incomplete Gamma Function ![\Gamma(a,x) = \int_x^\infty dt t^{(a-1)} \exp(-t)](https://www.gnu.org/software/gsl/doc/html/_images/math/6c244269dd7abddab7fe8c727f56ce46b423bc4b.png) for ![a](https://www.gnu.org/software/gsl/doc/html/_images/math/122c953b5f811e9257bc7845f38503ebe1dff7de.png) real and ![x \ge 0](https://www.gnu.org/software/gsl/doc/html/_images/math/fe351bceb81b6a8053c3711170955c252c9b2fc6.png).



- double `gsl_sf_gamma_inc_Q`(double *a*, double *x*)

- int `gsl_sf_gamma_inc_Q_e`(double *a*, double *x*, [gsl_sf_result](https://www.gnu.org/software/gsl/doc/html/specfunc.html#c.gsl_sf_result) * *result*)

  These routines compute the normalized incomplete Gamma Function ![Q(a,x) = 1/\Gamma(a) \int_x^\infty dt t^{(a-1)} \exp(-t)](https://www.gnu.org/software/gsl/doc/html/_images/math/f1c79a6a0999db4f98da724a7d3e0032be8e852d.png) for ![a > 0](https://www.gnu.org/software/gsl/doc/html/_images/math/681843f11d32939bc8a6fc88fe2837f905790fdb.png), ![x \ge 0](https://www.gnu.org/software/gsl/doc/html/_images/math/fe351bceb81b6a8053c3711170955c252c9b2fc6.png).



- double `gsl_sf_gamma_inc_P`(double *a*, double *x*)

- int `gsl_sf_gamma_inc_P_e`(double *a*, double *x*, [gsl_sf_result](https://www.gnu.org/software/gsl/doc/html/specfunc.html#c.gsl_sf_result) * *result*)

  These routines compute the complementary normalized incomplete Gamma Function ![P(a,x) = 1 - Q(a,x) = 1/\Gamma(a) \int_0^x dt t^{(a-1)} \exp(-t)](https://www.gnu.org/software/gsl/doc/html/_images/math/e2d9297ed030343f09d3caf1865759b28638a524.png) for ![a > 0](https://www.gnu.org/software/gsl/doc/html/_images/math/681843f11d32939bc8a6fc88fe2837f905790fdb.png), ![x \ge 0](https://www.gnu.org/software/gsl/doc/html/_images/math/fe351bceb81b6a8053c3711170955c252c9b2fc6.png).Note that Abramowitz & Stegun call ![P(a,x)](https://www.gnu.org/software/gsl/doc/html/_images/math/57d176e639dd52ef2f2cbdd74f1ed1696929434d.png) the incomplete gamma function (section 6.5).

### Beta Functions



- double `gsl_sf_beta`(double *a*, double *b*)

- int `gsl_sf_beta_e`(double *a*, double *b*, [gsl_sf_result](https://www.gnu.org/software/gsl/doc/html/specfunc.html#c.gsl_sf_result) * *result*)

  These routines compute the Beta Function, ![B(a,b) = \Gamma(a)\Gamma(b)/\Gamma(a+b)](https://www.gnu.org/software/gsl/doc/html/_images/math/9916da87786a0b630e1e6bbcb50b1c91b96b0f2f.png) subject to ![a](https://www.gnu.org/software/gsl/doc/html/_images/math/122c953b5f811e9257bc7845f38503ebe1dff7de.png) and ![b](https://www.gnu.org/software/gsl/doc/html/_images/math/88082f6488798774b9cd95a59c1ab4c6cafe909f.png)not being negative integers.



- double `gsl_sf_lnbeta`(double *a*, double *b*)

- int `gsl_sf_lnbeta_e`(double *a*, double *b*, [gsl_sf_result](https://www.gnu.org/software/gsl/doc/html/specfunc.html#c.gsl_sf_result) * *result*)

  These routines compute the logarithm of the Beta Function, ![\log(B(a,b))](https://www.gnu.org/software/gsl/doc/html/_images/math/60e817f57a02ac3013f9bbcfaa69a2b2ac5d1cdd.png) subject to ![a](https://www.gnu.org/software/gsl/doc/html/_images/math/122c953b5f811e9257bc7845f38503ebe1dff7de.png) and ![b](https://www.gnu.org/software/gsl/doc/html/_images/math/88082f6488798774b9cd95a59c1ab4c6cafe909f.png) not being negative integers.

### Incomplete Beta Function



- double `gsl_sf_beta_inc`(double *a*, double *b*, double *x*)

- int `gsl_sf_beta_inc_e`(double *a*, double *b*, double *x*, [gsl_sf_result](https://www.gnu.org/software/gsl/doc/html/specfunc.html#c.gsl_sf_result) * *result*)

  These routines compute the normalized incomplete Beta function ![I_x(a,b) = B_x(a,b) / B(a,b)](https://www.gnu.org/software/gsl/doc/html/_images/math/991f609168faf21f44c83deb87b519aafd872c5b.png)where![B_x(a,b) = \int_0^x t^{a-1} (1-t)^{b-1} dt](https://www.gnu.org/software/gsl/doc/html/_images/math/6edbf3ab45214ce6f9a3c2a98d977481623fce92.png)for ![0 \le x \le 1](https://www.gnu.org/software/gsl/doc/html/_images/math/fae6247d5d755eecbdf8fe9b611a68154f5c53ce.png). For ![a > 0](https://www.gnu.org/software/gsl/doc/html/_images/math/681843f11d32939bc8a6fc88fe2837f905790fdb.png), ![b > 0](https://www.gnu.org/software/gsl/doc/html/_images/math/22c2c2831e7a825fdc2fd2a9909a75acb01e9e4c.png) the value is computed using a continued fraction expansion. For all other values it is computed using the relation![I_x(a,b,x) = (1/a) x^a {}_2F_1(a,1-b,a+1,x)/B(a,b)](https://www.gnu.org/software/gsl/doc/html/_images/math/56ec3fa7135800f0cad1ccc477d9eabba8cdc870.png)

## Gegenbauer Functions

The Gegenbauer polynomials are defined in Abramowitz & Stegun, Chapter 22, where they are known as Ultraspherical polynomials. The functions described in this section are declared in the header file `gsl_sf_gegenbauer.h`.

- double `gsl_sf_gegenpoly_1`(double *lambda*, double *x*)

- double `gsl_sf_gegenpoly_2`(double *lambda*, double *x*)

- double `gsl_sf_gegenpoly_3`(double *lambda*, double *x*)

- int `gsl_sf_gegenpoly_1_e`(double *lambda*, double *x*, [gsl_sf_result](https://www.gnu.org/software/gsl/doc/html/specfunc.html#c.gsl_sf_result) * *result*)

- int `gsl_sf_gegenpoly_2_e`(double *lambda*, double *x*, [gsl_sf_result](https://www.gnu.org/software/gsl/doc/html/specfunc.html#c.gsl_sf_result) * *result*)

- int `gsl_sf_gegenpoly_3_e`(double *lambda*, double *x*, [gsl_sf_result](https://www.gnu.org/software/gsl/doc/html/specfunc.html#c.gsl_sf_result) * *result*)

  These functions evaluate the Gegenbauer polynomials ![C^{(\lambda)}_n(x)](https://www.gnu.org/software/gsl/doc/html/_images/math/3f631227668a9277f08b14767079225cf07eff7a.png) using explicit representations for ![n = 1, 2, 3](https://www.gnu.org/software/gsl/doc/html/_images/math/834b59e608f6bc5b7aa5248a1f83102ecbff39a0.png).

- double `gsl_sf_gegenpoly_n`(int *n*, double *lambda*, double *x*)

- int `gsl_sf_gegenpoly_n_e`(int *n*, double *lambda*, double *x*, [gsl_sf_result](https://www.gnu.org/software/gsl/doc/html/specfunc.html#c.gsl_sf_result) * *result*)

  These functions evaluate the Gegenbauer polynomial ![C^{(\lambda)}_n(x)](https://www.gnu.org/software/gsl/doc/html/_images/math/3f631227668a9277f08b14767079225cf07eff7a.png) for a specific value of `n`,`lambda`, `x` subject to ![\lambda > -1/2](https://www.gnu.org/software/gsl/doc/html/_images/math/9b90b73b20d8808fdd66a1fa5e91bb8d574495fd.png), ![n \ge 0](https://www.gnu.org/software/gsl/doc/html/_images/math/a4a2e96a00e9bf71dda87e79f4f76f7ad0bca7cf.png).

- int `gsl_sf_gegenpoly_array`(int *nmax*, double *lambda*, double *x*, double *result_array[]*)

  This function computes an array of Gegenbauer polynomials ![C^{(\lambda)}_n(x)](https://www.gnu.org/software/gsl/doc/html/_images/math/3f631227668a9277f08b14767079225cf07eff7a.png) for ![n = 0, 1, 2, \dots, nmax](https://www.gnu.org/software/gsl/doc/html/_images/math/654537ff15049b2844d807f893f70eb2084fb924.png), subject to ![\lambda > -1/2](https://www.gnu.org/software/gsl/doc/html/_images/math/9b90b73b20d8808fdd66a1fa5e91bb8d574495fd.png), ![nmax \ge 0](https://www.gnu.org/software/gsl/doc/html/_images/math/9231981f7fa9fa6385f065e7efa05013b2651585.png).

## Hermite Polynomials and Functions



The Hermite polynomials exist in two variants: the probabilists’ version ![He_n(x)](https://www.gnu.org/software/gsl/doc/html/_images/math/7ceaa11993ae7e55e30ad5da40a8d7247ffa4241.png) and the physicists’version ![H_n(x)](https://www.gnu.org/software/gsl/doc/html/_images/math/2ba97b7fe8fff12fa2e0a8ef0c4a1b31d13d402a.png). The are defined by the derivatives

![He_n(x) & = (-1)^n e^{x^2/2} \left({d \over dx}\right)^n e^{-x^2/2} \\ H_n(x) & = (-1)^n e^{x^2} \left({d \over dx}\right)^n e^{-x^2}](https://www.gnu.org/software/gsl/doc/html/_images/math/3eb9765274fa576c909be9d7f95e7967fc748908.png)

They are connected via

![He_n(x) & = 2^{-n/2} H_n \left( {x \over \sqrt{2}} \right) \\ H_n(x) & = 2^{n/2} He_n \left( \sqrt{2} x \right)](https://www.gnu.org/software/gsl/doc/html/_images/math/7be0ff56c56246a8c1acda4403a8e47607fdf6e9.png)

and satisfy the ordinary differential equations

![He_n^{\prime\prime}(x) - x He_n^{\prime}(x) + n He_n(x) & = 0 \\ H_n^{\prime\prime}(x) - 2x H_n^{\prime}(x) + 2n H_n(x) & = 0](https://www.gnu.org/software/gsl/doc/html/_images/math/d00c44464436c49e25118b07fe7d1e8ef6dc3134.png)

The closely related Hermite functions are defined by

![\psi_n(x) = \left( n! \sqrt{\pi} \right)^{-1/2} e^{-x^2/2} He_n \left( {\sqrt{2} x} \right)](https://www.gnu.org/software/gsl/doc/html/_images/math/cde63db973ef2f8be4731c43ce2d53a1ec9c3df8.png)

and satisfy the Schrödinger equation for a quantum mechanical harmonic oscillator

![\psi_n^{\prime\prime}(x) + (2n + 1 - x^2) \psi_n(x) = 0](https://www.gnu.org/software/gsl/doc/html/_images/math/d383b29ab33841d808c0427352dac054af9d5918.png)

Maybe most importantly, the Hermite functions ![\psi_n](https://www.gnu.org/software/gsl/doc/html/_images/math/b0e6c88df143828ae93965b2159475e090f3405d.png) are eigenfunctions of the (continuous) Fourier transform.

For further information see Abramowitz & Stegun, Chapter 22 and Szego, Gabor (1939, 1958, 1967), Orthogonal Polynomials, American Mathematical Society. The Hermite polynomials and functions are defined in the header file `gsl_sf_hermite.h`.

### Hermite Polynomials

- double `gsl_sf_hermite_prob`(const int *n*, const double *x*)

- int `gsl_sf_hermite_prob_e`(const int *n*, const double *x*, [gsl_sf_result](https://www.gnu.org/software/gsl/doc/html/specfunc.html#c.gsl_sf_result) * *result*)

  These routines evaluate the probabilists’ Hermite polynomial ![He_n(x)](https://www.gnu.org/software/gsl/doc/html/_images/math/7ceaa11993ae7e55e30ad5da40a8d7247ffa4241.png) of order `n` at position `x`.

- int `gsl_sf_hermite_prob_array`(const int *nmax*, const double *x*, double * *result_array*)

  This routine evaluates all probabilists’ Hermite polynomials ![He_n(x)](https://www.gnu.org/software/gsl/doc/html/_images/math/7ceaa11993ae7e55e30ad5da40a8d7247ffa4241.png) up to order `nmax` at position `x`. The results are stored in `result_array`.

- double `gsl_sf_hermite_prob_series`(const int *n*, const double *x*, const double * *a*)

- int `gsl_sf_hermite_prob_series_e`(const int *n*, const double *x*, const double * *a*, [gsl_sf_result](https://www.gnu.org/software/gsl/doc/html/specfunc.html#c.gsl_sf_result) * *result*)

  These routines evaluate the series ![\sum_{j=0}^n a_j He_j(x)](https://www.gnu.org/software/gsl/doc/html/_images/math/4fd664c186b20321f0b27587497b68968e254c80.png) with ![He_j](https://www.gnu.org/software/gsl/doc/html/_images/math/b6d127f7b018cd0ffa7269df010dc9c6c419dd46.png) being the ![j](https://www.gnu.org/software/gsl/doc/html/_images/math/7bc2ad053edeecad39b867d60495ac2caa4555dc.png)-th probabilists’ Hermite polynomial using the Clenshaw algorithm.

- double `gsl_sf_hermite_phys`(const int *n*, const double *x*)

- int `gsl_sf_hermite_phys_e`(const int *n*, const double *x*, [gsl_sf_result](https://www.gnu.org/software/gsl/doc/html/specfunc.html#c.gsl_sf_result) * *result*)

  These routines evaluate the physicists’ Hermite polynomial ![H_n(x)](https://www.gnu.org/software/gsl/doc/html/_images/math/2ba97b7fe8fff12fa2e0a8ef0c4a1b31d13d402a.png) of order `n` at position `x`.

- int `gsl_sf_hermite_phys_array`(const int *nmax*, const double *x*, double * *result_array*)

  This routine evaluates all physicists’ Hermite polynomials ![H_n](https://www.gnu.org/software/gsl/doc/html/_images/math/bbfd664cff8f4ec8447ae555c98138c1028cacc9.png) up to order `nmax` at position `x`. The results are stored in `result_array`.

- double `gsl_sf_hermite_phys_series`(const int *n*, const double *x*, const double * *a*)

- int `gsl_sf_hermite_phys_series_e`(const int *n*, const double *x*, const double * *a*, [gsl_sf_result](https://www.gnu.org/software/gsl/doc/html/specfunc.html#c.gsl_sf_result) * *result*)

  These routines evaluate the series ![\sum_{j=0}^n a_j H_j(x)](https://www.gnu.org/software/gsl/doc/html/_images/math/10a7bb619b7a385a17fd233d3defc6c1dc3b702d.png) with ![H_j](https://www.gnu.org/software/gsl/doc/html/_images/math/f8e8b91e44f16e182b4fb4afa26f9f307889587a.png) being the ![j](https://www.gnu.org/software/gsl/doc/html/_images/math/7bc2ad053edeecad39b867d60495ac2caa4555dc.png)-th physicists’ Hermite polynomial using the Clenshaw algorithm.

### Hermite Functions

- double `gsl_sf_hermite_func`(const int *n*, const double *x*)

- int `gsl_sf_hermite_func_e`(const int *n*, const double *x*, [gsl_sf_result](https://www.gnu.org/software/gsl/doc/html/specfunc.html#c.gsl_sf_result) * *result*)

  These routines evaluate the Hermite function ![\psi_n(x)](https://www.gnu.org/software/gsl/doc/html/_images/math/505a2deb7fe44f2507598b55784e54142d7cfc21.png) of order `n` at position `x`.

- int `gsl_sf_hermite_func_array`(const int *nmax*, const double *x*, double * *result_array*)

  This routine evaluates all Hermite functions ![\psi_n(x)](https://www.gnu.org/software/gsl/doc/html/_images/math/505a2deb7fe44f2507598b55784e54142d7cfc21.png) up to order `nmax` at position `x`. The results are stored in `result_array`.

- double `gsl_sf_hermite_func_series`(const int *n*, const double *x*, const double * *a*)

- int `gsl_sf_hermite_func_series_e`(const int *n*, const double *x*, const double * *a*, [gsl_sf_result](https://www.gnu.org/software/gsl/doc/html/specfunc.html#c.gsl_sf_result) * *result*)

  These routines evaluate the series ![\sum_{j=0}^n a_j \psi_j(x)](https://www.gnu.org/software/gsl/doc/html/_images/math/e03ad8b971290fd9332ce91c2d7c7069ab8e69d0.png) with ![\psi_j](https://www.gnu.org/software/gsl/doc/html/_images/math/38997fa2a068fe1129b129760a2685de41e67751.png) being the ![j](https://www.gnu.org/software/gsl/doc/html/_images/math/7bc2ad053edeecad39b867d60495ac2caa4555dc.png)-th Hermite function using the Clenshaw algorithm.

### Derivatives of Hermite Polynomials



- double `gsl_sf_hermite_prob_der`(const int *m*, const int *n*, const double *x*)

- int `gsl_sf_hermite_prob_der_e`(const int *m*, const int *n*, const double *x*, [gsl_sf_result](https://www.gnu.org/software/gsl/doc/html/specfunc.html#c.gsl_sf_result) * *result*)

  These routines evaluate the `m`-th derivative of the probabilists’ Hermite polynomial ![He_n(x)](https://www.gnu.org/software/gsl/doc/html/_images/math/7ceaa11993ae7e55e30ad5da40a8d7247ffa4241.png) of order `n` at position `x`.

- int `gsl_sf_hermite_prob_array_der`(const int *m*, const int *nmax*, const double *x*, double * *result_array*)

  This routine evaluates the `m`-th derivative of all probabilists’ Hermite polynomials ![He_n(x)](https://www.gnu.org/software/gsl/doc/html/_images/math/7ceaa11993ae7e55e30ad5da40a8d7247ffa4241.png) up to order `nmax` at position `x`. The results are stored in `result_array`.

- int `gsl_sf_hermite_prob_der_array`(const int *mmax*, const int *n*, const double *x*, double * *result_array*)

  This routine evaluates all derivatives (starting from 0) up to the `mmax`-th derivative of the probabilists’ Hermite polynomial of order `n` ![He_n(x)](https://www.gnu.org/software/gsl/doc/html/_images/math/7ceaa11993ae7e55e30ad5da40a8d7247ffa4241.png) at position `x`. The results are stored in `result_array`.

- double `gsl_sf_hermite_phys_der`(const int *m*, const int *n*, const double *x*)

- int `gsl_sf_hermite_phys_der_e`(const int *m*, const int *n*, const double *x*, [gsl_sf_result](https://www.gnu.org/software/gsl/doc/html/specfunc.html#c.gsl_sf_result) * *result*)

  These routines evaluate the `m`-th derivative of the physicists’ Hermite polynomial ![H_n(x)](https://www.gnu.org/software/gsl/doc/html/_images/math/2ba97b7fe8fff12fa2e0a8ef0c4a1b31d13d402a.png) of order `n` at position `x`.

- int `gsl_sf_hermite_phys_array_der`(const int *m*, const int *nmax*, const double *x*, double * *result_array*)

  This routine evaluates the `m`-th derivative of all physicists’ Hermite polynomials ![H_n](https://www.gnu.org/software/gsl/doc/html/_images/math/bbfd664cff8f4ec8447ae555c98138c1028cacc9.png) up to order `nmax` at position `x`. The results are stored in `result_array`.

- int `gsl_sf_hermite_phys_der_array`(const int *mmax*, const int *n*, const double *x*, double * *result_array*)

  This routine evaluates all derivatives (starting from 0) up to the `mmax`-th derivative of the physicists’ Hermite polynomial of order `n` ![H_n](https://www.gnu.org/software/gsl/doc/html/_images/math/bbfd664cff8f4ec8447ae555c98138c1028cacc9.png) at position `x`. The results are stored in `result_array`.

### Derivatives of Hermite Functions



- double `gsl_sf_hermite_func_der`(const int *m*, const int *n*, const double *x*)

- int `gsl_sf_hermite_func_der_e`(const int *m*, const int *n*, const double *x*, [gsl_sf_result](https://www.gnu.org/software/gsl/doc/html/specfunc.html#c.gsl_sf_result) * *result*)

  These routines evaluate the `m`-th derivative of the Hermite function ![\psi_n(x)](https://www.gnu.org/software/gsl/doc/html/_images/math/505a2deb7fe44f2507598b55784e54142d7cfc21.png) of order `n` at position `x`.

### Zeros of Hermite Polynomials and Hermite Functions

These routines calculate the ![s](https://www.gnu.org/software/gsl/doc/html/_images/math/dd69d21b40087e08fc6827d2f3481f57462046e7.png)-th zero of the Hermite Polynomial/Function of order ![n](https://www.gnu.org/software/gsl/doc/html/_images/math/a24554f1502cf204e7aa24a6c064962c3504de48.png). Since the zeros are symmetrical around zero, only positive zeros are calculated, ordered from smallest to largest, starting from index 1. Only for odd polynomial orders a zeroth zero exists, its value always being zero.

- double `gsl_sf_hermite_prob_zero`(const int *n*, const int *s*)

- int `gsl_sf_hermite_prob_zero_e`(const int *n*, const int *s*, [gsl_sf_result](https://www.gnu.org/software/gsl/doc/html/specfunc.html#c.gsl_sf_result) * *result*)

  These routines evaluate the `s`-th zero of the probabilists’ Hermite polynomial ![He_n(x)](https://www.gnu.org/software/gsl/doc/html/_images/math/7ceaa11993ae7e55e30ad5da40a8d7247ffa4241.png) of order `n`.

- double `gsl_sf_hermite_phys_zero`(const int *n*, const int *s*)

- int `gsl_sf_hermite_phys_zero_e`(const int *n*, const int *s*, [gsl_sf_result](https://www.gnu.org/software/gsl/doc/html/specfunc.html#c.gsl_sf_result) * *result*)

  These routines evaluate the `s`-th zero of the physicists’ Hermite polynomial ![H_n(x)](https://www.gnu.org/software/gsl/doc/html/_images/math/2ba97b7fe8fff12fa2e0a8ef0c4a1b31d13d402a.png) of order `n`.

- double `gsl_sf_hermite_func_zero`(const int *n*, const int *s*)

- int `gsl_sf_hermite_func_zero_e`(const int *n*, const int *s*, [gsl_sf_result](https://www.gnu.org/software/gsl/doc/html/specfunc.html#c.gsl_sf_result) * *result*)

  These routines evaluate the `s`-th zero of the Hermite function ![\psi_n(x)](https://www.gnu.org/software/gsl/doc/html/_images/math/505a2deb7fe44f2507598b55784e54142d7cfc21.png) of order `n`.

## Hypergeometric Functions

Hypergeometric functions are described in Abramowitz & Stegun, Chapters 13 and 15. These functions are declared in the header file `gsl_sf_hyperg.h`.

- double `gsl_sf_hyperg_0F1`(double *c*, double *x*)

- int `gsl_sf_hyperg_0F1_e`(double *c*, double *x*, [gsl_sf_result](https://www.gnu.org/software/gsl/doc/html/specfunc.html#c.gsl_sf_result) * *result*)

  These routines compute the hypergeometric function![{}_0F_1(c,x)](https://www.gnu.org/software/gsl/doc/html/_images/math/bcfa3e7f7d2880925b66a51d6f82af02a5da02ac.png)

- double `gsl_sf_hyperg_1F1_int`(int *m*, int *n*, double *x*)

- int `gsl_sf_hyperg_1F1_int_e`(int *m*, int *n*, double *x*, [gsl_sf_result](https://www.gnu.org/software/gsl/doc/html/specfunc.html#c.gsl_sf_result) * *result*)

  These routines compute the confluent hypergeometric function![{}_1F_1(m,n,x) = M(m,n,x)](https://www.gnu.org/software/gsl/doc/html/_images/math/38b515b1dec65a516e6ae1530ff61e623d00b380.png)for integer parameters `m`, `n`.

- double `gsl_sf_hyperg_1F1`(double *a*, double *b*, double *x*)

- int `gsl_sf_hyperg_1F1_e`(double *a*, double *b*, double *x*, [gsl_sf_result](https://www.gnu.org/software/gsl/doc/html/specfunc.html#c.gsl_sf_result) * *result*)

  These routines compute the confluent hypergeometric function![{}_1F_1(a,b,x) = M(a,b,x)](https://www.gnu.org/software/gsl/doc/html/_images/math/3494b70234133f38abf4af6555cf1bc404ccf524.png)for general parameters `a`, `b`.

- double `gsl_sf_hyperg_U_int`(int *m*, int *n*, double *x*)

- int `gsl_sf_hyperg_U_int_e`(int *m*, int *n*, double *x*, [gsl_sf_result](https://www.gnu.org/software/gsl/doc/html/specfunc.html#c.gsl_sf_result) * *result*)

  These routines compute the confluent hypergeometric function ![U(m,n,x)](https://www.gnu.org/software/gsl/doc/html/_images/math/7f824e6ec456b7cfbbb83600ea8e4777174fbec7.png) for integer parameters `m`, `n`.

- int `gsl_sf_hyperg_U_int_e10_e`(int *m*, int *n*, double *x*, [gsl_sf_result_e10](https://www.gnu.org/software/gsl/doc/html/specfunc.html#c.gsl_sf_result_e10) * *result*)

  This routine computes the confluent hypergeometric function ![U(m,n,x)](https://www.gnu.org/software/gsl/doc/html/_images/math/7f824e6ec456b7cfbbb83600ea8e4777174fbec7.png) for integer parameters `m`, `n` using the [`gsl_sf_result_e10`](https://www.gnu.org/software/gsl/doc/html/specfunc.html#c.gsl_sf_result_e10) type to return a result with extended range.

- double `gsl_sf_hyperg_U`(double *a*, double *b*, double *x*)

- int `gsl_sf_hyperg_U_e`(double *a*, double *b*, double *x*, [gsl_sf_result](https://www.gnu.org/software/gsl/doc/html/specfunc.html#c.gsl_sf_result) * *result*)

  These routines compute the confluent hypergeometric function ![U(a,b,x)](https://www.gnu.org/software/gsl/doc/html/_images/math/4b6600c24dc0c908af429b9a054adced72253f43.png).

- int `gsl_sf_hyperg_U_e10_e`(double *a*, double *b*, double *x*, [gsl_sf_result_e10](https://www.gnu.org/software/gsl/doc/html/specfunc.html#c.gsl_sf_result_e10) * *result*)

  This routine computes the confluent hypergeometric function ![U(a,b,x)](https://www.gnu.org/software/gsl/doc/html/_images/math/4b6600c24dc0c908af429b9a054adced72253f43.png) using the [`gsl_sf_result_e10`](https://www.gnu.org/software/gsl/doc/html/specfunc.html#c.gsl_sf_result_e10) type to return a result with extended range.

- double `gsl_sf_hyperg_2F1`(double *a*, double *b*, double *c*, double *x*)

- int `gsl_sf_hyperg_2F1_e`(double *a*, double *b*, double *c*, double *x*, [gsl_sf_result](https://www.gnu.org/software/gsl/doc/html/specfunc.html#c.gsl_sf_result) * *result*)

  These routines compute the Gauss hypergeometric function![{}_2F_1(a,b,c,x) = F(a,b,c,x)](https://www.gnu.org/software/gsl/doc/html/_images/math/9fc3a67c57f5d97b801a3ba0872461a3d612d85d.png)for ![|x| < 1](https://www.gnu.org/software/gsl/doc/html/_images/math/596e6c7f75eb9f289d30a08bc5e911aa74197ee0.png). If the arguments ![(a,b,c,x)](https://www.gnu.org/software/gsl/doc/html/_images/math/8e98137bc8bd71df98b8896a1e89e327660dc7bf.png) are too close to a singularity then the function can return the error code `GSL_EMAXITER` when the series approximation converges too slowly. This occurs in the region of ![x = 1](https://www.gnu.org/software/gsl/doc/html/_images/math/6e6e07a91c43bfbcc5dc9f231bf866e89faed555.png), ![c - a - b = m](https://www.gnu.org/software/gsl/doc/html/_images/math/af3360a176fef23821fdb97537ce8412b8490f54.png) for integer m.

- double `gsl_sf_hyperg_2F1_conj`(double *aR*, double *aI*, double *c*, double *x*)

- int `gsl_sf_hyperg_2F1_conj_e`(double *aR*, double *aI*, double *c*, double *x*, [gsl_sf_result](https://www.gnu.org/software/gsl/doc/html/specfunc.html#c.gsl_sf_result) * *result*)

  These routines compute the Gauss hypergeometric function![{}_2F_1(a_R + i a_I, aR - i aI, c, x)](https://www.gnu.org/software/gsl/doc/html/_images/math/7af6e83dc0305d2a95a455c7a719a6924dc70ce2.png)with complex parameters for ![|x| < 1](https://www.gnu.org/software/gsl/doc/html/_images/math/596e6c7f75eb9f289d30a08bc5e911aa74197ee0.png).

- double `gsl_sf_hyperg_2F1_renorm`(double *a*, double *b*, double *c*, double *x*)

- int `gsl_sf_hyperg_2F1_renorm_e`(double *a*, double *b*, double *c*, double *x*, [gsl_sf_result](https://www.gnu.org/software/gsl/doc/html/specfunc.html#c.gsl_sf_result) * *result*)

  These routines compute the renormalized Gauss hypergeometric function![{}_2F_1(a,b,c,x) / \Gamma(c)](https://www.gnu.org/software/gsl/doc/html/_images/math/8d80e62a24a2a2693fca56210af69981289189d7.png)for ![|x| < 1](https://www.gnu.org/software/gsl/doc/html/_images/math/596e6c7f75eb9f289d30a08bc5e911aa74197ee0.png).

- double `gsl_sf_hyperg_2F1_conj_renorm`(double *aR*, double *aI*, double *c*, double *x*)

- int `gsl_sf_hyperg_2F1_conj_renorm_e`(double *aR*, double *aI*, double *c*, double *x*, [gsl_sf_result](https://www.gnu.org/software/gsl/doc/html/specfunc.html#c.gsl_sf_result) * *result*)

  These routines compute the renormalized Gauss hypergeometric function![{}_2F_1(a_R + i a_I, a_R - i a_I, c, x) / \Gamma(c)](https://www.gnu.org/software/gsl/doc/html/_images/math/a9de4e4e748295af84ba331e5c33616e9a17fec8.png)for ![|x| < 1](https://www.gnu.org/software/gsl/doc/html/_images/math/596e6c7f75eb9f289d30a08bc5e911aa74197ee0.png).

- double `gsl_sf_hyperg_2F0`(double *a*, double *b*, double *x*)

- int `gsl_sf_hyperg_2F0_e`(double *a*, double *b*, double *x*, [gsl_sf_result](https://www.gnu.org/software/gsl/doc/html/specfunc.html#c.gsl_sf_result) * *result*)

  These routines compute the hypergeometric function![{}_2F_0(a,b,x)](https://www.gnu.org/software/gsl/doc/html/_images/math/8a2a943739dec2b03810e9d694b8dd4ed1502843.png)The series representation is a divergent hypergeometric series. However, for ![x < 0](https://www.gnu.org/software/gsl/doc/html/_images/math/93f2a3ce4b90fb07a7dc88c1bc60a6a43c93f45b.png) we have![{}_2F_0(a,b,x) = (-1/x)^a U(a,1+a-b,-1/x)](https://www.gnu.org/software/gsl/doc/html/_images/math/abbddf2698d1491a1fc04695c41257cfc65996e6.png)



## Laguerre Functions

The generalized Laguerre polynomials, sometimes referred to as associated Laguerre polynomials, are defined in terms of confluent hypergeometric functions as

![L^a_n(x) = {(a+1)_n \over n!} {}_1F_1(-n,a+1,x)](https://www.gnu.org/software/gsl/doc/html/_images/math/8d44a60a537879047cdeb4f2f1f924f7bfe2b199.png)

where ![(a)_n](https://www.gnu.org/software/gsl/doc/html/_images/math/5b7a65612977e37ac37178547ce8f38d10af8147.png) is the [Pochhammer symbol](https://www.gnu.org/software/gsl/doc/html/specfunc.html#pochhammer-symbol) (rising factorial). They are related to the plain Laguerre polynomials ![L_n(x)](https://www.gnu.org/software/gsl/doc/html/_images/math/b123242dbc6bc4d9cfd0c886281f55e895f544f9.png) by ![L^0_n(x) = L_n(x)](https://www.gnu.org/software/gsl/doc/html/_images/math/8a87c6ecb8b250d0f9ba106bcb67929864ab5af1.png) and ![L^k_n(x) = (-1)^k (d^k/dx^k) L_{(n+k)}(x)](https://www.gnu.org/software/gsl/doc/html/_images/math/ad26c13692d6830eae52c4f03d2a66c2bbe438ba.png) For more information see Abramowitz & Stegun, Chapter 22.

The functions described in this section are declared in the header file `gsl_sf_laguerre.h`.

- double `gsl_sf_laguerre_1`(double *a*, double *x*)

- double `gsl_sf_laguerre_2`(double *a*, double *x*)

- double `gsl_sf_laguerre_3`(double *a*, double *x*)

- int `gsl_sf_laguerre_1_e`(double *a*, double *x*, [gsl_sf_result](https://www.gnu.org/software/gsl/doc/html/specfunc.html#c.gsl_sf_result) * *result*)

- int `gsl_sf_laguerre_2_e`(double *a*, double *x*, [gsl_sf_result](https://www.gnu.org/software/gsl/doc/html/specfunc.html#c.gsl_sf_result) * *result*)

- int `gsl_sf_laguerre_3_e`(double *a*, double *x*, [gsl_sf_result](https://www.gnu.org/software/gsl/doc/html/specfunc.html#c.gsl_sf_result) * *result*)

  These routines evaluate the generalized Laguerre polynomials ![L^a_1(x)](https://www.gnu.org/software/gsl/doc/html/_images/math/3c9de28306415c095819d4866665bd09b83f2af7.png), ![L^a_2(x)](https://www.gnu.org/software/gsl/doc/html/_images/math/f646a830b64ca0adc7e786e35853ea41082e53ac.png), ![L^a_3(x)](https://www.gnu.org/software/gsl/doc/html/_images/math/492ac2ec35dfe232df7dc70afec5080da9eb96ef.png) using explicit representations.

- double `gsl_sf_laguerre_n`(const int *n*, const double *a*, const double *x*)

- int `gsl_sf_laguerre_n_e`(int *n*, double *a*, double *x*, [gsl_sf_result](https://www.gnu.org/software/gsl/doc/html/specfunc.html#c.gsl_sf_result) * *result*)

  These routines evaluate the generalized Laguerre polynomials ![L^a_n(x)](https://www.gnu.org/software/gsl/doc/html/_images/math/19e2c845d4653718c7f43ac05a1d2e1aea4d4750.png) for ![a > -1](https://www.gnu.org/software/gsl/doc/html/_images/math/adf9a4a0bef1d278ccb19fc5388f03d9853aea6c.png), ![n \ge 0](https://www.gnu.org/software/gsl/doc/html/_images/math/a4a2e96a00e9bf71dda87e79f4f76f7ad0bca7cf.png).

## Lambert W Functions

Lambert’s W functions, ![W(x)](https://www.gnu.org/software/gsl/doc/html/_images/math/7c30eca5a07d6837f588ebd9eaea953b9a14ac41.png), are defined to be solutions of the equation ![W(x) \exp(W(x)) = x](https://www.gnu.org/software/gsl/doc/html/_images/math/49c299c78b082bead85c647b3e7d20889fa797de.png). This function has multiple branches for ![x < 0](https://www.gnu.org/software/gsl/doc/html/_images/math/93f2a3ce4b90fb07a7dc88c1bc60a6a43c93f45b.png); however, it has only two real-valued branches. We define ![W_0(x)](https://www.gnu.org/software/gsl/doc/html/_images/math/faa96290620004db7b04ca0f2ce15c977368c812.png) to be the principal branch, where ![W > -1](https://www.gnu.org/software/gsl/doc/html/_images/math/29e0d85b337245bc561f6cd61cec719c5c037867.png) for ![x < 0](https://www.gnu.org/software/gsl/doc/html/_images/math/93f2a3ce4b90fb07a7dc88c1bc60a6a43c93f45b.png), and ![W_{-1}(x)](https://www.gnu.org/software/gsl/doc/html/_images/math/38103f208d6a39a424b54f989eac94931d4dcbb1.png) to be the other real branch, where ![W < -1](https://www.gnu.org/software/gsl/doc/html/_images/math/0eb6c71accd3c8495f8b102d2a90fd1a69acb4fc.png) for ![x < 0](https://www.gnu.org/software/gsl/doc/html/_images/math/93f2a3ce4b90fb07a7dc88c1bc60a6a43c93f45b.png). The Lambert functions are declared in the header file `gsl_sf_lambert.h`.

- double `gsl_sf_lambert_W0`(double *x*)

- int `gsl_sf_lambert_W0_e`(double *x*, [gsl_sf_result](https://www.gnu.org/software/gsl/doc/html/specfunc.html#c.gsl_sf_result) * *result*)

  These compute the principal branch of the Lambert W function, ![W_0(x)](https://www.gnu.org/software/gsl/doc/html/_images/math/faa96290620004db7b04ca0f2ce15c977368c812.png).

- double `gsl_sf_lambert_Wm1`(double *x*)

- int `gsl_sf_lambert_Wm1_e`(double *x*, [gsl_sf_result](https://www.gnu.org/software/gsl/doc/html/specfunc.html#c.gsl_sf_result) * *result*)

  These compute the secondary real-valued branch of the Lambert W function, ![W_{-1}(x)](https://www.gnu.org/software/gsl/doc/html/_images/math/38103f208d6a39a424b54f989eac94931d4dcbb1.png).

## Legendre Functions and Spherical Harmonics

The Legendre Functions and Legendre Polynomials are described in Abramowitz & Stegun, Chapter 8. These functions are declared in the header file `gsl_sf_legendre.h`.

### Legendre Polynomials

- double `gsl_sf_legendre_P1`(double *x*)

- double `gsl_sf_legendre_P2`(double *x*)

- double `gsl_sf_legendre_P3`(double *x*)

- int `gsl_sf_legendre_P1_e`(double *x*, [gsl_sf_result](https://www.gnu.org/software/gsl/doc/html/specfunc.html#c.gsl_sf_result) * *result*)

- int `gsl_sf_legendre_P2_e`(double *x*, [gsl_sf_result](https://www.gnu.org/software/gsl/doc/html/specfunc.html#c.gsl_sf_result) * *result*)

- int `gsl_sf_legendre_P3_e`(double *x*, [gsl_sf_result](https://www.gnu.org/software/gsl/doc/html/specfunc.html#c.gsl_sf_result) * *result*)

  These functions evaluate the Legendre polynomials ![P_l(x)](https://www.gnu.org/software/gsl/doc/html/_images/math/7dbb28d5ba6482ac53afe63180d058245cc61fe2.png) using explicit representations for ![l = 1, 2, 3](https://www.gnu.org/software/gsl/doc/html/_images/math/1de1d28b54b4b98b58f8e55ea3a607a4825575c6.png).

- double `gsl_sf_legendre_Pl`(int *l*, double *x*)

- int `gsl_sf_legendre_Pl_e`(int *l*, double *x*, [gsl_sf_result](https://www.gnu.org/software/gsl/doc/html/specfunc.html#c.gsl_sf_result) * *result*)

  These functions evaluate the Legendre polynomial ![P_l(x)](https://www.gnu.org/software/gsl/doc/html/_images/math/7dbb28d5ba6482ac53afe63180d058245cc61fe2.png) for a specific value of `l`, `x` subject to ![l \ge 0](https://www.gnu.org/software/gsl/doc/html/_images/math/911176bca9e219aa2993a1c97d05c5ed2d848ca5.png) and ![|x| \le 1](https://www.gnu.org/software/gsl/doc/html/_images/math/9bc759fdcccfd43a9a4d9673e82e34c393ebb9ba.png).

- int `gsl_sf_legendre_Pl_array`(int *lmax*, double *x*, double *result_array[]*)

- int `gsl_sf_legendre_Pl_deriv_array`(int *lmax*, double *x*, double *result_array[]*, double *result_deriv_array[]*)

  These functions compute arrays of Legendre polynomials ![P_l(x)](https://www.gnu.org/software/gsl/doc/html/_images/math/7dbb28d5ba6482ac53afe63180d058245cc61fe2.png) and derivatives ![dP_l(x)/dx](https://www.gnu.org/software/gsl/doc/html/_images/math/2a4def7e7978a4db8a4097162f2c39603cb18f75.png) for ![l = 0, \dots, lmax](https://www.gnu.org/software/gsl/doc/html/_images/math/34ebda78fd6b692f4556e06867ebb6184875cdf9.png) and ![|x| \le 1](https://www.gnu.org/software/gsl/doc/html/_images/math/9bc759fdcccfd43a9a4d9673e82e34c393ebb9ba.png).

- double `gsl_sf_legendre_Q0`(double *x*)

- int `gsl_sf_legendre_Q0_e`(double *x*, [gsl_sf_result](https://www.gnu.org/software/gsl/doc/html/specfunc.html#c.gsl_sf_result) * *result*)

  These routines compute the Legendre function ![Q_0(x)](https://www.gnu.org/software/gsl/doc/html/_images/math/e4519936094b06e3bd4c152e743bfd3daf5dffca.png) for ![x > -1](https://www.gnu.org/software/gsl/doc/html/_images/math/d16b9945a3b63aefe250aef2de3c86dfa323def5.png) and ![x \ne 1](https://www.gnu.org/software/gsl/doc/html/_images/math/a15a602a7c1b14af061aad3ea81845bcf8b3048e.png).

- double `gsl_sf_legendre_Q1`(double *x*)

- int `gsl_sf_legendre_Q1_e`(double *x*, [gsl_sf_result](https://www.gnu.org/software/gsl/doc/html/specfunc.html#c.gsl_sf_result) * *result*)

  These routines compute the Legendre function ![Q_1(x)](https://www.gnu.org/software/gsl/doc/html/_images/math/cadca37e62915fc86f31080f431c0229ef10429b.png) for ![x > -1](https://www.gnu.org/software/gsl/doc/html/_images/math/d16b9945a3b63aefe250aef2de3c86dfa323def5.png) and ![x \ne 1](https://www.gnu.org/software/gsl/doc/html/_images/math/a15a602a7c1b14af061aad3ea81845bcf8b3048e.png).

- double `gsl_sf_legendre_Ql`(int *l*, double *x*)

- int `gsl_sf_legendre_Ql_e`(int *l*, double *x*, [gsl_sf_result](https://www.gnu.org/software/gsl/doc/html/specfunc.html#c.gsl_sf_result) * *result*)

  These routines compute the Legendre function ![Q_l(x)](https://www.gnu.org/software/gsl/doc/html/_images/math/5cff1133688afcbd9b268d0a1b12d535484a1cf1.png) for ![x > -1](https://www.gnu.org/software/gsl/doc/html/_images/math/d16b9945a3b63aefe250aef2de3c86dfa323def5.png), ![x \ne 1](https://www.gnu.org/software/gsl/doc/html/_images/math/a15a602a7c1b14af061aad3ea81845bcf8b3048e.png) and ![l \ge 0](https://www.gnu.org/software/gsl/doc/html/_images/math/911176bca9e219aa2993a1c97d05c5ed2d848ca5.png).

### Associated Legendre Polynomials and Spherical Harmonics

The following functions compute the associated Legendre polynomials ![P_l^m(x)](https://www.gnu.org/software/gsl/doc/html/_images/math/a3d22c7040038fa537778a60100862dda88eb3f7.png) which are solutions of the differential equation

![(1 - x^2) {d^2 \over dx^2} P_l^m(x) - 2x {d \over dx} P_l^m(x) + \left( l(l+1) - {m^2 \over 1 - x^2} \right) P_l^m(x) = 0](https://www.gnu.org/software/gsl/doc/html/_images/math/26b0e90b1bdc164ce1ad6704f29cbcfe11113b4d.png)

where the degree ![l](https://www.gnu.org/software/gsl/doc/html/_images/math/a7cc908d8b1de27d7caaa2eab1f302b760d3dbf0.png) and order ![m](https://www.gnu.org/software/gsl/doc/html/_images/math/bcf5f439a74f5012a813f0a8e3e583d709b4c51e.png) satisfy ![0 \le l](https://www.gnu.org/software/gsl/doc/html/_images/math/b739d0906d5de9befad54d1cf2180b3e5c7de5cc.png) and ![0 \le m \le l](https://www.gnu.org/software/gsl/doc/html/_images/math/79d543ee5e359767db6c8b71405117092e222d64.png). The functions ![P_l^m(x)](https://www.gnu.org/software/gsl/doc/html/_images/math/a3d22c7040038fa537778a60100862dda88eb3f7.png) grow combinatorially with ![l](https://www.gnu.org/software/gsl/doc/html/_images/math/a7cc908d8b1de27d7caaa2eab1f302b760d3dbf0.png) and can overflow for ![l](https://www.gnu.org/software/gsl/doc/html/_images/math/a7cc908d8b1de27d7caaa2eab1f302b760d3dbf0.png) larger than about 150. Alternatively, one may calculate normalized associated Legendre polynomials. There are a number of different normalization conventions, and these functions can be stably computed up to degree and order 2700. The following normalizations are provided:

- Schmidt semi-normalization

  Schmidt semi-normalized associated Legendre polynomials are often used in the magnetics community and are defined as

  ![S_l^0(x) &= P_l^0(x) \\ S_l^m(x) &= (-1)^m \sqrt{2 {(l-m)! \over (l+m)!}} P_l^m(x), m > 0](https://www.gnu.org/software/gsl/doc/html/_images/math/d9db78d57a4b76c332e582f10f5661471949061a.png)

  The factor of ![(-1)^m](https://www.gnu.org/software/gsl/doc/html/_images/math/b1c7d755d3dafc5731f0faf9cc70bd7ecdecae33.png) is called the Condon-Shortley phase factor and can be excluded if desired by setting the parameter `csphase = 1` in the functions below.

- Spherical Harmonic Normalization

  The associated Legendre polynomials suitable for calculating spherical harmonics are defined as

  ![Y_l^m(x) = (-1)^m \sqrt{{2l + 1 \over 4 \pi} {(l-m)! \over (l+m)!}} P_l^m(x)](https://www.gnu.org/software/gsl/doc/html/_images/math/a8bde5c4565f545deb4a64f554c0e7dca16f7474.png)

  where again the phase factor ![(-1)^m](https://www.gnu.org/software/gsl/doc/html/_images/math/b1c7d755d3dafc5731f0faf9cc70bd7ecdecae33.png) can be included or excluded if desired.

- Full Normalization

  The fully normalized associated Legendre polynomials are defined as

  ![N_l^m(x) = (-1)^m \sqrt{(l + {1 \over 2}) {(l-m)! \over (l+m)!}} P_l^m(x)](https://www.gnu.org/software/gsl/doc/html/_images/math/9f5493831d58d7d3f0df3f2ba45e0a634df282e3.png)

  and have the property

  ![\int_{-1}^1 N_l^m(x)^2 dx = 1](https://www.gnu.org/software/gsl/doc/html/_images/math/9ca08937337169e63929f163f1c4280d1663b3bc.png)

The normalized associated Legendre routines below use a recurrence relation which is stable up to a degree and order of about 2700. Beyond this, the computed functions could suffer from underflow leading to incorrect results. Routines are provided to compute first and second derivatives ![dP_l^m(x)/dx](https://www.gnu.org/software/gsl/doc/html/_images/math/5b409e36a0681652e1e6a5d9d03ac3c5ffc0bd26.png) and ![d^2 P_l^m(x)/dx^2](https://www.gnu.org/software/gsl/doc/html/_images/math/590b15fd005073e67cea8f73ca44dd298a1a3612.png) as well as their alternate versions ![d P_l^m(\cos{\theta})/d\theta](https://www.gnu.org/software/gsl/doc/html/_images/math/40faa60e45c351c1dbaa3232e02970fa859a151b.png) and ![d^2 P_l^m(\cos{\theta})/d\theta^2](https://www.gnu.org/software/gsl/doc/html/_images/math/439de1fc69b3ee9c9cc1d28fff5ee818e624c4bc.png). While there is a simple scaling relationship between the two forms, the derivatives involving ![\theta](https://www.gnu.org/software/gsl/doc/html/_images/math/413c70072fb314747d7585b05b96ad9763469dfe.png) are heavily used in spherical harmonic expansions and so these routines are also provided.

In the functions below, a parameter of type [`gsl_sf_legendre_t`](https://www.gnu.org/software/gsl/doc/html/specfunc.html#c.gsl_sf_legendre_t) specifies the type of normalization to use. The possible values are

- `gsl_sf_legendre_t`

  ValueDescription`GSL_SF_LEGENDRE_NONE`The unnormalized associated Legendre polynomials ![P_l^m(x)](https://www.gnu.org/software/gsl/doc/html/_images/math/a3d22c7040038fa537778a60100862dda88eb3f7.png)`GSL_SF_LEGENDRE_SCHMIDT`The Schmidt semi-normalized associated Legendre polynomials ![S_l^m(x)](https://www.gnu.org/software/gsl/doc/html/_images/math/934e626275629973c925085ce4a508c5211694fe.png)`GSL_SF_LEGENDRE_SPHARM`The spherical harmonic associated Legendre polynomials ![Y_l^m(x)](https://www.gnu.org/software/gsl/doc/html/_images/math/4dba96f4e8e05e0e5de80db28fa454c93eb22f24.png)`GSL_SF_LEGENDRE_FULL`The fully normalized associated Legendre polynomials ![N_l^m(x)](https://www.gnu.org/software/gsl/doc/html/_images/math/b52b4c370a41cbb8fd269d94d33a58002c60c299.png)

- int `gsl_sf_legendre_array`(const [gsl_sf_legendre_t](https://www.gnu.org/software/gsl/doc/html/specfunc.html#c.gsl_sf_legendre_t) *norm*, const size_t *lmax*, const double *x*, double *result_array[]*)

- int `gsl_sf_legendre_array_e`(const [gsl_sf_legendre_t](https://www.gnu.org/software/gsl/doc/html/specfunc.html#c.gsl_sf_legendre_t) *norm*, const size_t *lmax*, const double *x*, const double *csphase*, double *result_array[]*)

  These functions calculate all normalized associated Legendre polynomials for ![0 \le l \le lmax](https://www.gnu.org/software/gsl/doc/html/_images/math/1dacc9b66c5a63bca3a0f895d58cb7a52745d5a2.png)and ![0 \le m \le l](https://www.gnu.org/software/gsl/doc/html/_images/math/79d543ee5e359767db6c8b71405117092e222d64.png) for ![|x| \le 1](https://www.gnu.org/software/gsl/doc/html/_images/math/9bc759fdcccfd43a9a4d9673e82e34c393ebb9ba.png). The `norm` parameter specifies which normalization is used. The normalized ![P_l^m(x)](https://www.gnu.org/software/gsl/doc/html/_images/math/a3d22c7040038fa537778a60100862dda88eb3f7.png) values are stored in `result_array`, whose minimum size can be obtained from calling [`gsl_sf_legendre_array_n()`](https://www.gnu.org/software/gsl/doc/html/specfunc.html#c.gsl_sf_legendre_array_n). The array index of ![P_l^m(x)](https://www.gnu.org/software/gsl/doc/html/_images/math/a3d22c7040038fa537778a60100862dda88eb3f7.png) is obtained from calling`gsl_sf_legendre_array_index(l, m)`. To include or exclude the Condon-Shortley phase factor of ![(-1)^m](https://www.gnu.org/software/gsl/doc/html/_images/math/b1c7d755d3dafc5731f0faf9cc70bd7ecdecae33.png), set the parameter `csphase` to either ![-1](https://www.gnu.org/software/gsl/doc/html/_images/math/e26aecc9561355dd0e6fbf2742d0ff0d6b58c864.png) or ![1](https://www.gnu.org/software/gsl/doc/html/_images/math/585581da152cce359927e7b9eeadef5ca7c3a42d.png) respectively in the `_e` function. This factor is excluded by default.

- int `gsl_sf_legendre_deriv_array`(const [gsl_sf_legendre_t](https://www.gnu.org/software/gsl/doc/html/specfunc.html#c.gsl_sf_legendre_t) *norm*, const size_t *lmax*, const double *x*, double *result_array[]*, double *result_deriv_array[]*)

- int `gsl_sf_legendre_deriv_array_e`(const [gsl_sf_legendre_t](https://www.gnu.org/software/gsl/doc/html/specfunc.html#c.gsl_sf_legendre_t) *norm*, const size_t *lmax*, const double *x*, const double *csphase*, double *result_array[]*, double *result_deriv_array[]*)

  These functions calculate all normalized associated Legendre functions and their first derivatives up to degree `lmax` for ![|x| < 1](https://www.gnu.org/software/gsl/doc/html/_images/math/596e6c7f75eb9f289d30a08bc5e911aa74197ee0.png). The parameter `norm` specifies the normalization used. The normalized ![P_l^m(x)](https://www.gnu.org/software/gsl/doc/html/_images/math/a3d22c7040038fa537778a60100862dda88eb3f7.png) values and their derivatives ![dP_l^m(x)/dx](https://www.gnu.org/software/gsl/doc/html/_images/math/5b409e36a0681652e1e6a5d9d03ac3c5ffc0bd26.png) are stored in `result_array` and `result_deriv_array` respectively. To include or exclude the Condon-Shortley phase factor of ![(-1)^m](https://www.gnu.org/software/gsl/doc/html/_images/math/b1c7d755d3dafc5731f0faf9cc70bd7ecdecae33.png), set the parameter `csphase` to either ![-1](https://www.gnu.org/software/gsl/doc/html/_images/math/e26aecc9561355dd0e6fbf2742d0ff0d6b58c864.png) or ![1](https://www.gnu.org/software/gsl/doc/html/_images/math/585581da152cce359927e7b9eeadef5ca7c3a42d.png) respectively in the `_e`function. This factor is excluded by default.

- int `gsl_sf_legendre_deriv_alt_array`(const [gsl_sf_legendre_t](https://www.gnu.org/software/gsl/doc/html/specfunc.html#c.gsl_sf_legendre_t) *norm*, const size_t *lmax*, const double *x*, double *result_array[]*, double *result_deriv_array[]*)

- int `gsl_sf_legendre_deriv_alt_array_e`(const [gsl_sf_legendre_t](https://www.gnu.org/software/gsl/doc/html/specfunc.html#c.gsl_sf_legendre_t) *norm*, const size_t *lmax*, const double *x*, const double *csphase*, double *result_array[]*, double *result_deriv_array[]*)

  These functions calculate all normalized associated Legendre functions and their (alternate) first derivatives up to degree `lmax` for ![|x| < 1](https://www.gnu.org/software/gsl/doc/html/_images/math/596e6c7f75eb9f289d30a08bc5e911aa74197ee0.png). The normalized ![P_l^m(x)](https://www.gnu.org/software/gsl/doc/html/_images/math/a3d22c7040038fa537778a60100862dda88eb3f7.png) values and their derivatives ![dP_l^m(\cos{\theta})/d\theta](https://www.gnu.org/software/gsl/doc/html/_images/math/e7913baa4f362fe067305693b52b35a2219f1927.png) are stored in `result_array` and `result_deriv_array` respectively. To include or exclude the Condon-Shortley phase factor of ![(-1)^m](https://www.gnu.org/software/gsl/doc/html/_images/math/b1c7d755d3dafc5731f0faf9cc70bd7ecdecae33.png), set the parameter `csphase` to either ![-1](https://www.gnu.org/software/gsl/doc/html/_images/math/e26aecc9561355dd0e6fbf2742d0ff0d6b58c864.png)or ![1](https://www.gnu.org/software/gsl/doc/html/_images/math/585581da152cce359927e7b9eeadef5ca7c3a42d.png) respectively in the `_e` function. This factor is excluded by default.

- int `gsl_sf_legendre_deriv2_array`(const [gsl_sf_legendre_t](https://www.gnu.org/software/gsl/doc/html/specfunc.html#c.gsl_sf_legendre_t) *norm*, const size_t *lmax*, const double *x*, double *result_array[]*, double *result_deriv_array[]*, double *result_deriv2_array[]*)

- int `gsl_sf_legendre_deriv2_array_e`(const [gsl_sf_legendre_t](https://www.gnu.org/software/gsl/doc/html/specfunc.html#c.gsl_sf_legendre_t) *norm*, const size_t *lmax*, const double *x*, const double *csphase*, double *result_array[]*, double *result_deriv_array[]*, double *result_deriv2_array[]*)

  These functions calculate all normalized associated Legendre functions and their first and second derivatives up to degree `lmax` for ![|x| < 1](https://www.gnu.org/software/gsl/doc/html/_images/math/596e6c7f75eb9f289d30a08bc5e911aa74197ee0.png). The parameter `norm` specifies the normalization used. The normalized ![P_l^m(x)](https://www.gnu.org/software/gsl/doc/html/_images/math/a3d22c7040038fa537778a60100862dda88eb3f7.png), their first derivatives ![dP_l^m(x)/dx](https://www.gnu.org/software/gsl/doc/html/_images/math/5b409e36a0681652e1e6a5d9d03ac3c5ffc0bd26.png), and their second derivatives ![d^2 P_l^m(x)/dx^2](https://www.gnu.org/software/gsl/doc/html/_images/math/590b15fd005073e67cea8f73ca44dd298a1a3612.png) are stored in `result_array`, `result_deriv_array`, and `result_deriv2_array` respectively. To include or exclude the Condon-Shortley phase factor of ![(-1)^m](https://www.gnu.org/software/gsl/doc/html/_images/math/b1c7d755d3dafc5731f0faf9cc70bd7ecdecae33.png), set the parameter `csphase` to either ![-1](https://www.gnu.org/software/gsl/doc/html/_images/math/e26aecc9561355dd0e6fbf2742d0ff0d6b58c864.png) or ![1](https://www.gnu.org/software/gsl/doc/html/_images/math/585581da152cce359927e7b9eeadef5ca7c3a42d.png) respectively in the `_e` function. This factor is excluded by default.

- int `gsl_sf_legendre_deriv2_alt_array`(const [gsl_sf_legendre_t](https://www.gnu.org/software/gsl/doc/html/specfunc.html#c.gsl_sf_legendre_t) *norm*, const size_t *lmax*, const double *x*, double *result_array[]*, double *result_deriv_array[]*, double *result_deriv2_array[]*)

- int `gsl_sf_legendre_deriv2_alt_array_e`(const [gsl_sf_legendre_t](https://www.gnu.org/software/gsl/doc/html/specfunc.html#c.gsl_sf_legendre_t) *norm*, const size_t *lmax*, const double *x*, const double *csphase*, double *result_array[]*, double *result_deriv_array[]*, double *result_deriv2_array[]*)

  These functions calculate all normalized associated Legendre functions and their (alternate) first and second derivatives up to degree `lmax` for ![|x| < 1](https://www.gnu.org/software/gsl/doc/html/_images/math/596e6c7f75eb9f289d30a08bc5e911aa74197ee0.png). The parameter `norm` specifies the normalization used. The normalized ![P_l^m(x)](https://www.gnu.org/software/gsl/doc/html/_images/math/a3d22c7040038fa537778a60100862dda88eb3f7.png), their first derivatives ![dP_l^m(\cos{\theta})/d\theta](https://www.gnu.org/software/gsl/doc/html/_images/math/e7913baa4f362fe067305693b52b35a2219f1927.png), and their second derivatives ![d^2 P_l^m(\cos{\theta})/d\theta^2](https://www.gnu.org/software/gsl/doc/html/_images/math/439de1fc69b3ee9c9cc1d28fff5ee818e624c4bc.png) are stored in `result_array`, `result_deriv_array`, and `result_deriv2_array` respectively. To include or exclude the Condon-Shortley phase factor of ![(-1)^m](https://www.gnu.org/software/gsl/doc/html/_images/math/b1c7d755d3dafc5731f0faf9cc70bd7ecdecae33.png), set the parameter `csphase` to either ![-1](https://www.gnu.org/software/gsl/doc/html/_images/math/e26aecc9561355dd0e6fbf2742d0ff0d6b58c864.png) or ![1](https://www.gnu.org/software/gsl/doc/html/_images/math/585581da152cce359927e7b9eeadef5ca7c3a42d.png) respectively in the `_e` function. This factor is excluded by default.

- size_t `gsl_sf_legendre_array_n`(const size_t *lmax*)

  This function returns the minimum array size for maximum degree `lmax` needed for the array versions of the associated Legendre functions. Size is calculated as the total number of ![P_l^m(x)](https://www.gnu.org/software/gsl/doc/html/_images/math/a3d22c7040038fa537778a60100862dda88eb3f7.png)functions, plus extra space for precomputing multiplicative factors used in the recurrence relations.

- size_t `gsl_sf_legendre_array_index`(const size_t *l*, const size_t *m*)

  This function returns the index into `result_array`, `result_deriv_array`, or `result_deriv2_array` corresponding to ![P_l^m(x)](https://www.gnu.org/software/gsl/doc/html/_images/math/a3d22c7040038fa537778a60100862dda88eb3f7.png), ![P_l^{'m}(x)](https://www.gnu.org/software/gsl/doc/html/_images/math/4837960bab4ed5fd403a7f23c9156c29fbf7d855.png), or ![P_l^{''m}(x)](https://www.gnu.org/software/gsl/doc/html/_images/math/40cda1181fd5a4c27a7b545a068ab2d284ae7336.png). The index is given by ![l(l+1)/2 + m](https://www.gnu.org/software/gsl/doc/html/_images/math/f071cdbd95b3d6d9c9b55df6d4cb143e81ceac47.png).

- double `gsl_sf_legendre_Plm`(int *l*, int *m*, double *x*)

- int `gsl_sf_legendre_Plm_e`(int *l*, int *m*, double *x*, [gsl_sf_result](https://www.gnu.org/software/gsl/doc/html/specfunc.html#c.gsl_sf_result) * *result*)

  These routines compute the associated Legendre polynomial ![P_l^m(x)](https://www.gnu.org/software/gsl/doc/html/_images/math/a3d22c7040038fa537778a60100862dda88eb3f7.png) for ![m \ge 0](https://www.gnu.org/software/gsl/doc/html/_images/math/1d8730c81b20b9e10b32c1cd86a12449f4d02db4.png), ![l \ge m](https://www.gnu.org/software/gsl/doc/html/_images/math/8563adb18e020a4ef24ce6a86840e457a12c8597.png), and ![|x| \le 1](https://www.gnu.org/software/gsl/doc/html/_images/math/9bc759fdcccfd43a9a4d9673e82e34c393ebb9ba.png).

- double `gsl_sf_legendre_sphPlm`(int *l*, int *m*, double *x*)

- int `gsl_sf_legendre_sphPlm_e`(int *l*, int *m*, double *x*, [gsl_sf_result](https://www.gnu.org/software/gsl/doc/html/specfunc.html#c.gsl_sf_result) * *result*)

  These routines compute the normalized associated Legendre polynomial ![\sqrt{(2l+1)/(4\pi)} \sqrt{(l-m)!/(l+m)!} P_l^m(x)](https://www.gnu.org/software/gsl/doc/html/_images/math/fa400862e73e609b82b8d8873faaed2f3047480e.png) suitable for use in spherical harmonics. The parameters must satisfy ![m \ge 0](https://www.gnu.org/software/gsl/doc/html/_images/math/1d8730c81b20b9e10b32c1cd86a12449f4d02db4.png), ![l \ge m](https://www.gnu.org/software/gsl/doc/html/_images/math/8563adb18e020a4ef24ce6a86840e457a12c8597.png), and ![|x| \le 1](https://www.gnu.org/software/gsl/doc/html/_images/math/9bc759fdcccfd43a9a4d9673e82e34c393ebb9ba.png). These routines avoid the overflows that occur for the standard normalization of ![P_l^m(x)](https://www.gnu.org/software/gsl/doc/html/_images/math/a3d22c7040038fa537778a60100862dda88eb3f7.png).

- int `gsl_sf_legendre_Plm_array`(int *lmax*, int *m*, double *x*, double *result_array[]*)

- int `gsl_sf_legendre_Plm_deriv_array`(int *lmax*, int *m*, double *x*, double *result_array[]*, double *result_deriv_array[]*)

  These functions are now deprecated and will be removed in a future release; see [`gsl_sf_legendre_array()`](https://www.gnu.org/software/gsl/doc/html/specfunc.html#c.gsl_sf_legendre_array) and [`gsl_sf_legendre_deriv_array()`](https://www.gnu.org/software/gsl/doc/html/specfunc.html#c.gsl_sf_legendre_deriv_array).

- int `gsl_sf_legendre_sphPlm_array`(int *lmax*, int *m*, double *x*, double *result_array[]*)

- int `gsl_sf_legendre_sphPlm_deriv_array`(int *lmax*, int *m*, double *x*, double *result_array[]*, double *result_deriv_array[]*)

  These functions are now deprecated and will be removed in a future release; see [`gsl_sf_legendre_array()`](https://www.gnu.org/software/gsl/doc/html/specfunc.html#c.gsl_sf_legendre_array) and [`gsl_sf_legendre_deriv_array()`](https://www.gnu.org/software/gsl/doc/html/specfunc.html#c.gsl_sf_legendre_deriv_array).

- int `gsl_sf_legendre_array_size`(const int *lmax*, const int *m*)

  This function is now deprecated and will be removed in a future release.

### Conical Functions

The Conical Functions ![P^\mu_{-(1/2)+i\lambda}(x)](https://www.gnu.org/software/gsl/doc/html/_images/math/92a103398830717987210d238b79a84ca9b84319.png) and ![Q^\mu_{-(1/2)+i\lambda}](https://www.gnu.org/software/gsl/doc/html/_images/math/018360ec84729d87c5a13a107222f573b1869f00.png) are described in Abramowitz & Stegun, Section 8.12.

- double `gsl_sf_conicalP_half`(double *lambda*, double *x*)

- int `gsl_sf_conicalP_half_e`(double *lambda*, double *x*, [gsl_sf_result](https://www.gnu.org/software/gsl/doc/html/specfunc.html#c.gsl_sf_result) * *result*)

  These routines compute the irregular Spherical Conical Function ![P^{1/2}_{-1/2 + i \lambda}(x)](https://www.gnu.org/software/gsl/doc/html/_images/math/375e9fd8ff54ff8fe2efc061010f3bfede7e388c.png) for ![x > -1](https://www.gnu.org/software/gsl/doc/html/_images/math/d16b9945a3b63aefe250aef2de3c86dfa323def5.png).

- double `gsl_sf_conicalP_mhalf`(double *lambda*, double *x*)

- int `gsl_sf_conicalP_mhalf_e`(double *lambda*, double *x*, [gsl_sf_result](https://www.gnu.org/software/gsl/doc/html/specfunc.html#c.gsl_sf_result) * *result*)

  These routines compute the regular Spherical Conical Function ![P^{-1/2}_{-1/2 + i \lambda}(x)](https://www.gnu.org/software/gsl/doc/html/_images/math/c32523cdbc2725394bdb29f1f5a2994bc42bdf39.png) for ![x > -1](https://www.gnu.org/software/gsl/doc/html/_images/math/d16b9945a3b63aefe250aef2de3c86dfa323def5.png).

- double `gsl_sf_conicalP_0`(double *lambda*, double *x*)

- int `gsl_sf_conicalP_0_e`(double *lambda*, double *x*, [gsl_sf_result](https://www.gnu.org/software/gsl/doc/html/specfunc.html#c.gsl_sf_result) * *result*)

  These routines compute the conical function ![P^0_{-1/2 + i \lambda}(x)](https://www.gnu.org/software/gsl/doc/html/_images/math/e443c41838f78e4b9b2aab1e147382024b364c9b.png) for ![x > -1](https://www.gnu.org/software/gsl/doc/html/_images/math/d16b9945a3b63aefe250aef2de3c86dfa323def5.png).

- double `gsl_sf_conicalP_1`(double *lambda*, double *x*)

- int `gsl_sf_conicalP_1_e`(double *lambda*, double *x*, [gsl_sf_result](https://www.gnu.org/software/gsl/doc/html/specfunc.html#c.gsl_sf_result) * *result*)

  These routines compute the conical function ![P^1_{-1/2 + i \lambda}(x)](https://www.gnu.org/software/gsl/doc/html/_images/math/719a3ea03844f3b08a12f3d3bf3e6aad8c6ab1bc.png) for ![x > -1](https://www.gnu.org/software/gsl/doc/html/_images/math/d16b9945a3b63aefe250aef2de3c86dfa323def5.png).

- double `gsl_sf_conicalP_sph_reg`(int *l*, double *lambda*, double *x*)

- int `gsl_sf_conicalP_sph_reg_e`(int *l*, double *lambda*, double *x*, [gsl_sf_result](https://www.gnu.org/software/gsl/doc/html/specfunc.html#c.gsl_sf_result) * *result*)

  These routines compute the Regular Spherical Conical Function ![P^{-1/2-l}_{-1/2 + i \lambda}(x)](https://www.gnu.org/software/gsl/doc/html/_images/math/b460136a738605a0dd0cbea44dac35c05f5c0bb0.png) for ![x > -1](https://www.gnu.org/software/gsl/doc/html/_images/math/d16b9945a3b63aefe250aef2de3c86dfa323def5.png) and ![l \ge -1](https://www.gnu.org/software/gsl/doc/html/_images/math/391d51ba08ea93cac708fad80403579821ef4190.png).

- double `gsl_sf_conicalP_cyl_reg`(int *m*, double *lambda*, double *x*)

- int `gsl_sf_conicalP_cyl_reg_e`(int *m*, double *lambda*, double *x*, [gsl_sf_result](https://www.gnu.org/software/gsl/doc/html/specfunc.html#c.gsl_sf_result) * *result*)

  These routines compute the Regular Cylindrical Conical Function ![P^{-m}_{-1/2 + i \lambda}(x)](https://www.gnu.org/software/gsl/doc/html/_images/math/8761b06046cae0bfd0cbf4e71ed6eaa5500344b3.png) for ![x > -1](https://www.gnu.org/software/gsl/doc/html/_images/math/d16b9945a3b63aefe250aef2de3c86dfa323def5.png) and ![m \ge -1](https://www.gnu.org/software/gsl/doc/html/_images/math/430e65431e132a1ab444fc7e226d27116a3f4d09.png).

### Radial Functions for Hyperbolic Space

The following spherical functions are specializations of Legendre functions which give the regular eigenfunctions of the Laplacian on a 3-dimensional hyperbolic space ![H^3](https://www.gnu.org/software/gsl/doc/html/_images/math/ac1282bd60f88fc1875628172c30bd28b0adf1cc.png). Of particular interest is the flat limit, ![\lambda \to \infty](https://www.gnu.org/software/gsl/doc/html/_images/math/7e78ec6a81407951fcaec214c8ab435d5b3dcf0f.png), ![\eta \to 0](https://www.gnu.org/software/gsl/doc/html/_images/math/b9f7f2f213e2788f0e51a5c8ecf7d7243c60c1d2.png), ![\lambda\eta](https://www.gnu.org/software/gsl/doc/html/_images/math/f5aeab2a1d3502428b8cc09f8f71a6bd95b1de40.png) fixed.

- double `gsl_sf_legendre_H3d_0`(double *lambda*, double *eta*)

- int `gsl_sf_legendre_H3d_0_e`(double *lambda*, double *eta*, [gsl_sf_result](https://www.gnu.org/software/gsl/doc/html/specfunc.html#c.gsl_sf_result) * *result*)

  These routines compute the zeroth radial eigenfunction of the Laplacian on the 3-dimensional hyperbolic space,![L^{H3d}_0(\lambda,\eta) := {\sin(\lambda\eta) \over \lambda\sinh(\eta)}](https://www.gnu.org/software/gsl/doc/html/_images/math/c09087183c03126353c49bbbe7e75a8a00c382dd.png)for ![\eta \ge 0](https://www.gnu.org/software/gsl/doc/html/_images/math/58384bf02bfbf61d8c7c5d7481bff5b52ffa8f62.png). In the flat limit this takes the form ![L^{H3d}_0(\lambda,\eta) = j_0(\lambda\eta)](https://www.gnu.org/software/gsl/doc/html/_images/math/425507b74e68b3d5a9541f653ca7673e76daa9d5.png).

  这些程序计算三维超几何空间的拉普拉斯算子的零次径向本征函数，$$L^{H3d}_0(\lambda,\eta) := {\sin(\lambda\eta) \over \lambda\sinh(\eta)}$$for $$\eta \ge 0$$. 在平滑极限下其形式为 $$L^{H3d}_0(\lambda,\eta) = j_0(\lambda\eta)$$.

- double `gsl_sf_legendre_H3d_1`(double *lambda*, double *eta*)

- int `gsl_sf_legendre_H3d_1_e`(double *lambda*, double *eta*, [gsl_sf_result](https://www.gnu.org/software/gsl/doc/html/specfunc.html#c.gsl_sf_result) * *result*)

  These routines compute the first radial eigenfunction of the Laplacian on the 3-dimensional hyperbolic space,![L^{H3d}_1(\lambda,\eta) := {1\over\sqrt{\lambda^2 + 1}} {\left(\sin(\lambda \eta)\over \lambda \sinh(\eta)\right)} \left(\coth(\eta) - \lambda \cot(\lambda\eta)\right)](https://www.gnu.org/software/gsl/doc/html/_images/math/5097c3db905c0261de082f29e5abf390d2a43e6f.png)for ![\eta \ge 0](https://www.gnu.org/software/gsl/doc/html/_images/math/58384bf02bfbf61d8c7c5d7481bff5b52ffa8f62.png) In the flat limit this takes the form ![L^{H3d}_1(\lambda,\eta) = j_1(\lambda\eta)](https://www.gnu.org/software/gsl/doc/html/_images/math/a9a7782002c9775ec109f9d21c384bec21bf141b.png).

- double `gsl_sf_legendre_H3d`(int *l*, double *lambda*, double *eta*)

- int `gsl_sf_legendre_H3d_e`(int *l*, double *lambda*, double *eta*, [gsl_sf_result](https://www.gnu.org/software/gsl/doc/html/specfunc.html#c.gsl_sf_result) * *result*)

  These routines compute the `l`-th radial eigenfunction of the Laplacian on the 3-dimensional hyperbolic space ![\eta \ge 0](https://www.gnu.org/software/gsl/doc/html/_images/math/58384bf02bfbf61d8c7c5d7481bff5b52ffa8f62.png) and ![l \ge 0](https://www.gnu.org/software/gsl/doc/html/_images/math/911176bca9e219aa2993a1c97d05c5ed2d848ca5.png). In the flat limit this takes the form ![L^{H3d}_l(\lambda,\eta) = j_l(\lambda\eta)](https://www.gnu.org/software/gsl/doc/html/_images/math/ef33f5052efc450c10cad91b38359046f357eb3e.png).

- int `gsl_sf_legendre_H3d_array`(int *lmax*, double *lambda*, double *eta*, double *result_array[]*)

  This function computes an array of radial eigenfunctions ![L^{H3d}_l( \lambda, \eta)](https://www.gnu.org/software/gsl/doc/html/_images/math/582a8d2e9e486de46e92d979b6db277f4f3410de.png) for ![0 \le l \le lmax](https://www.gnu.org/software/gsl/doc/html/_images/math/1dacc9b66c5a63bca3a0f895d58cb7a52745d5a2.png).

## Logarithm and Related Functions

Information on the properties of the Logarithm function can be found in Abramowitz & Stegun, Chapter 4. The functions described in this section are declared in the header file `gsl_sf_log.h`.

- double `gsl_sf_log`(double *x*)

- int `gsl_sf_log_e`(double *x*, [gsl_sf_result](https://www.gnu.org/software/gsl/doc/html/specfunc.html#c.gsl_sf_result) * *result*)

  These routines compute the logarithm of `x`, ![\log(x)](https://www.gnu.org/software/gsl/doc/html/_images/math/027153da3d62de1e09d2fa15790d2775800fcf76.png), for ![x > 0](https://www.gnu.org/software/gsl/doc/html/_images/math/fe20eedd0aef40a657ad04835ff74fa8a73efb56.png).

- double `gsl_sf_log_abs`(double *x*)

- int `gsl_sf_log_abs_e`(double *x*, [gsl_sf_result](https://www.gnu.org/software/gsl/doc/html/specfunc.html#c.gsl_sf_result) * *result*)

  These routines compute the logarithm of the magnitude of `x`, ![\log(|x|)](https://www.gnu.org/software/gsl/doc/html/_images/math/679126d0b212a63e190cd4247e8c15a7b095c8f7.png), for ![x \ne 0](https://www.gnu.org/software/gsl/doc/html/_images/math/41cbe33fb96924eb79bc5c017b30f8dc9cbbb606.png).

- int `gsl_sf_complex_log_e`(double *zr*, double *zi*, [gsl_sf_result](https://www.gnu.org/software/gsl/doc/html/specfunc.html#c.gsl_sf_result) * *lnr*, [gsl_sf_result](https://www.gnu.org/software/gsl/doc/html/specfunc.html#c.gsl_sf_result) * *theta*)

  This routine computes the complex logarithm of ![z = z_r + i z_i](https://www.gnu.org/software/gsl/doc/html/_images/math/2bd1785b0de8324fc35b8e3a6810476232e4ccf1.png). The results are returned as `lnr`, `theta` such that ![\exp(lnr + i \theta) = z_r + i z_i](https://www.gnu.org/software/gsl/doc/html/_images/math/179540bfba0bde22198e941a453be8c3d5e09e7e.png), where ![\theta](https://www.gnu.org/software/gsl/doc/html/_images/math/413c70072fb314747d7585b05b96ad9763469dfe.png) lies in the range ![[-\pi,\pi]](https://www.gnu.org/software/gsl/doc/html/_images/math/c4dc04287ff1eead37e48f3c3e590053640c817f.png).

- double `gsl_sf_log_1plusx`(double *x*)

- int `gsl_sf_log_1plusx_e`(double *x*, [gsl_sf_result](https://www.gnu.org/software/gsl/doc/html/specfunc.html#c.gsl_sf_result) * *result*)

  These routines compute ![\log(1 + x)](https://www.gnu.org/software/gsl/doc/html/_images/math/1f79d47828946d3b89c402137afa443329647bef.png) for ![x > -1](https://www.gnu.org/software/gsl/doc/html/_images/math/d16b9945a3b63aefe250aef2de3c86dfa323def5.png) using an algorithm that is accurate for small `x`.

- double `gsl_sf_log_1plusx_mx`(double *x*)

- int `gsl_sf_log_1plusx_mx_e`(double *x*, [gsl_sf_result](https://www.gnu.org/software/gsl/doc/html/specfunc.html#c.gsl_sf_result) * *result*)

  These routines compute ![\log(1 + x) - x](https://www.gnu.org/software/gsl/doc/html/_images/math/73b5783145dc93300223e068f49affc9157cef83.png) for ![x > -1](https://www.gnu.org/software/gsl/doc/html/_images/math/d16b9945a3b63aefe250aef2de3c86dfa323def5.png) using an algorithm that is accurate for small `x`.

## Mathieu Functions

The routines described in this section compute the angular and radial Mathieu functions, and their characteristic values. Mathieu functions are the solutions of the following two differential equations:

![{{d^2 y}\over{d v^2}}& + (a - 2q\cos 2v)y  = 0 \\ {{d^2 f}\over{d u^2}}& - (a - 2q\cosh 2u)f  = 0](https://www.gnu.org/software/gsl/doc/html/_images/math/15c535459c089d69cdc1fbe939d8b1897c69079c.png)

The angular Mathieu functions ![ce_r(x,q)](https://www.gnu.org/software/gsl/doc/html/_images/math/3bbfb38bc8ca0fa96b2cb8ca596b9cc7eedfe18f.png), ![se_r(x,q)](https://www.gnu.org/software/gsl/doc/html/_images/math/cd428141c1d5760e59c9af64fb23fd685ad0f7f8.png) are the even and odd periodic solutions of the first equation, which is known as Mathieu’s equation. These exist only for the discrete sequence of characteristic values ![a = a_r(q)](https://www.gnu.org/software/gsl/doc/html/_images/math/f6229de96225a27f2d0944c4a0d9f2e344d40447.png) (even-periodic) and ![a = b_r(q)](https://www.gnu.org/software/gsl/doc/html/_images/math/ca2cc2640770ca52217b2d9197153e53e4cd9612.png) (odd-periodic).

The radial Mathieu functions ![Mc^{(j)}_{r}(z,q)](https://www.gnu.org/software/gsl/doc/html/_images/math/1862ff7fc183cc2b47199e73408833a95df9c503.png) and ![Ms^{(j)}_{r}(z,q)](https://www.gnu.org/software/gsl/doc/html/_images/math/d5f1241623ec259cdcbb691e23ffec27cc86994c.png) are the solutions of the second equation, which is referred to as Mathieu’s modified equation. The radial Mathieu functions of the first, second, third and fourth kind are denoted by the parameter ![j](https://www.gnu.org/software/gsl/doc/html/_images/math/7bc2ad053edeecad39b867d60495ac2caa4555dc.png), which takes the value 1, 2, 3 or 4.

For more information on the Mathieu functions, see Abramowitz and Stegun, Chapter 20. These functions are defined in the header file `gsl_sf_mathieu.h`.

### Mathieu Function Workspace

The Mathieu functions can be computed for a single order or for multiple orders, using array-based routines. The array-based routines require a preallocated workspace.

- `gsl_sf_mathieu_workspace`

  Workspace required for array-based routines

- [gsl_sf_mathieu_workspace](https://www.gnu.org/software/gsl/doc/html/specfunc.html#c.gsl_sf_mathieu_workspace) * `gsl_sf_mathieu_alloc`(size_t *n*, double *qmax*)

  This function returns a workspace for the array versions of the Mathieu routines. The arguments n and `qmax` specify the maximum order and ![q](https://www.gnu.org/software/gsl/doc/html/_images/math/eb51bb30e72e99a9fdffeb315d1aaa91a5859951.png)-value of Mathieu functions which can be computed with this workspace.

- void `gsl_sf_mathieu_free`([gsl_sf_mathieu_workspace](https://www.gnu.org/software/gsl/doc/html/specfunc.html#c.gsl_sf_mathieu_workspace) * *work*)

  This function frees the workspace `work`.

### Mathieu Function Characteristic Values



- int `gsl_sf_mathieu_a`(int *n*, double *q*)

- int `gsl_sf_mathieu_a_e`(int *n*, double *q*, [gsl_sf_result](https://www.gnu.org/software/gsl/doc/html/specfunc.html#c.gsl_sf_result) * *result*)

- int `gsl_sf_mathieu_b`(int *n*, double *q*)

- int `gsl_sf_mathieu_b_e`(int *n*, double *q*, [gsl_sf_result](https://www.gnu.org/software/gsl/doc/html/specfunc.html#c.gsl_sf_result) * *result*)

  These routines compute the characteristic values ![a_n(q)](https://www.gnu.org/software/gsl/doc/html/_images/math/943f1f5e3dab1569a60df9d8a5a44fc52efe76b2.png), ![b_n(q)](https://www.gnu.org/software/gsl/doc/html/_images/math/2a0d369296c6df16f34dd769e6e80e71ca19fcb1.png) of the Mathieu functions ![ce_n(q,x)](https://www.gnu.org/software/gsl/doc/html/_images/math/1eaed6961d5d1aaf493c31cbb04a9904bb6b32ac.png) and ![se_n(q,x)](https://www.gnu.org/software/gsl/doc/html/_images/math/6c3b45f5a7fcbb81e57d3e34b891119e5b5aea5d.png), respectively.

- int `gsl_sf_mathieu_a_array`(int *order_min*, int *order_max*, double *q*, [gsl_sf_mathieu_workspace](https://www.gnu.org/software/gsl/doc/html/specfunc.html#c.gsl_sf_mathieu_workspace)* *work*, double *result_array[]*)

- int `gsl_sf_mathieu_b_array`(int *order_min*, int *order_max*, double *q*, [gsl_sf_mathieu_workspace](https://www.gnu.org/software/gsl/doc/html/specfunc.html#c.gsl_sf_mathieu_workspace)* *work*, double *result_array[]*)

  These routines compute a series of Mathieu characteristic values ![a_n(q)](https://www.gnu.org/software/gsl/doc/html/_images/math/943f1f5e3dab1569a60df9d8a5a44fc52efe76b2.png), ![b_n(q)](https://www.gnu.org/software/gsl/doc/html/_images/math/2a0d369296c6df16f34dd769e6e80e71ca19fcb1.png) for ![n](https://www.gnu.org/software/gsl/doc/html/_images/math/a24554f1502cf204e7aa24a6c064962c3504de48.png) from `order_min` to `order_max` inclusive, storing the results in the array `result_array`.

### Angular Mathieu Functions



- int `gsl_sf_mathieu_ce`(int *n*, double *q*, double *x*)

- int `gsl_sf_mathieu_ce_e`(int *n*, double *q*, double *x*, [gsl_sf_result](https://www.gnu.org/software/gsl/doc/html/specfunc.html#c.gsl_sf_result) * *result*)

- int `gsl_sf_mathieu_se`(int *n*, double *q*, double *x*)

- int `gsl_sf_mathieu_se_e`(int *n*, double *q*, double *x*, [gsl_sf_result](https://www.gnu.org/software/gsl/doc/html/specfunc.html#c.gsl_sf_result) * *result*)

  These routines compute the angular Mathieu functions ![ce_n(q,x)](https://www.gnu.org/software/gsl/doc/html/_images/math/1eaed6961d5d1aaf493c31cbb04a9904bb6b32ac.png) and ![se_n(q,x)](https://www.gnu.org/software/gsl/doc/html/_images/math/6c3b45f5a7fcbb81e57d3e34b891119e5b5aea5d.png), respectively.

- int `gsl_sf_mathieu_ce_array`(int *nmin*, int *nmax*, double *q*, double *x*, [gsl_sf_mathieu_workspace](https://www.gnu.org/software/gsl/doc/html/specfunc.html#c.gsl_sf_mathieu_workspace)* *work*, double *result_array[]*)

- int `gsl_sf_mathieu_se_array`(int *nmin*, int *nmax*, double *q*, double *x*, [gsl_sf_mathieu_workspace](https://www.gnu.org/software/gsl/doc/html/specfunc.html#c.gsl_sf_mathieu_workspace)* *work*, double *result_array[]*)

  These routines compute a series of the angular Mathieu functions ![ce_n(q,x)](https://www.gnu.org/software/gsl/doc/html/_images/math/1eaed6961d5d1aaf493c31cbb04a9904bb6b32ac.png) and ![se_n(q,x)](https://www.gnu.org/software/gsl/doc/html/_images/math/6c3b45f5a7fcbb81e57d3e34b891119e5b5aea5d.png) of order ![n](https://www.gnu.org/software/gsl/doc/html/_images/math/a24554f1502cf204e7aa24a6c064962c3504de48.png) from `nmin` to `nmax` inclusive, storing the results in the array `result_array`.

### Radial Mathieu Functions



- int `gsl_sf_mathieu_Mc`(int *j*, int *n*, double *q*, double *x*)

- int `gsl_sf_mathieu_Mc_e`(int *j*, int *n*, double *q*, double *x*, [gsl_sf_result](https://www.gnu.org/software/gsl/doc/html/specfunc.html#c.gsl_sf_result) * *result*)

- int `gsl_sf_mathieu_Ms`(int *j*, int *n*, double *q*, double *x*)

- int `gsl_sf_mathieu_Ms_e`(int *j*, int *n*, double *q*, double *x*, [gsl_sf_result](https://www.gnu.org/software/gsl/doc/html/specfunc.html#c.gsl_sf_result) * *result*)

  These routines compute the radial `j`-th kind Mathieu functions ![Mc_n^{(j)}(q,x)](https://www.gnu.org/software/gsl/doc/html/_images/math/a4eb0a0b24b84105b8b24174da4793b64a86d943.png) and ![Ms_n^{(j)}(q,x)](https://www.gnu.org/software/gsl/doc/html/_images/math/ccdbcfd397ef66eeb29f8f2685c94a925e451dc6.png)of order `n`.The allowed values of `j` are 1 and 2. The functions for ![j = 3,4](https://www.gnu.org/software/gsl/doc/html/_images/math/5514f50a23bcf94d778ef9ee64cb09e9d6a68a70.png) can be computed as ![M_n^{(3)} = M_n^{(1)} + iM_n^{(2)}](https://www.gnu.org/software/gsl/doc/html/_images/math/111595e79b3d7f23680e9cb59a327a1e78303f4e.png) and ![M_n^{(4)} = M_n^{(1)} - iM_n^{(2)}](https://www.gnu.org/software/gsl/doc/html/_images/math/4aa43c0847277bdbfee6e8150fdb2deb03de1d06.png), where ![M_n^{(j)} = Mc_n^{(j)}](https://www.gnu.org/software/gsl/doc/html/_images/math/4c23b98123644896c75db65d2158e12ba92b1d54.png) or ![Ms_n^{(j)}](https://www.gnu.org/software/gsl/doc/html/_images/math/b23e09bba32cd22d6718078945aa3942972a3300.png).

- int `gsl_sf_mathieu_Mc_array`(int *j*, int *nmin*, int *nmax*, double *q*, double *x*, [gsl_sf_mathieu_workspace](https://www.gnu.org/software/gsl/doc/html/specfunc.html#c.gsl_sf_mathieu_workspace)* *work*, double *result_array[]*)

- int `gsl_sf_mathieu_Ms_array`(int *j*, int *nmin*, int *nmax*, double *q*, double *x*, [gsl_sf_mathieu_workspace](https://www.gnu.org/software/gsl/doc/html/specfunc.html#c.gsl_sf_mathieu_workspace)* *work*, double *result_array[]*)

  These routines compute a series of the radial Mathieu functions of kind `j`, with order from `nmin` to `nmax` inclusive, storing the results in the array `result_array`.

## Power Function

The following functions are equivalent to the function [`gsl_pow_int()`](https://www.gnu.org/software/gsl/doc/html/math.html#c.gsl_pow_int) with an error estimate. These functions are declared in the header file `gsl_sf_pow_int.h`.

- double `gsl_sf_pow_int`(double *x*, int *n*)

- int `gsl_sf_pow_int_e`(double *x*, int *n*, [gsl_sf_result](https://www.gnu.org/software/gsl/doc/html/specfunc.html#c.gsl_sf_result) * *result*)

  These routines compute the power ![x^n](https://www.gnu.org/software/gsl/doc/html/_images/math/ca5dc9582574e5893064ffb5d9f66578c649ffb7.png) for integer `n`. The power is computed using the minimum number of multiplications. For example, ![x^8](https://www.gnu.org/software/gsl/doc/html/_images/math/b77674817382f8d63ad5d32d451ab476be4c354e.png) is computed as ![((x^2)^2)^2](https://www.gnu.org/software/gsl/doc/html/_images/math/6b8eafacce03777d28e9fb5fd8b9d4bbfc9d5c9e.png), requiring only 3 multiplications. For reasons of efficiency, these functions do not check for overflow or underflow conditions. The following is a simple example:`#include <gsl/gsl_sf_pow_int.h> /* compute 3.0**12 */ double y = gsl_sf_pow_int(3.0, 12); `

## Psi (Digamma) Function

The polygamma functions of order ![n](https://www.gnu.org/software/gsl/doc/html/_images/math/a24554f1502cf204e7aa24a6c064962c3504de48.png) are defined by

![\psi^{(n)}(x) = \left(d \over dx\right)^n \psi(x) = \left(d \over dx\right)^{n+1} \log(\Gamma(x))](https://www.gnu.org/software/gsl/doc/html/_images/math/5b70c5c4de8d79486f18a53b1bbb7e666e22ec5d.png)

where ![\psi(x) = \Gamma'(x)/\Gamma(x)](https://www.gnu.org/software/gsl/doc/html/_images/math/a0e47fd143801133c9e847f0652c09f6d828f5eb.png) is known as the digamma function. These functions are declared in the header file `gsl_sf_psi.h`.

### Digamma Function

- double `gsl_sf_psi_int`(int *n*)

- int `gsl_sf_psi_int_e`(int *n*, [gsl_sf_result](https://www.gnu.org/software/gsl/doc/html/specfunc.html#c.gsl_sf_result) * *result*)

  These routines compute the digamma function ![\psi(n)](https://www.gnu.org/software/gsl/doc/html/_images/math/81cdd27ffdc9b464d0eb45d5c0a8982c9481f4fd.png) for positive integer `n`. The digamma function is also called the Psi function.

- double `gsl_sf_psi`(double *x*)

- int `gsl_sf_psi_e`(double *x*, [gsl_sf_result](https://www.gnu.org/software/gsl/doc/html/specfunc.html#c.gsl_sf_result) * *result*)

  These routines compute the digamma function ![\psi(x)](https://www.gnu.org/software/gsl/doc/html/_images/math/007a100a79a8a0287befbd8a6ae59ef4277f0cc1.png) for general `x`, ![x \ne 0](https://www.gnu.org/software/gsl/doc/html/_images/math/41cbe33fb96924eb79bc5c017b30f8dc9cbbb606.png).

- double `gsl_sf_psi_1piy`(double *y*)

- int `gsl_sf_psi_1piy_e`(double *y*, [gsl_sf_result](https://www.gnu.org/software/gsl/doc/html/specfunc.html#c.gsl_sf_result) * *result*)

  These routines compute the real part of the digamma function on the line ![1 + i y](https://www.gnu.org/software/gsl/doc/html/_images/math/83afdf1daaa6b3394cbb7c22081cb03d122dbdc0.png), ![\Re[\psi(1 + i y)]](https://www.gnu.org/software/gsl/doc/html/_images/math/1202ef69019272953c99c3bbc7c0325ec379416a.png).

### Trigamma Function

- double `gsl_sf_psi_1_int`(int *n*)

- int `gsl_sf_psi_1_int_e`(int *n*, [gsl_sf_result](https://www.gnu.org/software/gsl/doc/html/specfunc.html#c.gsl_sf_result) * *result*)

  These routines compute the Trigamma function ![\psi'(n)](https://www.gnu.org/software/gsl/doc/html/_images/math/3adf343e470216043139785ce6afc38566a12340.png) for positive integer ![n](https://www.gnu.org/software/gsl/doc/html/_images/math/a24554f1502cf204e7aa24a6c064962c3504de48.png).

- double `gsl_sf_psi_1`(double *x*)

- int `gsl_sf_psi_1_e`(double *x*, [gsl_sf_result](https://www.gnu.org/software/gsl/doc/html/specfunc.html#c.gsl_sf_result) * *result*)

  These routines compute the Trigamma function ![\psi'(x)](https://www.gnu.org/software/gsl/doc/html/_images/math/a73e4c0ef1a777a531deb32d16b9ba3323f875f7.png) for general `x`.

### Polygamma Function

- double `gsl_sf_psi_n`(int *n*, double *x*)

- int `gsl_sf_psi_n_e`(int *n*, double *x*, [gsl_sf_result](https://www.gnu.org/software/gsl/doc/html/specfunc.html#c.gsl_sf_result) * *result*)

  These routines compute the polygamma function ![\psi^{(n)}(x)](https://www.gnu.org/software/gsl/doc/html/_images/math/5f49fb664707f7985679a30fd9d6410edf72b1db.png) for ![n \ge 0](https://www.gnu.org/software/gsl/doc/html/_images/math/a4a2e96a00e9bf71dda87e79f4f76f7ad0bca7cf.png), ![x > 0](https://www.gnu.org/software/gsl/doc/html/_images/math/fe20eedd0aef40a657ad04835ff74fa8a73efb56.png).

## Synchrotron Functions

The functions described in this section are declared in the header file `gsl_sf_synchrotron.h`.

- double `gsl_sf_synchrotron_1`(double *x*)

- int `gsl_sf_synchrotron_1_e`(double *x*, [gsl_sf_result](https://www.gnu.org/software/gsl/doc/html/specfunc.html#c.gsl_sf_result) * *result*)

  These routines compute the first synchrotron function ![x \int_x^\infty dt K_{5/3}(t)](https://www.gnu.org/software/gsl/doc/html/_images/math/bc560286ac450990cde75404d27d5f51256dbadd.png) for ![x \ge 0](https://www.gnu.org/software/gsl/doc/html/_images/math/fe351bceb81b6a8053c3711170955c252c9b2fc6.png).

- double `gsl_sf_synchrotron_2`(double *x*)

- int `gsl_sf_synchrotron_2_e`(double *x*, [gsl_sf_result](https://www.gnu.org/software/gsl/doc/html/specfunc.html#c.gsl_sf_result) * *result*)

  These routines compute the second synchrotron function ![x K_{2/3}(x)](https://www.gnu.org/software/gsl/doc/html/_images/math/873f77b2034ce3551894b77038a53662823dc073.png) for ![x \ge 0](https://www.gnu.org/software/gsl/doc/html/_images/math/fe351bceb81b6a8053c3711170955c252c9b2fc6.png).

## Transport Functions

The transport functions ![J(n,x)](https://www.gnu.org/software/gsl/doc/html/_images/math/c737003ef6b3668293b53f02fe95591324d0fe06.png) are defined by the integral representations

![J(n,x) = \int_0^x t^n e^t /(e^t - 1)^2 dt](https://www.gnu.org/software/gsl/doc/html/_images/math/c00f3e782f72b9382853d98e0821821789b206b6.png)

They are declared in the header file `gsl_sf_transport.h`.

- double `gsl_sf_transport_2`(double *x*)

- int `gsl_sf_transport_2_e`(double *x*, [gsl_sf_result](https://www.gnu.org/software/gsl/doc/html/specfunc.html#c.gsl_sf_result) * *result*)

  These routines compute the transport function ![J(2,x)](https://www.gnu.org/software/gsl/doc/html/_images/math/219d4891ce494161bc37002116067842a28736cb.png).

- double `gsl_sf_transport_3`(double *x*)

- int `gsl_sf_transport_3_e`(double *x*, [gsl_sf_result](https://www.gnu.org/software/gsl/doc/html/specfunc.html#c.gsl_sf_result) * *result*)

  These routines compute the transport function ![J(3,x)](https://www.gnu.org/software/gsl/doc/html/_images/math/66e5f3fdcbf4490995b4175e9e3f6abdbd09fed7.png).

- double `gsl_sf_transport_4`(double *x*)

- int `gsl_sf_transport_4_e`(double *x*, [gsl_sf_result](https://www.gnu.org/software/gsl/doc/html/specfunc.html#c.gsl_sf_result) * *result*)

  These routines compute the transport function ![J(4,x)](https://www.gnu.org/software/gsl/doc/html/_images/math/56963c892577065040646941ca452699a5c7050b.png).

- double `gsl_sf_transport_5`(double *x*)

- int `gsl_sf_transport_5_e`(double *x*, [gsl_sf_result](https://www.gnu.org/software/gsl/doc/html/specfunc.html#c.gsl_sf_result) * *result*)

  These routines compute the transport function ![J(5,x)](https://www.gnu.org/software/gsl/doc/html/_images/math/805cd0949b593f54fb2ce53e3278420429eca68c.png).

## Trigonometric Functions

The library includes its own trigonometric functions in order to provide consistency across platforms and reliable error estimates. These functions are declared in the header file `gsl_sf_trig.h`.

### Circular Trigonometric Functions



- double `gsl_sf_sin`(double *x*)

- int `gsl_sf_sin_e`(double *x*, [gsl_sf_result](https://www.gnu.org/software/gsl/doc/html/specfunc.html#c.gsl_sf_result) * *result*)

  These routines compute the sine function ![\sin(x)](https://www.gnu.org/software/gsl/doc/html/_images/math/35753c28a39bfd7cd9d559fb437d98a25407c6b1.png).



- double `gsl_sf_cos`(double *x*)

- int `gsl_sf_cos_e`(double *x*, [gsl_sf_result](https://www.gnu.org/software/gsl/doc/html/specfunc.html#c.gsl_sf_result) * *result*)

  These routines compute the cosine function ![\cos(x)](https://www.gnu.org/software/gsl/doc/html/_images/math/a560c690575185243b01795d2179a1492dbb5b02.png).



- double `gsl_sf_hypot`(double *x*, double *y*)

- int `gsl_sf_hypot_e`(double *x*, double *y*, [gsl_sf_result](https://www.gnu.org/software/gsl/doc/html/specfunc.html#c.gsl_sf_result) * *result*)

  These routines compute the hypotenuse function ![\sqrt{x^2 + y^2}](https://www.gnu.org/software/gsl/doc/html/_images/math/ce8d2988f650814df86140df0a851354750c0d4c.png) avoiding overflow and underflow.



- double `gsl_sf_sinc`(double *x*)

- int `gsl_sf_sinc_e`(double *x*, [gsl_sf_result](https://www.gnu.org/software/gsl/doc/html/specfunc.html#c.gsl_sf_result) * *result*)

  These routines compute ![\sinc(x) = \sin(\pi x) / (\pi x)](https://www.gnu.org/software/gsl/doc/html/_images/math/6275c833134f266211b490ce3bbe337cbf0383e9.png) for any value of `x`.

### Trigonometric Functions for Complex Arguments



- int `gsl_sf_complex_sin_e`(double *zr*, double *zi*, [gsl_sf_result](https://www.gnu.org/software/gsl/doc/html/specfunc.html#c.gsl_sf_result) * *szr*, [gsl_sf_result](https://www.gnu.org/software/gsl/doc/html/specfunc.html#c.gsl_sf_result) * *szi*)

  This function computes the complex sine, ![\sin(z_r + i z_i)](https://www.gnu.org/software/gsl/doc/html/_images/math/345b70dbacba65c78ea381e9a19b91a233e0d9c9.png) storing the real and imaginary parts in `szr`, `szi`.



- int `gsl_sf_complex_cos_e`(double *zr*, double *zi*, [gsl_sf_result](https://www.gnu.org/software/gsl/doc/html/specfunc.html#c.gsl_sf_result) * *czr*, [gsl_sf_result](https://www.gnu.org/software/gsl/doc/html/specfunc.html#c.gsl_sf_result) * *czi*)

  This function computes the complex cosine, ![\cos(z_r + i z_i)](https://www.gnu.org/software/gsl/doc/html/_images/math/8155f7272554cff656c6e1afb52de7ccb2669aa8.png) storing the real and imaginary parts in `czr`, `czi`.



- int `gsl_sf_complex_logsin_e`(double *zr*, double *zi*, [gsl_sf_result](https://www.gnu.org/software/gsl/doc/html/specfunc.html#c.gsl_sf_result) * *lszr*, [gsl_sf_result](https://www.gnu.org/software/gsl/doc/html/specfunc.html#c.gsl_sf_result) * *lszi*)

  This function computes the logarithm of the complex sine, ![\log(\sin(z_r + i z_i))](https://www.gnu.org/software/gsl/doc/html/_images/math/a1d0b333ac3fe555041bad3051ebaec27c088fc4.png) storing the real and imaginary parts in `lszr`, `lszi`.

### Hyperbolic Trigonometric Functions



- double `gsl_sf_lnsinh`(double *x*)

- int `gsl_sf_lnsinh_e`(double *x*, [gsl_sf_result](https://www.gnu.org/software/gsl/doc/html/specfunc.html#c.gsl_sf_result) * *result*)

  These routines compute ![\log(\sinh(x))](https://www.gnu.org/software/gsl/doc/html/_images/math/10bc3484e066b5c1474a6a6152c596c22096d665.png) for ![x > 0](https://www.gnu.org/software/gsl/doc/html/_images/math/fe20eedd0aef40a657ad04835ff74fa8a73efb56.png).



- double `gsl_sf_lncosh`(double *x*)

- int `gsl_sf_lncosh_e`(double *x*, [gsl_sf_result](https://www.gnu.org/software/gsl/doc/html/specfunc.html#c.gsl_sf_result) * *result*)

  These routines compute ![\log(\cosh(x))](https://www.gnu.org/software/gsl/doc/html/_images/math/1fa8b1b8d475e4291d3c68c03d197e745bc32e75.png) for any `x`.

### Conversion Functions



- int `gsl_sf_polar_to_rect`(double *r*, double *theta*, [gsl_sf_result](https://www.gnu.org/software/gsl/doc/html/specfunc.html#c.gsl_sf_result) * *x*, [gsl_sf_result](https://www.gnu.org/software/gsl/doc/html/specfunc.html#c.gsl_sf_result) * *y*)

  This function converts the polar coordinates (`r`, `theta`) to rectilinear coordinates (`x`, `y`), ![x = r\cos(\theta)](https://www.gnu.org/software/gsl/doc/html/_images/math/41e37eec630dbdcc736f052dfab95a4005f9ea72.png), ![y = r\sin(\theta)](https://www.gnu.org/software/gsl/doc/html/_images/math/473994f793f24105493d9474e6c4e4a06b47075f.png).

- int `gsl_sf_rect_to_polar`(double *x*, double *y*, [gsl_sf_result](https://www.gnu.org/software/gsl/doc/html/specfunc.html#c.gsl_sf_result) * *r*, [gsl_sf_result](https://www.gnu.org/software/gsl/doc/html/specfunc.html#c.gsl_sf_result) * *theta*)

  This function converts the rectilinear coordinates (`x`, `y`) to polar coordinates (`r`, `theta`), such that ![x = r\cos(\theta)](https://www.gnu.org/software/gsl/doc/html/_images/math/41e37eec630dbdcc736f052dfab95a4005f9ea72.png), ![y = r\sin(\theta)](https://www.gnu.org/software/gsl/doc/html/_images/math/473994f793f24105493d9474e6c4e4a06b47075f.png). The argument `theta` lies in the range ![[-\pi, \pi]](https://www.gnu.org/software/gsl/doc/html/_images/math/228eaab332fea24f5302d0e746f958f53c359816.png).

### Restriction Functions



- double `gsl_sf_angle_restrict_symm`(double *theta*)

- int `gsl_sf_angle_restrict_symm_e`(double * *theta*)

  These routines force the angle `theta` to lie in the range ![(-\pi,\pi]](https://www.gnu.org/software/gsl/doc/html/_images/math/2c4cbfed45ff565dab518bf9e6dcdb1558490157.png).Note that the mathematical value of ![\pi](https://www.gnu.org/software/gsl/doc/html/_images/math/e619909e385649c40e3610ab5d7f9f8643e8b14b.png) is slightly greater than `M_PI`, so the machine numbers `M_PI` and `-M_PI` are included in the range.

- double `gsl_sf_angle_restrict_pos`(double *theta*)

- int `gsl_sf_angle_restrict_pos_e`(double * *theta*)

  These routines force the angle `theta` to lie in the range ![[0, 2\pi)](https://www.gnu.org/software/gsl/doc/html/_images/math/02a7e538500116dd51b2009667087a97c04b94c5.png).Note that the mathematical value of ![2\pi](https://www.gnu.org/software/gsl/doc/html/_images/math/059fbd5d7d78643be5cf4a37b45b0f8140e2d605.png) is slightly greater than `2*M_PI`, so the machine number `2*M_PI` is included in the range.

### Trigonometric Functions With Error Estimates

- int `gsl_sf_sin_err_e`(double *x*, double *dx*, [gsl_sf_result](https://www.gnu.org/software/gsl/doc/html/specfunc.html#c.gsl_sf_result) * *result*)

  This routine computes the sine of an angle `x` with an associated absolute error `dx`, ![\sin(x \pm dx)](https://www.gnu.org/software/gsl/doc/html/_images/math/6b1c6693e1413728b368e8c57cb9f53cb4336b8e.png). Note that this function is provided in the error-handling form only since its purpose is to compute the propagated error.

- int `gsl_sf_cos_err_e`(double *x*, double *dx*, [gsl_sf_result](https://www.gnu.org/software/gsl/doc/html/specfunc.html#c.gsl_sf_result) * *result*)

  This routine computes the cosine of an angle `x` with an associated absolute error `dx`, ![\cos(x \pm dx)](https://www.gnu.org/software/gsl/doc/html/_images/math/aecc22b3505a14adfe334e8236ee6105ea6615e2.png). Note that this function is provided in the error-handling form only since its purpose is to compute the propagated error.

## Zeta Functions

The Riemann zeta function is defined in Abramowitz & Stegun, Section 23.2. The functions described in this section are declared in the header file `gsl_sf_zeta.h`.

### Riemann Zeta Function

The Riemann zeta function is defined by the infinite sum

![\zeta(s) = \sum_{k=1}^\infty k^{-s}](https://www.gnu.org/software/gsl/doc/html/_images/math/4bafbaf9387636fea0aec513ccc71658a8705c21.png)

- double `gsl_sf_zeta_int`(int *n*)

- int `gsl_sf_zeta_int_e`(int *n*, [gsl_sf_result](https://www.gnu.org/software/gsl/doc/html/specfunc.html#c.gsl_sf_result) * *result*)

  These routines compute the Riemann zeta function ![\zeta(n)](https://www.gnu.org/software/gsl/doc/html/_images/math/6e95a82fa88afe07e0bc5bc1ad4058d12ad8d9d8.png) for integer `n`, ![n \ne 1](https://www.gnu.org/software/gsl/doc/html/_images/math/dd56bdc91c3ccd978e89c874cf0d4c3de8b56aaa.png).

- double `gsl_sf_zeta`(double *s*)

- int `gsl_sf_zeta_e`(double *s*, [gsl_sf_result](https://www.gnu.org/software/gsl/doc/html/specfunc.html#c.gsl_sf_result) * *result*)

  These routines compute the Riemann zeta function ![\zeta(s)](https://www.gnu.org/software/gsl/doc/html/_images/math/578516529b1ebc9108740411102af0e4bd2560ea.png) for arbitrary `s`, ![s \ne 1](https://www.gnu.org/software/gsl/doc/html/_images/math/f8f52a732fd6f03db2d92e2efaa00ee3304f95ae.png).

### Riemann Zeta Function Minus One

For large positive argument, the Riemann zeta function approaches one. In this region the fractional part is interesting, and therefore we need a function to evaluate it explicitly.

- double `gsl_sf_zetam1_int`(int *n*)

- int `gsl_sf_zetam1_int_e`(int *n*, [gsl_sf_result](https://www.gnu.org/software/gsl/doc/html/specfunc.html#c.gsl_sf_result) * *result*)

  These routines compute ![\zeta(n) - 1](https://www.gnu.org/software/gsl/doc/html/_images/math/09e1667faaa691f41c26af671cdd4c3aee28697f.png) for integer `n`, ![n \ne 1](https://www.gnu.org/software/gsl/doc/html/_images/math/dd56bdc91c3ccd978e89c874cf0d4c3de8b56aaa.png).

- double `gsl_sf_zetam1`(double *s*)

- int `gsl_sf_zetam1_e`(double *s*, [gsl_sf_result](https://www.gnu.org/software/gsl/doc/html/specfunc.html#c.gsl_sf_result) * *result*)

  These routines compute ![\zeta(s) - 1](https://www.gnu.org/software/gsl/doc/html/_images/math/4892d51287b937f55b8d2901057aa36117a758b5.png) for arbitrary `s`, ![s \ne 1](https://www.gnu.org/software/gsl/doc/html/_images/math/f8f52a732fd6f03db2d92e2efaa00ee3304f95ae.png).

### Hurwitz Zeta Function

The Hurwitz zeta function is defined by

![\zeta(s,q) = \sum_0^\infty (k+q)^{-s}](https://www.gnu.org/software/gsl/doc/html/_images/math/dc82b96ff62566ebdb4eb42a7edf2b3589f6909c.png)

- double `gsl_sf_hzeta`(double *s*, double *q*)

- int `gsl_sf_hzeta_e`(double *s*, double *q*, [gsl_sf_result](https://www.gnu.org/software/gsl/doc/html/specfunc.html#c.gsl_sf_result) * *result*)

  These routines compute the Hurwitz zeta function ![\zeta(s,q)](https://www.gnu.org/software/gsl/doc/html/_images/math/5dd1b1a409c2e3ce69160dc6c98663cd4854ba2b.png) for ![s > 1](https://www.gnu.org/software/gsl/doc/html/_images/math/6f0998d9e1aedc1393282d6ff8038339837f0475.png), ![q > 0](https://www.gnu.org/software/gsl/doc/html/_images/math/13beb6f1c51ffc7d09c8c8b90260bcab716f46e2.png).

### Eta Function

The eta function is defined by

![\eta(s) = (1-2^{1-s}) \zeta(s)](https://www.gnu.org/software/gsl/doc/html/_images/math/7aade218a1708ed53a3280a2826d784f564b8388.png)

- double `gsl_sf_eta_int`(int *n*)

- int `gsl_sf_eta_int_e`(int *n*, [gsl_sf_result](https://www.gnu.org/software/gsl/doc/html/specfunc.html#c.gsl_sf_result) * *result*)

  These routines compute the eta function ![\eta(n)](https://www.gnu.org/software/gsl/doc/html/_images/math/cae2c18e702ef14d37fb636e9472c0d763562778.png) for integer `n`.

- double `gsl_sf_eta`(double *s*)

- int `gsl_sf_eta_e`(double *s*, [gsl_sf_result](https://www.gnu.org/software/gsl/doc/html/specfunc.html#c.gsl_sf_result) * *result*)

  These routines compute the eta function ![\eta(s)](https://www.gnu.org/software/gsl/doc/html/_images/math/21b545f1fa002873a9e86afeabd5fd1f2b56ca82.png) for arbitrary `s`.

## Examples

The following example demonstrates the use of the error handling form of the special functions, in this case to compute the Bessel function ![J_0(5.0)](https://www.gnu.org/software/gsl/doc/html/_images/math/79b4a8e75e07d9fc4bcbb411c25179e70b864399.png),

```
#include <stdio.h>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_sf_bessel.h>

int
main (void)
{
  double x = 5.0;
  gsl_sf_result result;

  double expected = -0.17759677131433830434739701;

  int status = gsl_sf_bessel_J0_e (x, &result);

  printf ("status  = %s\n", gsl_strerror(status));
  printf ("J0(5.0) = %.18f\n"
          "      +/- % .18f\n",
          result.val, result.err);
  printf ("exact   = %.18f\n", expected);
  return status;
}
```

Here are the results of running the program,

```
status  = success
J0(5.0) = -0.177596771314338264
      +/-  0.000000000000000193
exact   = -0.177596771314338292
```

The next program computes the same quantity using the natural form of the function. In this case the error term `result.err` and return status are not accessible.

```
#include <stdio.h>
#include <gsl/gsl_sf_bessel.h>

int
main (void)
{
  double x = 5.0;
  double expected = -0.17759677131433830434739701;

  double y = gsl_sf_bessel_J0 (x);

  printf ("J0(5.0) = %.18f\n", y);
  printf ("exact   = %.18f\n", expected);
  return 0;
}
```

The results of the function are the same,

```
J0(5.0) = -0.177596771314338264
exact   = -0.177596771314338292
```

## References and Further Reading

The library follows the conventions of the following book where possible,

- Handbook of Mathematical Functions, edited by Abramowitz & Stegun, Dover, ISBN 0486612724.

The following papers contain information on the algorithms used to compute the special functions,

- Allan J. MacLeod, MISCFUN: A software package to compute uncommon special functions. ACM Trans. Math. Soft., vol.: 22, 1996, 288–301
- G.N. Watson, A Treatise on the Theory of Bessel Functions, 2nd Edition (Cambridge University Press, 1944).
- G. Nemeth, Mathematical Approximations of Special Functions, Nova Science Publishers, ISBN 1-56072-052-2
- B.C. Carlson, Special Functions of Applied Mathematics (1977)
- N. M. Temme, Special Functions: An Introduction to the Classical Functions of Mathematical Physics (1996), ISBN 978-0471113133.
- W.J. Thompson, Atlas for Computing Mathematical Functions, John Wiley & Sons, New York (1997).
- Y.Y. Luke, Algorithms for the Computation of Mathematical Functions, Academic Press, New York (1977).
- S. A. Holmes and W. E. Featherstone, A unified approach to the Clenshaw summation and the recursive computation of very high degree and order normalised associated Legendre functions, Journal of Geodesy, 76, pg. 279-299, 2002.