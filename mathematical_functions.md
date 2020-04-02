# Mathematical Functions

This chapter describes basic mathematical functions. Some of these functions are present in system libraries, but the alternative versions given here can be used as a substitute when the system functions are not available.

The functions and macros described in this chapter are defined in the header file `gsl_math.h`.



## Mathematical Constants

The library ensures that the standard BSD mathematical constants are defined. For reference, here is a list of the constants:

| `M_E`        | The base of exponentials, ![e](https://www.gnu.org/software/gsl/doc/html/_images/math/9a7ff5d0469d9a20a6a7da1b2c1ab9b1014d53c9.png) |
| ------------ | ------------------------------------------------------------ |
| `M_LOG2E`    | The base-2 logarithm of ![e](https://www.gnu.org/software/gsl/doc/html/_images/math/9a7ff5d0469d9a20a6a7da1b2c1ab9b1014d53c9.png), ![\log_2 (e)](https://www.gnu.org/software/gsl/doc/html/_images/math/3e4536235dc561151e9100418c509cf6ec7eaf44.png) |
| `M_LOG10E`   | The base-10 logarithm of ![e](https://www.gnu.org/software/gsl/doc/html/_images/math/9a7ff5d0469d9a20a6a7da1b2c1ab9b1014d53c9.png), ![\log_{10} (e)](https://www.gnu.org/software/gsl/doc/html/_images/math/08d24dd899b44691682c216724c9479fe7f71803.png) |
| `M_SQRT2`    | The square root of two, ![\sqrt 2](https://www.gnu.org/software/gsl/doc/html/_images/math/dbc6bd271a82149b54abaf21ba8e7dce06d0f27c.png) |
| `M_SQRT1_2`  | The square root of one-half, ![\sqrt{1/2}](https://www.gnu.org/software/gsl/doc/html/_images/math/0974de440d2a16dd32af8e2c9b02855bcb6673a3.png) |
| `M_SQRT3`    | The square root of three, ![\sqrt 3](https://www.gnu.org/software/gsl/doc/html/_images/math/747d29fce59b708e06e9c8e8e1cfab682d05cbef.png) |
| `M_PI`       | The constant pi, ![\pi](https://www.gnu.org/software/gsl/doc/html/_images/math/e619909e385649c40e3610ab5d7f9f8643e8b14b.png) |
| `M_PI_2`     | Pi divided by two, ![\pi/2](https://www.gnu.org/software/gsl/doc/html/_images/math/dc93cceb6e7be84e1d3470cd32db7aa763298782.png) |
| `M_PI_4`     | Pi divided by four, ![\pi/4](https://www.gnu.org/software/gsl/doc/html/_images/math/43c43fc5b07afe32ae8a42b6edc70240d606a959.png) |
| `M_SQRTPI`   | The square root of pi, ![\sqrt\pi](https://www.gnu.org/software/gsl/doc/html/_images/math/a2ce252654bd9a2e1e3e679cd7305bbc67e35b67.png) |
| `M_2_SQRTPI` | Two divided by the square root of pi, ![2/\sqrt\pi](https://www.gnu.org/software/gsl/doc/html/_images/math/824a45cdf917efc64ce60f3d7a780fa3276b2e97.png) |
| `M_1_PI`     | The reciprocal of pi, ![1/\pi](https://www.gnu.org/software/gsl/doc/html/_images/math/5527809086314c0f9baf85da3c4c883075eaea84.png) |
| `M_2_PI`     | Twice the reciprocal of pi, ![2/\pi](https://www.gnu.org/software/gsl/doc/html/_images/math/26f4dfcb945d580e1e8e88569db63fb8d67b9c7b.png) |
| `M_LN10`     | The natural logarithm of ten, ![\ln(10)](https://www.gnu.org/software/gsl/doc/html/_images/math/7d01a53757b4fbbd6a89fdddb25186c50e85c0b3.png) |
| `M_LN2`      | The natural logarithm of two, ![\ln(2)](https://www.gnu.org/software/gsl/doc/html/_images/math/550d78bce5747b1ee1ded694dc87981eca5230d5.png) |
| `M_LNPI`     | The natural logarithm of pi, ![\ln(\pi)](https://www.gnu.org/software/gsl/doc/html/_images/math/9a64777a97f3fab686e7777fac01d471c211b7f2.png) |
| `M_EULER`    | Euler’s constant, ![\gamma](https://www.gnu.org/software/gsl/doc/html/_images/math/3ad522562ec56d0efcab960b91eeeda7e974cdf2.png) |



## Infinities and Not-a-number

- `GSL_POSINF`

  This macro contains the IEEE representation of positive infinity, ![+\infty](https://www.gnu.org/software/gsl/doc/html/_images/math/6be41aa884c6c25ee576afe5476eea358adc2931.png). It is computed from the expression `+1.0/0.0`.

- `GSL_NEGINF`

  This macro contains the IEEE representation of negative infinity, ![-\infty](https://www.gnu.org/software/gsl/doc/html/_images/math/342dc8da9165373fac8f4b9024eea6d6e76c53d9.png). It is computed from the expression `-1.0/0.0`.



- `GSL_NAN`

  This macro contains the IEEE representation of the Not-a-Number symbol, `NaN`. It is computed from the ratio `0.0/0.0`.

- int `gsl_isnan`(const double *x*)

  This function returns 1 if `x` is not-a-number.

- int `gsl_isinf`(const double *x*)

  This function returns ![+1](https://www.gnu.org/software/gsl/doc/html/_images/math/a7bc7432fa358292afd63292ed1a8acd0553012b.png) if `x` is positive infinity, ![-1](https://www.gnu.org/software/gsl/doc/html/_images/math/e26aecc9561355dd0e6fbf2742d0ff0d6b58c864.png) if `x` is negative infinity and 0 otherwise. [[1\]](https://www.gnu.org/software/gsl/doc/html/math.html#f1)

- int `gsl_finite`(const double *x*)

  This function returns 1 if `x` is a real number, and 0 if it is infinite or not-a-number.

## Elementary Functions

The following routines provide portable implementations of functions found in the BSD math library. When native versions are not available the functions described here can be used instead. The substitution can be made automatically if you use `autoconf` to compile your application (see [Portability functions](https://www.gnu.org/software/gsl/doc/html/usage.html#portability-functions)).



- double `gsl_log1p`(const double *x*)

  This function computes the value of ![\log(1+x)](https://www.gnu.org/software/gsl/doc/html/_images/math/2b09c0e82a24902ab2522d2763cae03560f7074d.png) in a way that is accurate for small `x`. It provides an alternative to the BSD math function `log1p(x)`.



- double `gsl_expm1`(const double *x*)

  This function computes the value of ![\exp(x)-1](https://www.gnu.org/software/gsl/doc/html/_images/math/40d15d7dfac571e38c04c5f83ac3cb15c43ecff5.png) in a way that is accurate for small `x`. It provides an alternative to the BSD math function `expm1(x)`.



- double `gsl_hypot`(const double *x*, const double *y*)

  This function computes the value of ![\sqrt{x^2 + y^2}](https://www.gnu.org/software/gsl/doc/html/_images/math/ce8d2988f650814df86140df0a851354750c0d4c.png) in a way that avoids overflow. It provides an alternative to the BSD math function `hypot(x,y)`.



- double `gsl_hypot3`(const double *x*, const double *y*, const double *z*)

  This function computes the value of ![\sqrt{x^2 + y^2 + z^2}](https://www.gnu.org/software/gsl/doc/html/_images/math/323db571b2702ba9b7a914fc773d85c542a0ff97.png) in a way that avoids overflow.



- double `gsl_acosh`(const double *x*)

  This function computes the value of ![\arccosh{(x)}](https://www.gnu.org/software/gsl/doc/html/_images/math/5fa0751db7c465784cc11049a02bc2af37d23c85.png). It provides an alternative to the standard math function `acosh(x)`.



- double `gsl_asinh`(const double *x*)

  This function computes the value of ![\arcsinh{(x)}](https://www.gnu.org/software/gsl/doc/html/_images/math/d7606577b734d9c9e2c0c1eb20019006400d59b4.png). It provides an alternative to the standard math function `asinh(x)`.



- double `gsl_atanh`(const double *x*)

  This function computes the value of ![\arctanh{(x)}](https://www.gnu.org/software/gsl/doc/html/_images/math/93050438540a58dbef1a345bc4d8fbe7f6bfe211.png). It provides an alternative to the standard math function `atanh(x)`.



- double `gsl_ldexp`(double *x*, int *e*)

  This function computes the value of ![x * 2^e](https://www.gnu.org/software/gsl/doc/html/_images/math/f10e65fc6ac0148184a687098b60631c9d73a148.png). It provides an alternative to the standard math function `ldexp(x,e)`.



- double `gsl_frexp`(double *x*, int * *e*)

  This function splits the number `x` into its normalized fraction ![f](https://www.gnu.org/software/gsl/doc/html/_images/math/d947f897010fa03fc9854d6af00264296a50930c.png) and exponent ![e](https://www.gnu.org/software/gsl/doc/html/_images/math/9a7ff5d0469d9a20a6a7da1b2c1ab9b1014d53c9.png), such that ![x = f * 2^e](https://www.gnu.org/software/gsl/doc/html/_images/math/4265da59f7bdd86698eb37b66c1fb695c7c95092.png) and ![0.5 <= f < 1](https://www.gnu.org/software/gsl/doc/html/_images/math/8f5325f4a80d53f6c88cbff0444d5604b8bffa0b.png). The function returns ![f](https://www.gnu.org/software/gsl/doc/html/_images/math/d947f897010fa03fc9854d6af00264296a50930c.png) and stores the exponent in ![e](https://www.gnu.org/software/gsl/doc/html/_images/math/9a7ff5d0469d9a20a6a7da1b2c1ab9b1014d53c9.png). If ![x](https://www.gnu.org/software/gsl/doc/html/_images/math/1c5715c8edd0bbc330baa0f0494f1f7951619754.png) is zero, both ![f](https://www.gnu.org/software/gsl/doc/html/_images/math/d947f897010fa03fc9854d6af00264296a50930c.png) and ![e](https://www.gnu.org/software/gsl/doc/html/_images/math/9a7ff5d0469d9a20a6a7da1b2c1ab9b1014d53c9.png) are set to zero. This function provides an alternative to the standard math function `frexp(x, e)`.

## Small integer powers

A common complaint about the standard C library is its lack of a function for calculating (small) integer powers. GSL provides some simple functions to fill this gap. For reasons of efficiency, these functions do not check for overflow or underflow conditions.

- double `gsl_pow_int`(double *x*, int *n*)

- double `gsl_pow_uint`(double *x*, unsigned int *n*)

  These routines computes the power ![x^n](https://www.gnu.org/software/gsl/doc/html/_images/math/ca5dc9582574e5893064ffb5d9f66578c649ffb7.png) for integer `n`. The power is computed efficiently—for example, ![x^8](https://www.gnu.org/software/gsl/doc/html/_images/math/b77674817382f8d63ad5d32d451ab476be4c354e.png) is computed as ![((x^2)^2)^2](https://www.gnu.org/software/gsl/doc/html/_images/math/6b8eafacce03777d28e9fb5fd8b9d4bbfc9d5c9e.png), requiring only 3 multiplications. A version of this function which also computes the numerical error in the result is available as [`gsl_sf_pow_int_e()`](https://www.gnu.org/software/gsl/doc/html/specfunc.html#c.gsl_sf_pow_int_e).

- double `gsl_pow_2`(const double *x*)

- double `gsl_pow_3`(const double *x*)

- double `gsl_pow_4`(const double *x*)

- double `gsl_pow_5`(const double *x*)

- double `gsl_pow_6`(const double *x*)

- double `gsl_pow_7`(const double *x*)

- double `gsl_pow_8`(const double *x*)

- double `gsl_pow_9`(const double *x*)

  These functions can be used to compute small integer powers ![x^2](https://www.gnu.org/software/gsl/doc/html/_images/math/85f34d6e6902d5250565503ca969ae7846fabdba.png), ![x^3](https://www.gnu.org/software/gsl/doc/html/_images/math/2572f46f9499ab927e2beb944b98db6c30567910.png), etc. efficiently. The functions will be inlined when `HAVE_INLINE` is defined, so that use of these functions should be as efficient as explicitly writing the corresponding product expression:`#include <gsl/gsl_math.h> double y = gsl_pow_4 (3.141)  /* compute 3.141**4 */ `

## Testing the Sign of Numbers

- `GSL_SIGN`(x)

  This macro returns the sign of `x`. It is defined as `((x) >= 0 ? 1 : -1)`. Note that with this definition the sign of zero is positive (regardless of its IEEE sign bit).

## Testing for Odd and Even Numbers

- `GSL_IS_ODD`(n)

  This macro evaluates to 1 if `n` is odd and 0 if `n` is even. The argument `n` must be of integer type.

- `GSL_IS_EVEN`(n)

  This macro is the opposite of [`GSL_IS_ODD`](https://www.gnu.org/software/gsl/doc/html/math.html#c.GSL_IS_ODD). It evaluates to 1 if `n` is even and 0 if `n` is odd. The argument `n` must be of integer type.

## Maximum and Minimum functions

Note that the following macros perform multiple evaluations of their arguments, so they should not be used with arguments that have side effects (such as a call to a random number generator).



- `GSL_MAX`(a, b)

  This macro returns the maximum of `a` and `b`. It is defined as `((a) > (b) ? (a):(b))`.



- `GSL_MIN`(a, b)

  This macro returns the minimum of `a` and `b`. It is defined as `((a) < (b) ? (a):(b))`.

- extern inline double `GSL_MAX_DBL`(double *a*, double *b*)

  This function returns the maximum of the double precision numbers `a` and `b` using an inline function. The use of a function allows for type checking of the arguments as an extra safety feature. On platforms where inline functions are not available the macro [`GSL_MAX`](https://www.gnu.org/software/gsl/doc/html/math.html#c.GSL_MAX) will be automatically substituted.

- extern inline double `GSL_MIN_DBL`(double *a*, double *b*)

  This function returns the minimum of the double precision numbers `a` and `b` using an inline function. The use of a function allows for type checking of the arguments as an extra safety feature. On platforms where inline functions are not available the macro [`GSL_MIN`](https://www.gnu.org/software/gsl/doc/html/math.html#c.GSL_MIN) will be automatically substituted.

- extern inline int `GSL_MAX_INT`(int *a*, int *b*)

- extern inline int `GSL_MIN_INT`(int *a*, int *b*)

  These functions return the maximum or minimum of the integers `a` and `b` using an inline function. On platforms where inline functions are not available the macros [`GSL_MAX`](https://www.gnu.org/software/gsl/doc/html/math.html#c.GSL_MAX) or [`GSL_MIN`](https://www.gnu.org/software/gsl/doc/html/math.html#c.GSL_MIN)will be automatically substituted.

- extern inline long double `GSL_MAX_LDBL`(long double *a*, long double *b*)

- extern inline long double `GSL_MIN_LDBL`(long double *a*, long double *b*)

  These functions return the maximum or minimum of the long doubles `a` and `b` using an inline function. On platforms where inline functions are not available the macros [`GSL_MAX`](https://www.gnu.org/software/gsl/doc/html/math.html#c.GSL_MAX) or [`GSL_MIN`](https://www.gnu.org/software/gsl/doc/html/math.html#c.GSL_MIN)will be automatically substituted.

## Approximate Comparison of Floating Point Numbers

It is sometimes useful to be able to compare two floating point numbers approximately, to allow for rounding and truncation errors. The following function implements the approximate floating-point comparison algorithm proposed by D.E. Knuth in Section 4.2.2 of “Seminumerical Algorithms” (3rd edition).



- int `gsl_fcmp`(double *x*, double *y*, double *epsilon*)

  This function determines whether `x` and `y` are approximately equal to a relative accuracy `epsilon`.The relative accuracy is measured using an interval of size ![2 \delta](https://www.gnu.org/software/gsl/doc/html/_images/math/9ec07d1b6fd692c29c7cd0aee52e724f153d6da4.png), where ![\delta = 2^k \epsilon](https://www.gnu.org/software/gsl/doc/html/_images/math/cf42b9ac4fad0ddba7dd198c568ca34d3e78c3f5.png) and ![k](https://www.gnu.org/software/gsl/doc/html/_images/math/90645b5c2abc13d7de104fd87abeb2f19406e5f8.png) is the maximum base-2 exponent of ![x](https://www.gnu.org/software/gsl/doc/html/_images/math/1c5715c8edd0bbc330baa0f0494f1f7951619754.png) and ![y](https://www.gnu.org/software/gsl/doc/html/_images/math/b3555326e18542ce3a7ddc86abbd1aa293204e4e.png) as computed by the function `frexp()`.If ![x](https://www.gnu.org/software/gsl/doc/html/_images/math/1c5715c8edd0bbc330baa0f0494f1f7951619754.png) and ![y](https://www.gnu.org/software/gsl/doc/html/_images/math/b3555326e18542ce3a7ddc86abbd1aa293204e4e.png) lie within this interval, they are considered approximately equal and the function returns 0. Otherwise if ![x < y](https://www.gnu.org/software/gsl/doc/html/_images/math/44d74e421cb7c2fdf7a5d4a5d3d435fa5cb1093b.png), the function returns ![-1](https://www.gnu.org/software/gsl/doc/html/_images/math/e26aecc9561355dd0e6fbf2742d0ff0d6b58c864.png), or if ![x > y](https://www.gnu.org/software/gsl/doc/html/_images/math/4569a01f975ddbb6732f359394b14c1189cf3092.png), the function returns ![+1](https://www.gnu.org/software/gsl/doc/html/_images/math/a7bc7432fa358292afd63292ed1a8acd0553012b.png).Note that ![x](https://www.gnu.org/software/gsl/doc/html/_images/math/1c5715c8edd0bbc330baa0f0494f1f7951619754.png) and ![y](https://www.gnu.org/software/gsl/doc/html/_images/math/b3555326e18542ce3a7ddc86abbd1aa293204e4e.png) are compared to relative accuracy, so this function is not suitable for testing whether a value is approximately zero.The implementation is based on the package `fcmp` by T.C. Belding.

### Footnotes

[[1]](https://www.gnu.org/software/gsl/doc/html/math.html#id1)Note that the C99 standard only requires the system `isinf()` function to return a non-zero value, without the sign of the infinity. The implementation in some earlier versions of GSL used the system `isinf()` function and may have this behavior on some platforms. Therefore, it is advisable to test the sign of `x` separately, if needed, rather than relying the sign of the return value from [`gsl_isinf()`](https://www.gnu.org/software/gsl/doc/html/math.html#c.gsl_isinf).

