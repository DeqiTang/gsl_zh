# Complex Numbers

The functions described in this chapter provide support for complex numbers. The algorithms take care to avoid unnecessary intermediate underflows and overflows, allowing the functions to be evaluated over as much of the complex plane as possible.

For multiple-valued functions the branch cuts have been chosen to follow the conventions of Abramowitz and Stegun. The functions return principal values which are the same as those in GNU Calc, which in turn are the same as those in “Common Lisp, The Language (Second Edition)” [[1\]](https://www.gnu.org/software/gsl/doc/html/complex.html#f1) and the HP-28/48 series of calculators.

The complex types are defined in the header file `gsl_complex.h`, while the corresponding complex functions and arithmetic operations are defined in `gsl_complex_math.h`.



## Representation of complex numbers

Complex numbers are represented using the type `gsl_complex`. The internal representation of this type may vary across platforms and should not be accessed directly. The functions and macros described below allow complex numbers to be manipulated in a portable way.

For reference, the default form of the `gsl_complex` type is given by the following struct:

```
typedef struct
{
  double dat[2];
} gsl_complex;
```

The real and imaginary part are stored in contiguous elements of a two element array. This eliminates any padding between the real and imaginary parts, `dat[0]` and `dat[1]`, allowing the struct to be mapped correctly onto packed complex arrays.

- gsl_complex `gsl_complex_rect`(double *x*, double *y*)

  This function uses the rectangular Cartesian components ![(x,y)](https://www.gnu.org/software/gsl/doc/html/_images/math/7ecf269c4c3ea45e59922da2237627de837dbc5f.png) to return the complex number ![z = x + i y](https://www.gnu.org/software/gsl/doc/html/_images/math/040bcc5e5b070ab2f73a5ccbe84f3fce1aabc65e.png). An inline version of this function is used when `HAVE_INLINE` is defined.

- gsl_complex `gsl_complex_polar`(double *r*, double *theta*)

  This function returns the complex number ![z = r \exp(i \theta) = r (\cos(\theta) + i \sin(\theta))](https://www.gnu.org/software/gsl/doc/html/_images/math/cb9e970effd53d28aeeb490be17705a42ad3484c.png) from the polar representation (`r`, `theta`).

- `GSL_REAL`(z)

- `GSL_IMAG`(z)

  These macros return the real and imaginary parts of the complex number `z`.

- `GSL_SET_COMPLEX`(zp, x, y)

  This macro uses the Cartesian components (`x`, `y`) to set the real and imaginary parts of the complex number pointed to by `zp`. For example:`GSL_SET_COMPLEX(&z, 3, 4) `sets ![z](https://www.gnu.org/software/gsl/doc/html/_images/math/4c9e578e25d40c5b0c5cad91a1e2f34e58b5be5e.png) to be ![3 + 4i](https://www.gnu.org/software/gsl/doc/html/_images/math/7c556d38361e52d0fc0b99072d8eef0ab3582b1c.png).

- `GSL_SET_REAL`(zp, x)

- `GSL_SET_IMAG`(zp, y)

  These macros allow the real and imaginary parts of the complex number pointed to by `zp` to be set independently.

## Properties of complex numbers



- double `gsl_complex_arg`(gsl_complex *z*)

  This function returns the argument of the complex number `z`, ![\arg(z)](https://www.gnu.org/software/gsl/doc/html/_images/math/4e987623320522533af6916fd597784c19cd21df.png), where ![-\pi < \arg(z) <= \pi](https://www.gnu.org/software/gsl/doc/html/_images/math/7917e9eac29126dba6156e03115199684bfbdc22.png).



- double `gsl_complex_abs`(gsl_complex *z*)

  This function returns the magnitude of the complex number `z`, ![|z|](https://www.gnu.org/software/gsl/doc/html/_images/math/e21e85b701e7cd55df38706041140e741b3ab735.png).

- double `gsl_complex_abs2`(gsl_complex *z*)

  This function returns the squared magnitude of the complex number `z`, ![|z|^2](https://www.gnu.org/software/gsl/doc/html/_images/math/d49a27f7462df28deafdb13d71eb343bbc47d7b3.png).

- double `gsl_complex_logabs`(gsl_complex *z*)

  This function returns the natural logarithm of the magnitude of the complex number `z`, ![\log|z|](https://www.gnu.org/software/gsl/doc/html/_images/math/45eedc748dd68c161b2ecdf1f976fed4538a39fa.png). It allows an accurate evaluation of ![\log|z|](https://www.gnu.org/software/gsl/doc/html/_images/math/45eedc748dd68c161b2ecdf1f976fed4538a39fa.png) when ![|z|](https://www.gnu.org/software/gsl/doc/html/_images/math/e21e85b701e7cd55df38706041140e741b3ab735.png) is close to one. The direct evaluation of `log(gsl_complex_abs(z))` would lead to a loss of precision in this case.



## Complex arithmetic operators

- gsl_complex `gsl_complex_add`(gsl_complex *a*, gsl_complex *b*)

  This function returns the sum of the complex numbers `a` and `b`, ![z=a+b](https://www.gnu.org/software/gsl/doc/html/_images/math/982c9a8b117a69f5ec10fdc9503cdffccf0e28cf.png).

- gsl_complex `gsl_complex_sub`(gsl_complex *a*, gsl_complex *b*)

  This function returns the difference of the complex numbers `a` and `b`, ![z=a-b](https://www.gnu.org/software/gsl/doc/html/_images/math/220e9dd427a240cbb17c3fa910429966c7684204.png).

- gsl_complex `gsl_complex_mul`(gsl_complex *a*, gsl_complex *b*)

  This function returns the product of the complex numbers `a` and `b`, ![z=ab](https://www.gnu.org/software/gsl/doc/html/_images/math/6011addb64f618ad2a89ba17a03b7a8af1ddbaaa.png).

- gsl_complex `gsl_complex_div`(gsl_complex *a*, gsl_complex *b*)

  This function returns the quotient of the complex numbers `a` and `b`, ![z=a/b](https://www.gnu.org/software/gsl/doc/html/_images/math/2942d718132347e1d85747db237318ed53546b8e.png).

- gsl_complex `gsl_complex_add_real`(gsl_complex *a*, double *x*)

  This function returns the sum of the complex number `a` and the real number `x`, ![z=a+x](https://www.gnu.org/software/gsl/doc/html/_images/math/e60b951e63fd9050b2e2342b223cc2c148c6130a.png).

- gsl_complex `gsl_complex_sub_real`(gsl_complex *a*, double *x*)

  This function returns the difference of the complex number `a` and the real number `x`, ![z=a-x](https://www.gnu.org/software/gsl/doc/html/_images/math/05a1c7c887eb40a19bbcf9442746560eff43cceb.png).

- gsl_complex `gsl_complex_mul_real`(gsl_complex *a*, double *x*)

  This function returns the product of the complex number `a` and the real number `x`, ![z=ax](https://www.gnu.org/software/gsl/doc/html/_images/math/a8a2cbb20023a8e6da57ecf02259766b067ded61.png).

- gsl_complex `gsl_complex_div_real`(gsl_complex *a*, double *x*)

  This function returns the quotient of the complex number `a` and the real number `x`, ![z=a/x](https://www.gnu.org/software/gsl/doc/html/_images/math/39418e9ce3c2c1050d08234d3c80e239ff611c97.png).

- gsl_complex `gsl_complex_add_imag`(gsl_complex *a*, double *y*)

  This function returns the sum of the complex number `a` and the imaginary number ![iy](https://www.gnu.org/software/gsl/doc/html/_images/math/4a7de0628a5dbec916dfa40ac8e749ebb1a1fc6d.png), ![z=a+iy](https://www.gnu.org/software/gsl/doc/html/_images/math/86a5743fcd9cc8a731042f12fe96f093b7206601.png).

- gsl_complex `gsl_complex_sub_imag`(gsl_complex *a*, double *y*)

  This function returns the difference of the complex number `a` and the imaginary number ![iy](https://www.gnu.org/software/gsl/doc/html/_images/math/4a7de0628a5dbec916dfa40ac8e749ebb1a1fc6d.png), ![z=a-iy](https://www.gnu.org/software/gsl/doc/html/_images/math/0123334294a942f857a4f068e4f14a2e35cad287.png).

- gsl_complex `gsl_complex_mul_imag`(gsl_complex *a*, double *y*)

  This function returns the product of the complex number `a` and the imaginary number ![iy](https://www.gnu.org/software/gsl/doc/html/_images/math/4a7de0628a5dbec916dfa40ac8e749ebb1a1fc6d.png), ![z=a*(iy)](https://www.gnu.org/software/gsl/doc/html/_images/math/6036689ce613378cc7e6fad9e07880182f1a9708.png).

- gsl_complex `gsl_complex_div_imag`(gsl_complex *a*, double *y*)

  This function returns the quotient of the complex number `a` and the imaginary number ![iy](https://www.gnu.org/software/gsl/doc/html/_images/math/4a7de0628a5dbec916dfa40ac8e749ebb1a1fc6d.png), ![z=a/(iy)](https://www.gnu.org/software/gsl/doc/html/_images/math/44135db4322530533f10b1356e5d55be63bcba13.png).



- gsl_complex `gsl_complex_conjugate`(gsl_complex *z*)

  This function returns the complex conjugate of the complex number `z`, ![z^* = x - i y](https://www.gnu.org/software/gsl/doc/html/_images/math/4295789ac626fbd1804de8c52630d8583a4f2a7d.png).

- gsl_complex `gsl_complex_inverse`(gsl_complex *z*)

  This function returns the inverse, or reciprocal, of the complex number `z`, ![1/z = (x - i y)/(x^2 + y^2)](https://www.gnu.org/software/gsl/doc/html/_images/math/e98d18d5b399d39faf9a5b19b9833954e2cb9c07.png).

- gsl_complex `gsl_complex_negative`(gsl_complex *z*)

  This function returns the negative of the complex number `z`, ![-z = (-x) + i(-y)](https://www.gnu.org/software/gsl/doc/html/_images/math/6ddc62463caccd9adc7c68a0ffa3e5e4a61f2a5c.png).

## Elementary Complex Functions



- gsl_complex `gsl_complex_sqrt`(gsl_complex *z*)

  This function returns the square root of the complex number `z`, ![\sqrt z](https://www.gnu.org/software/gsl/doc/html/_images/math/1e93bf4c46b55634caae6f651ba3bda0466c398e.png). The branch cut is the negative real axis. The result always lies in the right half of the complex plane.

- gsl_complex `gsl_complex_sqrt_real`(double *x*)

  This function returns the complex square root of the real number `x`, where `x` may be negative.



- gsl_complex `gsl_complex_pow`(gsl_complex *z*, gsl_complex *a*)

  The function returns the complex number `z` raised to the complex power `a`, ![z^a](https://www.gnu.org/software/gsl/doc/html/_images/math/45e05d55d25801d91fafe6e6be32fa8e7425055d.png). This is computed as ![\exp(\log(z)*a)](https://www.gnu.org/software/gsl/doc/html/_images/math/6123e97b37743a646ce2d755a8158d95f115508e.png) using complex logarithms and complex exponentials.

- gsl_complex `gsl_complex_pow_real`(gsl_complex *z*, double *x*)

  This function returns the complex number `z` raised to the real power `x`, ![z^x](https://www.gnu.org/software/gsl/doc/html/_images/math/955322a64b1d256608972925a62ca5954fa543ed.png).

- gsl_complex `gsl_complex_exp`(gsl_complex *z*)

  This function returns the complex exponential of the complex number `z`, ![\exp(z)](https://www.gnu.org/software/gsl/doc/html/_images/math/699391713a6b0b09a7ee484bab1d3cbfdaaa48c3.png).



- gsl_complex `gsl_complex_log`(gsl_complex *z*)

  This function returns the complex natural logarithm (base ![e](https://www.gnu.org/software/gsl/doc/html/_images/math/9a7ff5d0469d9a20a6a7da1b2c1ab9b1014d53c9.png)) of the complex number `z`, ![\log(z)](https://www.gnu.org/software/gsl/doc/html/_images/math/d048492e6a2d239de32450eb6d3cc09da5e0f530.png). The branch cut is the negative real axis.

- gsl_complex `gsl_complex_log10`(gsl_complex *z*)

  This function returns the complex base-10 logarithm of the complex number `z`, ![\log_{10} (z)](https://www.gnu.org/software/gsl/doc/html/_images/math/bd7d7e717491b4fcb0c5a12902cffb394e7dff70.png).

- gsl_complex `gsl_complex_log_b`(gsl_complex *z*, gsl_complex *b*)

  This function returns the complex base-`b` logarithm of the complex number `z`, ![\log_b(z)](https://www.gnu.org/software/gsl/doc/html/_images/math/56c6abebbe666c8bdcf26a5f6ac598d09541800e.png). This quantity is computed as the ratio ![\log(z)/\log(b)](https://www.gnu.org/software/gsl/doc/html/_images/math/82f11fed5082089a3c68b383f46169711f5996a6.png).



## Complex Trigonometric Functions



- gsl_complex `gsl_complex_sin`(gsl_complex *z*)

  This function returns the complex sine of the complex number `z`, ![\sin(z) = (\exp(iz) - \exp(-iz))/(2i)](https://www.gnu.org/software/gsl/doc/html/_images/math/abb8f35c9587ffc69e231fc2ea37d1c6da0e6fa0.png).



- gsl_complex `gsl_complex_cos`(gsl_complex *z*)

  This function returns the complex cosine of the complex number `z`, ![\cos(z) = (\exp(iz) + \exp(-iz))/2](https://www.gnu.org/software/gsl/doc/html/_images/math/46675efc532a592fce3f5b4a9acad2c3d6c864c0.png).



- gsl_complex `gsl_complex_tan`(gsl_complex *z*)

  This function returns the complex tangent of the complex number `z`, ![\tan(z) = \sin(z)/\cos(z)](https://www.gnu.org/software/gsl/doc/html/_images/math/7b71b42aece73564e3d5dca1335ab315fd29d9c8.png).

- gsl_complex `gsl_complex_sec`(gsl_complex *z*)

  This function returns the complex secant of the complex number `z`, ![\sec(z) = 1/\cos(z)](https://www.gnu.org/software/gsl/doc/html/_images/math/353f32e38cf6de496bb3bbea1acba97c12be248f.png).

- gsl_complex `gsl_complex_csc`(gsl_complex *z*)

  This function returns the complex cosecant of the complex number `z`, ![\csc(z) = 1/\sin(z)](https://www.gnu.org/software/gsl/doc/html/_images/math/64d9c0ac2d15aadab9997399befd7c1c08dca7c9.png).

- gsl_complex `gsl_complex_cot`(gsl_complex *z*)

  This function returns the complex cotangent of the complex number `z`, ![\cot(z) = 1/\tan(z)](https://www.gnu.org/software/gsl/doc/html/_images/math/158cebb8a7cd89975c899bc204ca68850c368400.png).



## Inverse Complex Trigonometric Functions

- gsl_complex `gsl_complex_arcsin`(gsl_complex *z*)

  This function returns the complex arcsine of the complex number `z`, ![\arcsin(z)](https://www.gnu.org/software/gsl/doc/html/_images/math/864f7521cc80b742d5618ea2b3c9058448b7098a.png). The branch cuts are on the real axis, less than ![-1](https://www.gnu.org/software/gsl/doc/html/_images/math/e26aecc9561355dd0e6fbf2742d0ff0d6b58c864.png) and greater than ![1](https://www.gnu.org/software/gsl/doc/html/_images/math/585581da152cce359927e7b9eeadef5ca7c3a42d.png).

- gsl_complex `gsl_complex_arcsin_real`(double *z*)

  This function returns the complex arcsine of the real number `z`, ![\arcsin(z)](https://www.gnu.org/software/gsl/doc/html/_images/math/864f7521cc80b742d5618ea2b3c9058448b7098a.png). For ![z](https://www.gnu.org/software/gsl/doc/html/_images/math/4c9e578e25d40c5b0c5cad91a1e2f34e58b5be5e.png) between ![-1](https://www.gnu.org/software/gsl/doc/html/_images/math/e26aecc9561355dd0e6fbf2742d0ff0d6b58c864.png)and ![1](https://www.gnu.org/software/gsl/doc/html/_images/math/585581da152cce359927e7b9eeadef5ca7c3a42d.png), the function returns a real value in the range ![[-\pi/2,\pi/2]](https://www.gnu.org/software/gsl/doc/html/_images/math/66cc619e3f53bbfb9e46f2126bab6db1cda6a056.png). For ![z](https://www.gnu.org/software/gsl/doc/html/_images/math/4c9e578e25d40c5b0c5cad91a1e2f34e58b5be5e.png) less than ![-1](https://www.gnu.org/software/gsl/doc/html/_images/math/e26aecc9561355dd0e6fbf2742d0ff0d6b58c864.png) the result has a real part of ![-\pi/2](https://www.gnu.org/software/gsl/doc/html/_images/math/2fd1133dcbbbb76109c2d0aefa1cda7cb9add5e5.png) and a positive imaginary part. For ![z](https://www.gnu.org/software/gsl/doc/html/_images/math/4c9e578e25d40c5b0c5cad91a1e2f34e58b5be5e.png) greater than ![1](https://www.gnu.org/software/gsl/doc/html/_images/math/585581da152cce359927e7b9eeadef5ca7c3a42d.png) the result has a real part of ![\pi/2](https://www.gnu.org/software/gsl/doc/html/_images/math/dc93cceb6e7be84e1d3470cd32db7aa763298782.png) and a negative imaginary part.

- gsl_complex `gsl_complex_arccos`(gsl_complex *z*)

  This function returns the complex arccosine of the complex number `z`, ![\arccos(z)](https://www.gnu.org/software/gsl/doc/html/_images/math/63b67a54694bafb079e27ac32db8e0b48fe3d82d.png). The branch cuts are on the real axis, less than ![-1](https://www.gnu.org/software/gsl/doc/html/_images/math/e26aecc9561355dd0e6fbf2742d0ff0d6b58c864.png) and greater than ![1](https://www.gnu.org/software/gsl/doc/html/_images/math/585581da152cce359927e7b9eeadef5ca7c3a42d.png).

- gsl_complex `gsl_complex_arccos_real`(double *z*)

  This function returns the complex arccosine of the real number `z`, ![\arccos(z)](https://www.gnu.org/software/gsl/doc/html/_images/math/63b67a54694bafb079e27ac32db8e0b48fe3d82d.png). For ![z](https://www.gnu.org/software/gsl/doc/html/_images/math/4c9e578e25d40c5b0c5cad91a1e2f34e58b5be5e.png) between ![-1](https://www.gnu.org/software/gsl/doc/html/_images/math/e26aecc9561355dd0e6fbf2742d0ff0d6b58c864.png) and ![1](https://www.gnu.org/software/gsl/doc/html/_images/math/585581da152cce359927e7b9eeadef5ca7c3a42d.png), the function returns a real value in the range ![[0,\pi]](https://www.gnu.org/software/gsl/doc/html/_images/math/65fbac805f66f8d00f94040e4e3b2f84d28705a5.png). For ![z](https://www.gnu.org/software/gsl/doc/html/_images/math/4c9e578e25d40c5b0c5cad91a1e2f34e58b5be5e.png) less than ![-1](https://www.gnu.org/software/gsl/doc/html/_images/math/e26aecc9561355dd0e6fbf2742d0ff0d6b58c864.png) the result has a real part of ![\pi](https://www.gnu.org/software/gsl/doc/html/_images/math/e619909e385649c40e3610ab5d7f9f8643e8b14b.png) and a negative imaginary part. For ![z](https://www.gnu.org/software/gsl/doc/html/_images/math/4c9e578e25d40c5b0c5cad91a1e2f34e58b5be5e.png) greater than ![1](https://www.gnu.org/software/gsl/doc/html/_images/math/585581da152cce359927e7b9eeadef5ca7c3a42d.png) the result is purely imaginary and positive.

- gsl_complex `gsl_complex_arctan`(gsl_complex *z*)

  This function returns the complex arctangent of the complex number `z`, ![\arctan(z)](https://www.gnu.org/software/gsl/doc/html/_images/math/f74c979e3a9511564c4b1b145ecfeee6b1e5f52c.png). The branch cuts are on the imaginary axis, below ![-i](https://www.gnu.org/software/gsl/doc/html/_images/math/6f549eec571d5bac35a3be924a58077839f0cb74.png) and above ![i](https://www.gnu.org/software/gsl/doc/html/_images/math/1cc632900aae8b2837a1f383619c0ad753be7d29.png).

- gsl_complex `gsl_complex_arcsec`(gsl_complex *z*)

  This function returns the complex arcsecant of the complex number `z`, ![\arcsec(z) = \arccos(1/z)](https://www.gnu.org/software/gsl/doc/html/_images/math/02d59402e7dba5833c48cb65ca3661c0be58559b.png).

- gsl_complex `gsl_complex_arcsec_real`(double *z*)

  This function returns the complex arcsecant of the real number `z`, ![\arcsec(z) = \arccos(1/z)](https://www.gnu.org/software/gsl/doc/html/_images/math/02d59402e7dba5833c48cb65ca3661c0be58559b.png).

- gsl_complex `gsl_complex_arccsc`(gsl_complex *z*)

  This function returns the complex arccosecant of the complex number `z`, ![\arccsc(z) = \arcsin(1/z)](https://www.gnu.org/software/gsl/doc/html/_images/math/21f4838c8cfbb1536d751ddc6e9013aaee394037.png).

- gsl_complex `gsl_complex_arccsc_real`(double *z*)

  This function returns the complex arccosecant of the real number `z`, ![\arccsc(z) = \arcsin(1/z)](https://www.gnu.org/software/gsl/doc/html/_images/math/21f4838c8cfbb1536d751ddc6e9013aaee394037.png).

- gsl_complex `gsl_complex_arccot`(gsl_complex *z*)

  This function returns the complex arccotangent of the complex number `z`, ![\arccot(z) = \arctan(1/z)](https://www.gnu.org/software/gsl/doc/html/_images/math/cf6a72e57f0004c0a136855af094a3a9098ca001.png).



## Complex Hyperbolic Functions

- gsl_complex `gsl_complex_sinh`(gsl_complex *z*)

  This function returns the complex hyperbolic sine of the complex number `z`, ![\sinh(z) = (\exp(z) - \exp(-z))/2](https://www.gnu.org/software/gsl/doc/html/_images/math/1f7237d06473ac31ebad5affeac17d2e44ec4298.png).

- gsl_complex `gsl_complex_cosh`(gsl_complex *z*)

  This function returns the complex hyperbolic cosine of the complex number `z`, ![\cosh(z) = (\exp(z) + \exp(-z))/2](https://www.gnu.org/software/gsl/doc/html/_images/math/ccbd6cbbb77122003fccb039ee145aa991cc3c49.png).

- gsl_complex `gsl_complex_tanh`(gsl_complex *z*)

  This function returns the complex hyperbolic tangent of the complex number `z`, ![\tanh(z) = \sinh(z)/\cosh(z)](https://www.gnu.org/software/gsl/doc/html/_images/math/732c1b637078035c3e5bec897040b79a09b81a12.png).

- gsl_complex `gsl_complex_sech`(gsl_complex *z*)

  This function returns the complex hyperbolic secant of the complex number `z`, ![\sech(z) = 1/\cosh(z)](https://www.gnu.org/software/gsl/doc/html/_images/math/192b5d426b8d653d2ef99315cd8b6165c3b49b71.png).

- gsl_complex `gsl_complex_csch`(gsl_complex *z*)

  This function returns the complex hyperbolic cosecant of the complex number `z`, ![\csch(z) = 1/\sinh(z)](https://www.gnu.org/software/gsl/doc/html/_images/math/b56274ee4b91999f28b707579de457abf61294e5.png).

- gsl_complex `gsl_complex_coth`(gsl_complex *z*)

  This function returns the complex hyperbolic cotangent of the complex number `z`, ![\coth(z) = 1/\tanh(z)](https://www.gnu.org/software/gsl/doc/html/_images/math/8101ad4bd2e6de7721a4b633cf339959d6b21838.png).



## Inverse Complex Hyperbolic Functions

- gsl_complex `gsl_complex_arcsinh`(gsl_complex *z*)

  This function returns the complex hyperbolic arcsine of the complex number `z`, ![\arcsinh(z)](https://www.gnu.org/software/gsl/doc/html/_images/math/8f5540e84fd888b6905753af8079fd7610228cf1.png). The branch cuts are on the imaginary axis, below ![-i](https://www.gnu.org/software/gsl/doc/html/_images/math/6f549eec571d5bac35a3be924a58077839f0cb74.png) and above ![i](https://www.gnu.org/software/gsl/doc/html/_images/math/1cc632900aae8b2837a1f383619c0ad753be7d29.png).

- gsl_complex `gsl_complex_arccosh`(gsl_complex *z*)

  This function returns the complex hyperbolic arccosine of the complex number `z`, ![\arccosh(z)](https://www.gnu.org/software/gsl/doc/html/_images/math/7c26272972c0a7c6f4a66153407443b8bdd8f802.png). The branch cut is on the real axis, less than ![1](https://www.gnu.org/software/gsl/doc/html/_images/math/585581da152cce359927e7b9eeadef5ca7c3a42d.png). Note that in this case we use the negative square root in formula 4.6.21 of Abramowitz & Stegun giving ![\arccosh(z)=\log(z-\sqrt{z^2-1})](https://www.gnu.org/software/gsl/doc/html/_images/math/7bd57f3c72b59c948914117370a594397be4fc52.png).

- gsl_complex `gsl_complex_arccosh_real`(double *z*)

  This function returns the complex hyperbolic arccosine of the real number `z`, ![\arccosh(z)](https://www.gnu.org/software/gsl/doc/html/_images/math/7c26272972c0a7c6f4a66153407443b8bdd8f802.png).

- gsl_complex `gsl_complex_arctanh`(gsl_complex *z*)

  This function returns the complex hyperbolic arctangent of the complex number `z`, ![\arctanh(z)](https://www.gnu.org/software/gsl/doc/html/_images/math/ebe0cc77e6b81598fe93bf53030a2f319f58dc97.png). The branch cuts are on the real axis, less than ![-1](https://www.gnu.org/software/gsl/doc/html/_images/math/e26aecc9561355dd0e6fbf2742d0ff0d6b58c864.png) and greater than ![1](https://www.gnu.org/software/gsl/doc/html/_images/math/585581da152cce359927e7b9eeadef5ca7c3a42d.png).

- gsl_complex `gsl_complex_arctanh_real`(double *z*)

  This function returns the complex hyperbolic arctangent of the real number `z`, ![\arctanh(z)](https://www.gnu.org/software/gsl/doc/html/_images/math/ebe0cc77e6b81598fe93bf53030a2f319f58dc97.png).

- gsl_complex `gsl_complex_arcsech`(gsl_complex *z*)

  This function returns the complex hyperbolic arcsecant of the complex number `z`, ![\arcsech(z) = \arccosh(1/z)](https://www.gnu.org/software/gsl/doc/html/_images/math/57f15c6b7cc9ca3e937ba210b94083588e49582d.png).

- gsl_complex `gsl_complex_arccsch`(gsl_complex *z*)

  This function returns the complex hyperbolic arccosecant of the complex number `z`, ![\arccsch(z) = \arcsinh(1/z)](https://www.gnu.org/software/gsl/doc/html/_images/math/bab8831a88374ae10f7b4ac3ff0e4483235ce56d.png).

- gsl_complex `gsl_complex_arccoth`(gsl_complex *z*)

  This function returns the complex hyperbolic arccotangent of the complex number `z`, ![\arccoth(z) = \arctanh(1/z)](https://www.gnu.org/software/gsl/doc/html/_images/math/5d7dabcf49975ef8a12e7878b1590d588c51de4e.png).

## References and Further Reading

The implementations of the elementary and trigonometric functions are based on the following papers,

- T. E. Hull, Thomas F. Fairgrieve, Ping Tak Peter Tang, “Implementing Complex Elementary Functions Using Exception Handling”, ACM Transactions on Mathematical Software, Volume 20 (1994), pp 215–244, Corrigenda, p553
- T. E. Hull, Thomas F. Fairgrieve, Ping Tak Peter Tang, “Implementing the complex arcsin and arccosine functions using exception handling”, ACM Transactions on Mathematical Software, Volume 23 (1997) pp 299–335

The general formulas and details of branch cuts can be found in the following books,

- Abramowitz and Stegun, Handbook of Mathematical Functions, “Circular Functions in Terms of Real and Imaginary Parts”, Formulas 4.3.55–58, “Inverse Circular Functions in Terms of Real and Imaginary Parts”, Formulas 4.4.37–39, “Hyperbolic Functions in Terms of Real and Imaginary Parts”, Formulas 4.5.49–52, “Inverse Hyperbolic Functions—relation to Inverse Circular Functions”, Formulas 4.6.14–19.
- Dave Gillespie, Calc Manual, Free Software Foundation, ISBN 1-882114-18-3

### Footnotes

[[1]](https://www.gnu.org/software/gsl/doc/html/complex.html#id1)Note that the first edition uses different definitions.