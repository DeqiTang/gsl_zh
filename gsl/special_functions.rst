

.. index:: special functions

********
特殊函数
********

This chapter describes the GSL special function library.  The library
includes routines for calculating the values of Airy functions, Bessel
functions, Clausen functions, Coulomb wave functions, Coupling
coefficients, the Dawson function, Debye functions, Dilogarithms,
Elliptic integrals, Jacobi elliptic functions, Error functions,
Exponential integrals, Fermi-Dirac functions, Gamma functions,
Gegenbauer functions, Hermite polynomials and functions, Hypergeometric functions, Laguerre functions,
Legendre functions and Spherical Harmonics, the Psi (Digamma) Function,
Synchrotron functions, Transport functions, Trigonometric functions and
Zeta functions.  Each routine also computes an estimate of the numerical
error in the calculated value of the function.

本章描述GSL特殊函数库。本库包含计算艾里函数，贝塞尔函数，克劳森函数，库伦波函数，耦合系数，道森函数，德拜函数，二重对数，椭圆积分，雅可比椭圆函数，误差函数，指数积分，费米狄拉克函数，伽马函数，盖根堡函数，厄米特多项式和函数，超几何函数，拉盖尔函数，勒让德函数和球谐函数，Psi(双伽马)函数，同步加速函数，输运函数，三角函数和ζ函数的值的程序。每个程序也计算函数所计算的值的数值误差的一个估值。


The functions in this chapter are declared in individual header files,
such as :file:`gsl_sf_airy.h`, :file:`gsl_sf_bessel.h`, etc.  The complete
set of header files can be included using the file :file:`gsl_sf.h`.

本章描述的函数被声明在单独的头文件, 比如\ :file:`gsl_sf_airy.h`\, \ :file:`gsl_sf_bessel.h`\等。头文件的完全集能通过文件\ :file:`gsl_sf.h`\进行包含。


用法
=====

The special functions are available in two calling conventions, a
*natural form* which returns the numerical value of the function and
an *error-handling form* which returns an error code.  The two types
of function provide alternative ways of accessing the same underlying
code.

特殊函数可通过两个调用约定进行获取，一个是\ *自然形式* \其返回函数的数值，另一个是\ *误差处理形式*\其返回一个误差代码。两个类型的函数提供了访问同样的底层代码的可选方式。

The *natural form* returns only the value of the function and can be
used directly in mathematical expressions.  For example, the following
function call will compute the value of the Bessel function :math:`J_0(x)`::

    double y = gsl_sf_bessel_J0 (x);

\ *自然形式* \仅返回函数的值而可以被直接用在数学表达式中。比如，下面的函数
调用将会计算Bessel函数\ :math:`J_0(x)`\::

    double y = gsl_sf_bessel_J0(x);

There is no way to access an error code or to estimate the error using
this method.  To allow access to this information the alternative
error-handling form stores the value and error in a modifiable argument::

    gsl_sf_result result;
    int status = gsl_sf_bessel_J0_e (x, &result);

使用这个方法是不可以获取误差代码或者估计误差值。要允许访问该信息，可选
的误差处理方法是将结果值和误差存在一个可修改参数中::
    
    gsl_sf_result result;
    int status = gsl_sf_bessel_J0_e (x, &result):

The error-handling functions have the suffix :code:`_e`. The returned
status value indicates error conditions such as overflow, underflow or
loss of precision.  If there are no errors the error-handling functions
return :code:`GSL_SUCCESS`.

gsl_sf_result 结构
========================

The error handling form of the special functions always calculate an
error estimate along with the value of the result.  Therefore,
structures are provided for amalgamating a value and error estimate.
These structures are declared in the header file :file:`gsl_sf_result.h`.

特俗函数的异常处理形式总是在计算结果时同时计算出误差估计值。因此返回结构是
为了将结果值与误差估计合并起来。
这些结构被声明在了头文件\ :file:`gsl_sf_result.h` \中。

The following struct contains value and error fields.

.. type:: gsl_sf_result

   ::

     typedef struct
     {
       double val;
       double err;
     } gsl_sf_result;

   The field :data:`val` contains the value and the field :data:`err` contains
   an estimate of the absolute error in the value.

In some cases, an overflow or underflow can be detected and handled by a
function.  In this case, it may be possible to return a scaling exponent
as well as an error/value pair in order to save the result from
exceeding the dynamic range of the built-in types.  The
following struct contains value and error fields as well
as an exponent field such that the actual result is obtained as
:code:`result * 10^(e10)`.

.. type:: gsl_sf_result_e10

   ::

     typedef struct
     {
       double val;
       double err;
       int    e10;
     } gsl_sf_result_e10;

Modes
=====

The goal of the library is to achieve double precision accuracy wherever
possible.  However the cost of evaluating some special functions to
double precision can be significant, particularly where very high order
terms are required.  In these cases a :code:`mode` argument, of type
:type:`gsl_mode_t` allows the
accuracy of the function to be reduced in order to improve performance.
The following precision levels are available for the mode argument,

.. type:: gsl_mode_t

   .. macro:: GSL_PREC_DOUBLE

      Double-precision, a relative accuracy of approximately :math:`2 * 10^{-16}`.

   .. macro:: GSL_PREC_SINGLE

      Single-precision, a relative accuracy of approximately :math:`10^{-7}`.

   .. macro:: GSL_PREC_APPROX

      Approximate values, a relative accuracy of approximately :math:`5 * 10^{-4}`.

The approximate mode provides the fastest evaluation at the lowest
accuracy.

Airy Functions and Derivatives
==============================
.. include:: specfunc-airy.rst

Bessel Functions
================
.. include:: specfunc-bessel.rst

Clausen Functions
=================
.. include:: specfunc-clausen.rst

Coulomb Functions
=================
.. include:: specfunc-coulomb.rst

Coupling Coefficients
=====================
.. include:: specfunc-coupling.rst

Dawson Function
===============
.. include:: specfunc-dawson.rst

Debye Functions
===============
.. include:: specfunc-debye.rst

.. _dilog-function:

Dilogarithm
===========
.. include:: specfunc-dilog.rst

Elementary Operations
=====================
.. include:: specfunc-elementary.rst

Elliptic Integrals
==================
.. include:: specfunc-ellint.rst

Elliptic Functions (Jacobi)
===========================
.. include:: specfunc-elljac.rst

Error Functions
===============
.. include:: specfunc-erf.rst

Exponential Functions
=====================
.. include:: specfunc-exp.rst

Exponential Integrals
=====================
.. include:: specfunc-expint.rst

Fermi-Dirac Function
====================
.. include:: specfunc-fermi-dirac.rst

Gamma and Beta Functions
========================
.. include:: specfunc-gamma.rst

Gegenbauer Functions
====================
.. include:: specfunc-gegenbauer.rst

Hermite Polynomials and Functions
=================================
.. include:: specfunc-hermite.rst

Hypergeometric Functions
========================
.. include:: specfunc-hyperg.rst

.. _laguerre-functions:

Laguerre Functions
==================
.. include:: specfunc-laguerre.rst

Lambert W Functions
===================
.. include:: specfunc-lambert.rst

Legendre Functions and Spherical Harmonics
==========================================
.. include:: specfunc-legendre.rst

Logarithm and Related Functions
===============================
.. include:: specfunc-log.rst

Mathieu Functions
=================
.. include:: specfunc-mathieu.rst

Power Function
==============
.. include:: specfunc-pow-int.rst

Psi (Digamma) Function
======================
.. include:: specfunc-psi.rst

Synchrotron Functions
=====================
.. include:: specfunc-synchrotron.rst

Transport Functions
===================
.. include:: specfunc-transport.rst

Trigonometric Functions
=======================
.. include:: specfunc-trig.rst

Zeta Functions
==============
.. include:: specfunc-zeta.rst

Examples
========

The following example demonstrates the use of the error handling form of
the special functions, in this case to compute the Bessel function
:math:`J_0(5.0)`,

.. include:: examples/specfun_e.c
   :code:

Here are the results of running the program,

.. include:: examples/specfun_e.txt
   :code:

The next program computes the same quantity using the natural form of
the function. In this case the error term :data:`result.err` and return
status are not accessible.

.. include:: examples/specfun.c
   :code:

The results of the function are the same,

.. include:: examples/specfun.txt
   :code:

References and Further Reading
==============================

The library follows the conventions of the following book where possible,

* Handbook of Mathematical Functions, edited by Abramowitz & Stegun,
  Dover,  ISBN 0486612724.

The following papers contain information on the algorithms used 
to compute the special functions,

.. index:: MISCFUN

* Allan J. MacLeod, MISCFUN: A software package to compute uncommon
  special functions.  ACM Trans. Math. Soft., vol.: 22,
  1996, 288--301

* Bunck, B. F., A fast algorithm for evaluation of normalized Hermite
  functions, BIT Numer. Math, 49: 281-295, 2009.

* G.N. Watson, A Treatise on the Theory of Bessel Functions,
  2nd Edition (Cambridge University Press, 1944).

* G. Nemeth, Mathematical Approximations of Special Functions,
  Nova Science Publishers, ISBN 1-56072-052-2

* B.C. Carlson, Special Functions of Applied Mathematics (1977)

* N. M. Temme, Special Functions: An Introduction to the Classical
  Functions of Mathematical Physics (1996), ISBN 978-0471113133.

* W.J. Thompson, Atlas for Computing Mathematical Functions, John Wiley & Sons,
  New York (1997).

* Y.Y. Luke, Algorithms for the Computation of Mathematical Functions, Academic
  Press, New York (1977).

* S. A. Holmes and W. E. Featherstone, A unified approach to the Clenshaw
  summation and the recursive computation of very high degree and order
  normalised associated Legendre functions, Journal of Geodesy, 76,
  pg. 279-299, 2002.
