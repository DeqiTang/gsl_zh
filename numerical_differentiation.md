# Numerical Differentiation数值微分

The functions described in this chapter compute numerical derivatives by finite differencing. An adaptive algorithm is used to find the best choice of finite difference and to estimate the error in the derivative. These functions are declared in the header file `gsl_deriv.h`.

本章描述的函数通过有限差分来计算数值导数。一个自适应算法被使用来寻找有限差分的最好选择并估计导数的误差。这些函数被声明在头文件`gsl_derive.h`中。

## Functions函数

- int `gsl_deriv_central`(const [gsl_function](https://www.gnu.org/software/gsl/doc/html/roots.html#c.gsl_function) * *f*, double *x*, double *h*, double * *result*, double * *abserr*)

  This function computes the numerical derivative of the function `f` at the point `x` using an adaptive central difference algorithm with a step-size of `h`. The derivative is returned in `result` and an estimate of its absolute error is returned in `abserr`.The initial value of `h` is used to estimate an optimal step-size, based on the scaling of the truncation error and round-off error in the derivative calculation. The derivative is computed using a 5-point rule for equally spaced abscissae at ![x - h](https://www.gnu.org/software/gsl/doc/html/_images/math/0983613ed90fec490adad5526090d1a3bb6560bd.png), ![x - h/2](https://www.gnu.org/software/gsl/doc/html/_images/math/63589fdd754ec9ada8b1774bfbfe05d8afc0c3ac.png), ![x](https://www.gnu.org/software/gsl/doc/html/_images/math/1c5715c8edd0bbc330baa0f0494f1f7951619754.png), ![x + h/2](https://www.gnu.org/software/gsl/doc/html/_images/math/f8e2566f1296887e95dacf80737f8c5e943f1e73.png), ![x+h](https://www.gnu.org/software/gsl/doc/html/_images/math/6418f9ab3f15b84c1b188e96dd3d26658af9a5e6.png), with an error estimate taken from the difference between the 5-point rule and the corresponding 3-point rule ![x-h](https://www.gnu.org/software/gsl/doc/html/_images/math/e0aa92a0b978b557ea2b64f1d5c18d3e679d3c0b.png), ![x](https://www.gnu.org/software/gsl/doc/html/_images/math/1c5715c8edd0bbc330baa0f0494f1f7951619754.png), ![x+h](https://www.gnu.org/software/gsl/doc/html/_images/math/6418f9ab3f15b84c1b188e96dd3d26658af9a5e6.png). Note that the value of the function at ![x](https://www.gnu.org/software/gsl/doc/html/_images/math/1c5715c8edd0bbc330baa0f0494f1f7951619754.png) does not contribute to the derivative calculation, so only 4-points are actually used.

  这个函数使用一个步长为`h`的自适应中心差分算法来计算函数`f`在点`x`的数值导数。导数被返回到`result`中，并且其绝对误差的一个估值被返回到`abserr`中。`h`的初值被用来估算一个最佳的步长，基于导数计算中截断误差和舍入误差的扩展。导数是通过使用横坐标上等距的$$x-h, x - h/2, x, x+h/2, x+h$$五点规则进行计算，同时误差估值从5点之间的差异以及相应的三点规则$$x-h, x , x+h$$中求得。注意函数在$$x$$上的值对导数计算没有贡献，因此仅有4个点被实际计算。

- int `gsl_deriv_forward`(const [gsl_function](https://www.gnu.org/software/gsl/doc/html/roots.html#c.gsl_function) * *f*, double *x*, double *h*, double * *result*, double * *abserr*)

  This function computes the numerical derivative of the function `f` at the point `x` using an adaptive forward difference algorithm with a step-size of `h`. The function is evaluated only at points greater than `x`, and never at `x` itself. The derivative is returned in `result` and an estimate of its absolute error is returned in `abserr`. This function should be used if ![f(x)](https://www.gnu.org/software/gsl/doc/html/_images/math/c7deb6ce5befe14135ebd23fa69801be7a796b15.png) has a discontinuity at `x`, or is undefined for values less than `x`.The initial value of `h` is used to estimate an optimal step-size, based on the scaling of the truncation error and round-off error in the derivative calculation. The derivative at ![x](https://www.gnu.org/software/gsl/doc/html/_images/math/1c5715c8edd0bbc330baa0f0494f1f7951619754.png) is computed using an “open” 4-point rule for equally spaced abscissae at ![x+h/4](https://www.gnu.org/software/gsl/doc/html/_images/math/d19ca822e823ac7780401fdede8b6ee0c27c0ef7.png), ![x + h/2](https://www.gnu.org/software/gsl/doc/html/_images/math/f8e2566f1296887e95dacf80737f8c5e943f1e73.png), ![x + 3h/4](https://www.gnu.org/software/gsl/doc/html/_images/math/ee0e234164a9ed0e1a03f343048827301698b529.png), ![x+h](https://www.gnu.org/software/gsl/doc/html/_images/math/6418f9ab3f15b84c1b188e96dd3d26658af9a5e6.png), with an error estimate taken from the difference between the 4-point rule and the corresponding 2-point rule ![x+h/2](https://www.gnu.org/software/gsl/doc/html/_images/math/2ef59434d2ae809dcf070d956f05b28a9485969e.png), ![x+h](https://www.gnu.org/software/gsl/doc/html/_images/math/6418f9ab3f15b84c1b188e96dd3d26658af9a5e6.png).

  这个函数使用一个步长为`h`的自适应前向差分算法来计算函数`f`在`x`上的数值导数。函数仅在大于`x`的电商被求值，并且永远不会在`x`点本身上。导数被返回在`result`中，且其绝对误差的一个估值被返回到`abserr`中。这个函数应当在$$f(x)$$在点`x`不连续或者在小于`x`时未被定义时被采用。初值`h`被用来估算一个最佳的步长，基于导数计算中截断误差和舍入误差的扩展。$$x$$点处的导数是通过一个横坐标上的等间距点$$x+h/4, x+h/2,x+3h/4,x+h$$的"open"4点规则进行计算，同事误差估值从4点的差异以及相应的2点规则$$x+h/2, x+h$$中求得。

- int `gsl_deriv_backward`(const [gsl_function](https://www.gnu.org/software/gsl/doc/html/roots.html#c.gsl_function) * *f*, double *x*, double *h*, double * *result*, double * *abserr*)

  This function computes the numerical derivative of the function `f` at the point `x` using an adaptive backward difference algorithm with a step-size of `h`. The function is evaluated only at points less than `x`, and never at `x` itself. The derivative is returned in `result` and an estimate of its absolute error is returned in `abserr`. This function should be used if ![f(x)](https://www.gnu.org/software/gsl/doc/html/_images/math/c7deb6ce5befe14135ebd23fa69801be7a796b15.png) has a discontinuity at `x`, or is undefined for values greater than `x`.This function is equivalent to calling [`gsl_deriv_forward()`](https://www.gnu.org/software/gsl/doc/html/diff.html#c.gsl_deriv_forward) with a negative step-size.

  这个函数使用一个步长为`h`的自适应后向差分算法来计算函数`f`在点`x`上的数值导数。函数仅在大于`x`的点上进行求值，并且永远不会在`x`点本身上。导数被返回到`result`中，且其绝对误差的一个估值被返回到`abserr`中。这个函数应当在$$f(x)$$在`x`不连续或者在大于`x`未定义时被使用。这个函数等价于调用通过一个负步长调用[`gsl_deriv_forward()`](https://www.gnu.org/software/gsl/doc/html/diff.html#c.gsl_deriv_forward) 。

## 示例

The following code estimates the derivative of the function ![f(x) = x^{3/2}](https://www.gnu.org/software/gsl/doc/html/_images/math/594cf07b54e7ca611902d8c02a3f06c46c2723f8.png) at ![x = 2](https://www.gnu.org/software/gsl/doc/html/_images/math/d12b4d167765d1113ff1ed159775f8317c4e8976.png) and at ![x = 0](https://www.gnu.org/software/gsl/doc/html/_images/math/d874e492e19ddb8e484e79e1f4b417bebaa83731.png). The function ![f(x)](https://www.gnu.org/software/gsl/doc/html/_images/math/c7deb6ce5befe14135ebd23fa69801be7a796b15.png) is undefined for ![x < 0](https://www.gnu.org/software/gsl/doc/html/_images/math/93f2a3ce4b90fb07a7dc88c1bc60a6a43c93f45b.png) so the derivative at ![x=0](https://www.gnu.org/software/gsl/doc/html/_images/math/df35f8ce9167aff62ff38e2d458feeaec4475441.png) is computed using [`gsl_deriv_forward()`](https://www.gnu.org/software/gsl/doc/html/diff.html#c.gsl_deriv_forward).

下面的代码估计了函数$$f(x)=x^{3/2}$$在$$x=2$$和$$x=0$$处的导数。函数$$f(x)$$在$$x<0$$时是未定义的，因此在$$x=0$$处的导数是通过[`gsl_deriv_forward()`](https://www.gnu.org/software/gsl/doc/html/diff.html#c.gsl_deriv_forward)进行计算。

```
#include <stdio.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_deriv.h>

double f (double x, void * params)
{
  (void)(params); /* avoid unused parameter warning */
  return pow (x, 1.5);
}

int
main (void)
{
  gsl_function F;
  double result, abserr;

  F.function = &f;
  F.params = 0;

  printf ("f(x) = x^(3/2)\n");

  gsl_deriv_central (&F, 2.0, 1e-8, &result, &abserr);
  printf ("x = 2.0\n");
  printf ("f'(x) = %.10f +/- %.10f\n", result, abserr);
  printf ("exact = %.10f\n\n", 1.5 * sqrt(2.0));

  gsl_deriv_forward (&F, 0.0, 1e-8, &result, &abserr);
  printf ("x = 0.0\n");
  printf ("f'(x) = %.10f +/- %.10f\n", result, abserr);
  printf ("exact = %.10f\n", 0.0);

  return 0;
}
```

Here is the output of the program,

下面是程序的输出，

```
f(x) = x^(3/2)
x = 2.0
f'(x) = 2.1213203120 +/- 0.0000005006
exact = 2.1213203436

x = 0.0
f'(x) = 0.0000000160 +/- 0.0000000339
exact = 0.0000000000
```

## References and Further Reading参考和进一步阅读

The algorithms used by these functions are described in the following sources:

这些函数使用到的算法在下面资源中被描述:

- Abramowitz and Stegun, *Handbook of Mathematical Functions*, Section 25.3.4, and Table 25.5 (Coefficients for Differentiation).
- S.D. Conte and Carl de Boor, *Elementary Numerical Analysis: An Algorithmic Approach*, McGraw-Hill, 1972.