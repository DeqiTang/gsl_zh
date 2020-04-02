# 小波变换

This chapter describes functions for performing Discrete Wavelet Transforms (DWTs). The library includes wavelets for real data in both one and two dimensions. The wavelet functions are declared in the header files `gsl_wavelet.h` and `gsl_wavelet2d.h`.

本章描述用于执行离散小波变换(DWTs)。本库包含包括一维和二维的实数的小波变换。小波函数被声明在头文件`gsl_wavelet.h`和`gsl_wavelet2d.h`中。

## 定义

The continuous wavelet transform and its inverse are defined by the relations,

![w(s, \tau) = \int_{-\infty}^\infty f(t) * \psi^*_{s,\tau}(t) dt](https://www.gnu.org/software/gsl/doc/html/_images/math/e05e828b1a786e212606c0f554e6dc39f3acb707.png)

and,

![f(t) = \int_0^\infty ds \int_{-\infty}^\infty w(s, \tau) * \psi_{s,\tau}(t) d\tau](https://www.gnu.org/software/gsl/doc/html/_images/math/bab436d0d1bf370fafcf04060f720bda0a27790a.png)

where the basis functions ![\psi_{s,\tau}](https://www.gnu.org/software/gsl/doc/html/_images/math/9ba110ff1af8da6586a7607e797865c44fc31d22.png) are obtained by scaling and translation from a single function, referred to as the {mother wavelet*.

连续小波变换和它的逆变换被以下关系定义，

$\omega (s,\tau) = \int_{-\infty}^{\infty}f(t)*\psi_{s,\tau}^*(t)dt$

和，

$f ( t ) = \int _ { 0 } ^ { \infty } d s \int _ { - \infty } ^ { \infty } w ( s , \tau ) * \psi _ { s , \tau } ( t ) d \tau$

其中基函数$$\psi_{s,\tau}$$是通过一个被称为小波母函数独立函数扩展和平移得到的。

The discrete version of the wavelet transform acts on equally-spaced samples, with fixed scaling and translation steps (![s](https://www.gnu.org/software/gsl/doc/html/_images/math/dd69d21b40087e08fc6827d2f3481f57462046e7.png), ![\tau](https://www.gnu.org/software/gsl/doc/html/_images/math/de9fbd420c3bcc6304590e768c308e743336d58d.png)). The frequency and time axes are sampled *dyadically* on scales of ![2^j](https://www.gnu.org/software/gsl/doc/html/_images/math/daabaa9bf289cf3da57f9cec49388d9fc78b73d6.png)through a level parameter ![j](https://www.gnu.org/software/gsl/doc/html/_images/math/7bc2ad053edeecad39b867d60495ac2caa4555dc.png).

The resulting family of functions ![\{\psi_{j,n}\}](https://www.gnu.org/software/gsl/doc/html/_images/math/44d42479b04466104aa0cfb964fad718532cebe7.png) constitutes an orthonormal basis for square-integrable signals. The discrete wavelet transform is an ![O(N)](https://www.gnu.org/software/gsl/doc/html/_images/math/cf9a0d26df85a5810610fdb1cd0e58492d00c419.png) algorithm, and is also referred to as the *fast wavelet transform*.

小波变换的离散版本作用在具有扩展和平移步$$(s,\tau)$$等间距的样品上。频率和时间轴通过二分体的方式在$$2^j$$尺度上通过一个层级参数$$j$$进行采样。

## 初始化

- `gsl_wavelet`

  This structure contains the filter coefficients defining the wavelet and any associated offset parameters.

- [gsl_wavelet](https://www.gnu.org/software/gsl/doc/html/dwt.html#c.gsl_wavelet) * `gsl_wavelet_alloc`(const [gsl_wavelet_type](https://www.gnu.org/software/gsl/doc/html/dwt.html#c.gsl_wavelet_type) * *T*, size_t *k*)

  This function allocates and initializes a wavelet object of type `T`. The parameter `k` selects the specific member of the wavelet family. A null pointer is returned if insufficient memory is available or if a unsupported member is selected.

- `gsl_wavelet`

  这个结构包含定义了小波和任何相关的补偿参数的过滤系数。

- [gsl_wavelet](https://www.gnu.org/software/gsl/doc/html/dwt.html#c.gsl_wavelet) * `gsl_wavelet_alloc`(const [gsl_wavelet_type](https://www.gnu.org/software/gsl/doc/html/dwt.html#c.gsl_wavelet_type) * *T*, size_t *k*)

  这个函数分配和初始化一个类型为`T`的小波对象。参数`k`选择小波家族的特定的成员。如果没有足够的内存可以获取或者不支持的成员被选择到就会返回一个空指针。

The following wavelet types are implemented:

- `gsl_wavelet_type` 

  `gsl_wavelet_daubechies` `gsl_wavelet_daubechies_centered`

  This is the Daubechies wavelet family of maximum phase with ![k/2](https://www.gnu.org/software/gsl/doc/html/_images/math/ec622d5416b57dfa17fc7b544fd223fb6f2c5632.png) vanishing moments. The implemented wavelets are ![k=4, 6, \dots, 20](https://www.gnu.org/software/gsl/doc/html/_images/math/0346e470a1962a3147689ed42ddd7a02a3db06ef.png), with `k` even.

  `gsl_wavelet_haar` `gsl_wavelet_haar_centered`

  This is the Haar wavelet. The only valid choice of ![k](https://www.gnu.org/software/gsl/doc/html/_images/math/90645b5c2abc13d7de104fd87abeb2f19406e5f8.png) for the Haar wavelet is ![k=2](https://www.gnu.org/software/gsl/doc/html/_images/math/07e922dd1b33506fc566e1d2d19ea41c6ae5c69b.png).

  `gsl_wavelet_bspline` `gsl_wavelet_bspline_centered`

  This is the biorthogonal B-spline wavelet family of order ![(i,j)](https://www.gnu.org/software/gsl/doc/html/_images/math/4b982cbeb3e7a1b695da41fa331848138a69542a.png). The implemented values of ![k = 100*i + j](https://www.gnu.org/software/gsl/doc/html/_images/math/2dddbeb8b272e20d88629b2a381e92931ab181aa.png) are 103, 105, 202, 204, 206, 208, 301, 303, 305 307, 309.

以下小波被实现了:

* `gsl_wavelet_type`

  `gsl_wavelet_daubechies` `gsl_wavelet_daubechies_centered`

  这是最大相的具有$$k/2$$的消失矩的Daubechies小波家族。被实现的小波是$$k = 4,6,...,20$$，$$k$$为偶数。

  `gsl_wavelet_haar` `gsl_wavelet_haar_centered`

  这是Haar小波。仅有的对于Harr小波的$$k$$的有效选择是$$k=2$$。

  `gsl_wavelet_bspline` `gsl_wavelet_bspline_centered`

  这是双正交的具有阶数为$$(i, j)$$B-样条小波家族。被实现的$$k=100*i+j$$值是103, 105, 202, 204, 206, 208, 301, 303, 305 307, 309。

The centered forms of the wavelets align the coefficients of the various sub-bands on edges. Thus the resulting visualization of the coefficients of the wavelet transform in the phase plane is easier to understand.

- const char * `gsl_wavelet_name`(const [gsl_wavelet](https://www.gnu.org/software/gsl/doc/html/dwt.html#c.gsl_wavelet) * *w*)

  This function returns a pointer to the name of the wavelet family for `w`.

- void `gsl_wavelet_free`([gsl_wavelet](https://www.gnu.org/software/gsl/doc/html/dwt.html#c.gsl_wavelet) * *w*)

  This function frees the wavelet object `w`.

- `gsl_wavelet_workspace`

  This structure contains scratch space of the same size as the input data and is used to hold intermediate results during the transform.

- [gsl_wavelet_workspace](https://www.gnu.org/software/gsl/doc/html/dwt.html#c.gsl_wavelet_workspace) * `gsl_wavelet_workspace_alloc`(size_t *n*)

  This function allocates a workspace for the discrete wavelet transform. To perform a one-dimensional transform on `n` elements, a workspace of size `n` must be provided. For two-dimensional transforms of `n`-by-`n` matrices it is sufficient to allocate a workspace of size `n`, since the transform operates on individual rows and columns. A null pointer is returned if insufficient memory is available.

- void `gsl_wavelet_workspace_free`([gsl_wavelet_workspace](https://www.gnu.org/software/gsl/doc/html/dwt.html#c.gsl_wavelet_workspace) * *work*)

  This function frees the allocated workspace `work`.

## Transform Functions

This sections describes the actual functions performing the discrete wavelet transform. Note that the transforms use periodic boundary conditions. If the signal is not periodic in the sample length then spurious coefficients will appear at the beginning and end of each level of the transform.



### Wavelet transforms in one dimension

- int `gsl_wavelet_transform`(const [gsl_wavelet](https://www.gnu.org/software/gsl/doc/html/dwt.html#c.gsl_wavelet) * *w*, double * *data*, size_t *stride*, size_t *n*, gsl_wavelet_direction *dir*, [gsl_wavelet_workspace](https://www.gnu.org/software/gsl/doc/html/dwt.html#c.gsl_wavelet_workspace) * *work*)

- int `gsl_wavelet_transform_forward`(const [gsl_wavelet](https://www.gnu.org/software/gsl/doc/html/dwt.html#c.gsl_wavelet) * *w*, double * *data*, size_t *stride*, size_t *n*, [gsl_wavelet_workspace](https://www.gnu.org/software/gsl/doc/html/dwt.html#c.gsl_wavelet_workspace) * *work*)

- int `gsl_wavelet_transform_inverse`(const [gsl_wavelet](https://www.gnu.org/software/gsl/doc/html/dwt.html#c.gsl_wavelet) * *w*, double * *data*, size_t *stride*, size_t *n*, [gsl_wavelet_workspace](https://www.gnu.org/software/gsl/doc/html/dwt.html#c.gsl_wavelet_workspace) * *work*)

  These functions compute in-place forward and inverse discrete wavelet transforms of length `n`with stride `stride` on the array `data`. The length of the transform `n` is restricted to powers of two. For the `transform` version of the function the argument `dir` can be either `forward` (![+1](https://www.gnu.org/software/gsl/doc/html/_images/math/a7bc7432fa358292afd63292ed1a8acd0553012b.png)) or `backward` (![-1](https://www.gnu.org/software/gsl/doc/html/_images/math/e26aecc9561355dd0e6fbf2742d0ff0d6b58c864.png)). A workspace `work` of length `n` must be provided.For the forward transform, the elements of the original array are replaced by the discrete wavelet transform ![f_i \rightarrow w_{j,k}](https://www.gnu.org/software/gsl/doc/html/_images/math/aa813843234799c48b45dedfa50f8f97db0cc511.png) in a packed triangular storage layout, where `j` is the index of the level ![j = 0 \dots J-1](https://www.gnu.org/software/gsl/doc/html/_images/math/96c7702eebca894f934737d662b75bd52d4cf15a.png) and `k` is the index of the coefficient within each level, ![k = 0 \dots 2^j - 1](https://www.gnu.org/software/gsl/doc/html/_images/math/6bb3c7410a77b118b9cde958b21e35b9540a867e.png). The total number of levels is ![J = \log_2(n)](https://www.gnu.org/software/gsl/doc/html/_images/math/4396aeddb3e5bd5b13f7787561f348fe4cf8639f.png). The output data has the following form,![(s_{-1,0}, d_{0,0}, d_{1,0}, d_{1,1}, d_{2,0},\cdots, d_{j,k},\cdots, d_{J-1,2^{J-1} - 1})](https://www.gnu.org/software/gsl/doc/html/_images/math/640df446fead6855e9c041886e6012bf1d76fca1.png)where the first element is the smoothing coefficient ![s_{-1,0}](https://www.gnu.org/software/gsl/doc/html/_images/math/122e1b5777b14527ec937f3e94770d1b191646ca.png), followed by the detail coefficients ![d_{j,k}](https://www.gnu.org/software/gsl/doc/html/_images/math/39c9eba6a2cddb09e26b918fd1e74cb7fb74dcd7.png) for each level ![j](https://www.gnu.org/software/gsl/doc/html/_images/math/7bc2ad053edeecad39b867d60495ac2caa4555dc.png). The backward transform inverts these coefficients to obtain the original data.These functions return a status of `GSL_SUCCESS` upon successful completion. [`GSL_EINVAL`](https://www.gnu.org/software/gsl/doc/html/err.html#c.GSL_EINVAL) is returned if `n` is not an integer power of 2 or if insufficient workspace is provided.



### Wavelet transforms in two dimension

The library provides functions to perform two-dimensional discrete wavelet transforms on square matrices. The matrix dimensions must be an integer power of two. There are two possible orderings of the rows and columns in the two-dimensional wavelet transform, referred to as the “standard” and “non-standard” forms.

The “standard” transform performs a complete discrete wavelet transform on the rows of the matrix, followed by a separate complete discrete wavelet transform on the columns of the resulting row-transformed matrix. This procedure uses the same ordering as a two-dimensional Fourier transform.

The “non-standard” transform is performed in interleaved passes on the rows and columns of the matrix for each level of the transform. The first level of the transform is applied to the matrix rows, and then to the matrix columns. This procedure is then repeated across the rows and columns of the data for the subsequent levels of the transform, until the full discrete wavelet transform is complete. The non-standard form of the discrete wavelet transform is typically used in image analysis.

The functions described in this section are declared in the header file `gsl_wavelet2d.h`.

- int `gsl_wavelet2d_transform`(const [gsl_wavelet](https://www.gnu.org/software/gsl/doc/html/dwt.html#c.gsl_wavelet) * *w*, double * *data*, size_t *tda*, size_t *size1*, size_t *size2*, gsl_wavelet_direction *dir*, [gsl_wavelet_workspace](https://www.gnu.org/software/gsl/doc/html/dwt.html#c.gsl_wavelet_workspace) * *work*)

- int `gsl_wavelet2d_transform_forward`(const [gsl_wavelet](https://www.gnu.org/software/gsl/doc/html/dwt.html#c.gsl_wavelet) * *w*, double * *data*, size_t *tda*, size_t *size1*, size_t *size2*, [gsl_wavelet_workspace](https://www.gnu.org/software/gsl/doc/html/dwt.html#c.gsl_wavelet_workspace) * *work*)

- int `gsl_wavelet2d_transform_inverse`(const [gsl_wavelet](https://www.gnu.org/software/gsl/doc/html/dwt.html#c.gsl_wavelet) * *w*, double * *data*, size_t *tda*, size_t *size1*, size_t *size2*, [gsl_wavelet_workspace](https://www.gnu.org/software/gsl/doc/html/dwt.html#c.gsl_wavelet_workspace) * *work*)

  These functions compute two-dimensional in-place forward and inverse discrete wavelet transforms in standard form on the array `data` stored in row-major form with dimensions `size1` and `size2` and physical row length `tda`. The dimensions must be equal (square matrix) and are restricted to powers of two. For the `transform` version of the function the argument `dir` can be either `forward` (![+1](https://www.gnu.org/software/gsl/doc/html/_images/math/a7bc7432fa358292afd63292ed1a8acd0553012b.png)) or `backward` (![-1](https://www.gnu.org/software/gsl/doc/html/_images/math/e26aecc9561355dd0e6fbf2742d0ff0d6b58c864.png)). A workspace `work` of the appropriate size must be provided. On exit, the appropriate elements of the array `data` are replaced by their two-dimensional wavelet transform.The functions return a status of `GSL_SUCCESS` upon successful completion. [`GSL_EINVAL`](https://www.gnu.org/software/gsl/doc/html/err.html#c.GSL_EINVAL) is returned if `size1` and `size2` are not equal and integer powers of 2, or if insufficient workspace is provided.

- int `gsl_wavelet2d_transform_matrix`(const [gsl_wavelet](https://www.gnu.org/software/gsl/doc/html/dwt.html#c.gsl_wavelet) * *w*, [gsl_matrix](https://www.gnu.org/software/gsl/doc/html/vectors.html#c.gsl_matrix) * *m*, gsl_wavelet_direction *dir*, [gsl_wavelet_workspace](https://www.gnu.org/software/gsl/doc/html/dwt.html#c.gsl_wavelet_workspace) * *work*)

- int `gsl_wavelet2d_transform_matrix_forward`(const [gsl_wavelet](https://www.gnu.org/software/gsl/doc/html/dwt.html#c.gsl_wavelet) * *w*, [gsl_matrix](https://www.gnu.org/software/gsl/doc/html/vectors.html#c.gsl_matrix) * *m*, [gsl_wavelet_workspace](https://www.gnu.org/software/gsl/doc/html/dwt.html#c.gsl_wavelet_workspace) * *work*)

- int `gsl_wavelet2d_transform_matrix_inverse`(const [gsl_wavelet](https://www.gnu.org/software/gsl/doc/html/dwt.html#c.gsl_wavelet) * *w*, [gsl_matrix](https://www.gnu.org/software/gsl/doc/html/vectors.html#c.gsl_matrix) * *m*, [gsl_wavelet_workspace](https://www.gnu.org/software/gsl/doc/html/dwt.html#c.gsl_wavelet_workspace) * *work*)

  These functions compute the two-dimensional in-place wavelet transform on a matrix `m`.

- int `gsl_wavelet2d_nstransform`(const [gsl_wavelet](https://www.gnu.org/software/gsl/doc/html/dwt.html#c.gsl_wavelet) * *w*, double * *data*, size_t *tda*, size_t *size1*, size_t *size2*, gsl_wavelet_direction *dir*, [gsl_wavelet_workspace](https://www.gnu.org/software/gsl/doc/html/dwt.html#c.gsl_wavelet_workspace) * *work*)

- int `gsl_wavelet2d_nstransform_forward`(const [gsl_wavelet](https://www.gnu.org/software/gsl/doc/html/dwt.html#c.gsl_wavelet) * *w*, double * *data*, size_t *tda*, size_t *size1*, size_t *size2*, [gsl_wavelet_workspace](https://www.gnu.org/software/gsl/doc/html/dwt.html#c.gsl_wavelet_workspace) * *work*)

- int `gsl_wavelet2d_nstransform_inverse`(const [gsl_wavelet](https://www.gnu.org/software/gsl/doc/html/dwt.html#c.gsl_wavelet) * *w*, double * *data*, size_t *tda*, size_t *size1*, size_t *size2*, [gsl_wavelet_workspace](https://www.gnu.org/software/gsl/doc/html/dwt.html#c.gsl_wavelet_workspace) * *work*)

  These functions compute the two-dimensional wavelet transform in non-standard form.

- int `gsl_wavelet2d_nstransform_matrix`(const [gsl_wavelet](https://www.gnu.org/software/gsl/doc/html/dwt.html#c.gsl_wavelet) * *w*, [gsl_matrix](https://www.gnu.org/software/gsl/doc/html/vectors.html#c.gsl_matrix) * *m*, gsl_wavelet_direction *dir*, [gsl_wavelet_workspace](https://www.gnu.org/software/gsl/doc/html/dwt.html#c.gsl_wavelet_workspace) * *work*)

- int `gsl_wavelet2d_nstransform_matrix_forward`(const [gsl_wavelet](https://www.gnu.org/software/gsl/doc/html/dwt.html#c.gsl_wavelet) * *w*, [gsl_matrix](https://www.gnu.org/software/gsl/doc/html/vectors.html#c.gsl_matrix) * *m*, [gsl_wavelet_workspace](https://www.gnu.org/software/gsl/doc/html/dwt.html#c.gsl_wavelet_workspace) * *work*)

- int `gsl_wavelet2d_nstransform_matrix_inverse`(const [gsl_wavelet](https://www.gnu.org/software/gsl/doc/html/dwt.html#c.gsl_wavelet) * *w*, [gsl_matrix](https://www.gnu.org/software/gsl/doc/html/vectors.html#c.gsl_matrix) * *m*, [gsl_wavelet_workspace](https://www.gnu.org/software/gsl/doc/html/dwt.html#c.gsl_wavelet_workspace) * *work*)

  These functions compute the non-standard form of the two-dimensional in-place wavelet transform on a matrix `m`.

## Examples

The following program demonstrates the use of the one-dimensional wavelet transform functions. It computes an approximation to an input signal (of length 256) using the 20 largest components of the wavelet transform, while setting the others to zero.

```
#include <stdio.h>
#include <math.h>
#include <gsl/gsl_sort.h>
#include <gsl/gsl_wavelet.h>

int
main (int argc, char **argv)
{
  (void)(argc); /* avoid unused parameter warning */
  int i, n = 256, nc = 20;
  double *orig_data = malloc (n * sizeof (double));
  double *data = malloc (n * sizeof (double));
  double *abscoeff = malloc (n * sizeof (double));
  size_t *p = malloc (n * sizeof (size_t));

  FILE * f;
  gsl_wavelet *w;
  gsl_wavelet_workspace *work;

  w = gsl_wavelet_alloc (gsl_wavelet_daubechies, 4);
  work = gsl_wavelet_workspace_alloc (n);

  f = fopen (argv[1], "r");
  for (i = 0; i < n; i++)
    {
      fscanf (f, "%lg", &orig_data[i]);
      data[i] = orig_data[i];
    }
  fclose (f);

  gsl_wavelet_transform_forward (w, data, 1, n, work);

  for (i = 0; i < n; i++)
    {
      abscoeff[i] = fabs (data[i]);
    }

  gsl_sort_index (p, abscoeff, 1, n);

  for (i = 0; (i + nc) < n; i++)
    data[p[i]] = 0;

  gsl_wavelet_transform_inverse (w, data, 1, n, work);

  for (i = 0; i < n; i++)
    {
      printf ("%g %g\n", orig_data[i], data[i]);
    }

  gsl_wavelet_free (w);
  gsl_wavelet_workspace_free (work);

  free (data);
  free (orig_data);
  free (abscoeff);
  free (p);
  return 0;
}
```

The output can be used with the GNU plotutils `graph` program:

```
$ ./a.out ecg.dat > dwt.txt
$ graph -T ps -x 0 256 32 -h 0.3 -a dwt.txt > dwt.ps
```

[Fig. 26](https://www.gnu.org/software/gsl/doc/html/dwt.html#fig-dwt) shows an original and compressed version of a sample ECG recording from the MIT-BIH Arrhythmia Database, part of the PhysioNet archive of public-domain of medical datasets.



![_images/dwt.png](https://www.gnu.org/software/gsl/doc/html/_images/dwt.png)

Fig. 26 Original (upper) and wavelet-compressed (lower) ECG signals, using the 20 largest components of the Daubechies(4) discrete wavelet transform.

## References and Further Reading

The mathematical background to wavelet transforms is covered in the original lectures by Daubechies,

- Ingrid Daubechies. Ten Lectures on Wavelets. *CBMS-NSF Regional Conference Series in Applied Mathematics* (1992), SIAM, ISBN 0898712742.

An easy to read introduction to the subject with an emphasis on the application of the wavelet transform in various branches of science is,

- Paul S. Addison. *The Illustrated Wavelet Transform Handbook*. Institute of Physics Publishing (2002), ISBN 0750306920.

For extensive coverage of signal analysis by wavelets, wavelet packets and local cosine bases see,

- S. G. Mallat. *A wavelet tour of signal processing* (Second edition). Academic Press (1999), ISBN 012466606X.

The concept of multiresolution analysis underlying the wavelet transform is described in,

- S. G. Mallat. Multiresolution Approximations and Wavelet Orthonormal Bases of L^2(R).*Transactions of the American Mathematical Society*, 315(1), 1989, 69–87.
- S. G. Mallat. A Theory for Multiresolution Signal Decomposition—The Wavelet Representation.*IEEE Transactions on Pattern Analysis and Machine Intelligence*, 11, 1989, 674–693.

The coefficients for the individual wavelet families implemented by the library can be found in the following papers,

- I. Daubechies. Orthonormal Bases of Compactly Supported Wavelets. *Communications on Pure and Applied Mathematics*, 41 (1988) 909–996.
- A. Cohen, I. Daubechies, and J.-C. Feauveau. Biorthogonal Bases of Compactly Supported Wavelets. *Communications on Pure and Applied Mathematics*, 45 (1992) 485–560.

The PhysioNet archive of physiological datasets can be found online at <http://www.physionet.org/>and is described in the following paper,

- Goldberger et al. PhysioBank, PhysioToolkit, and PhysioNet: Components of a New Research Resource for Complex Physiologic Signals. *Circulation* 101(23):e215-e220 2000.