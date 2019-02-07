# Sorting排序

This chapter describes functions for sorting data, both directly and indirectly (using an index). All the functions use the *heapsort* algorithm. Heapsort is an ![O(N \log N)](https://www.gnu.org/software/gsl/doc/html/_images/math/e0fa68606d706eb7042cf5d3e71f416939b46971.png) algorithm which operates in-place and does not require any additional storage. It also provides consistent performance, the running time for its worst-case (ordered data) being not significantly longer than the average and best cases. Note that the heapsort algorithm does not preserve the relative ordering of equal elements—it is an *unstable* sort. However the resulting order of equal elements will be consistent across different platforms when using these functions.

本章描述用于对数据进行排序的函数，包括直接和间接的(使用下标)。所有的函数使用堆排序算法。堆排序是一个$$O(NlogN)$$算法，其在原位进行操作并且不需要任何额外的存储。它也提供一致的性能，其最差的情形(有序数据)的运行时间并不会比平均或者最好的情形显著更长。注意堆排序算法不会保存相等元素的相对顺序—它是一个不稳定排序。然而相同元素的顺序在不同的平台上使用这些函数时会得到一致的结果。

## Sorting objects排序对象

The following function provides a simple alternative to the standard library function `qsort()`. It is intended for systems lacking `qsort()`, not as a replacement for it. The function `qsort()` should be used whenever possible, as it will be faster and can provide stable ordering of equal elements. Documentation for `qsort()` is available in the GNU C Library Reference Manual.

西面的函数提供了一个对标准库函数`qsort()`的简单替代。它为缺失`qsort()`的系统而准备，而不是`qsort()`的替代品。函数`qsort()`应当在可能的情况下首先被使用，因为其更快且能对相同元素提供稳定的排序。`qsort()`的文档在GNU C库参考手册中可以找到。

The functions described in this section are defined in the header file `gsl_heapsort.h`.

这节描述的函数被定义在头文件`gsl_heapsort.h`中。

- void `gsl_heapsort`(void * *array*, size_t *count*, size_t *size*, [gsl_comparison_fn_t](https://www.gnu.org/software/gsl/doc/html/sort.html#c.gsl_comparison_fn_t) *compare*)

  This function sorts the `count` elements of the array `array`, each of size `size`, into ascending order using the comparison function `compare`. The type of the comparison function is defined by`gsl_comparison_fn_t``int (*gsl_comparison_fn_t) (const void * a, const void * b) `A comparison function should return a negative integer if the first argument is less than the second argument, `0` if the two arguments are equal and a positive integer if the first argument is greater than the second argument.For example, the following function can be used to sort doubles into ascending numerical order.`int compare_doubles (const double * a, const double * b) {   if (*a > *b)     return 1;   else if (*a < *b)     return -1;   else     return 0; } `The appropriate function call to perform the sort is:`gsl_heapsort (array, count, sizeof(double), compare_doubles); `Note that unlike `qsort()` the heapsort algorithm cannot be made into a stable sort by pointer arithmetic. The trick of comparing pointers for equal elements in the comparison function does not work for the heapsort algorithm. The heapsort algorithm performs an internal rearrangement of the data which destroys its initial ordering.

  这个函数使用比较函数`compare`对每个元素的大小为`size`的数组`array`的`count`个元素进行升序排序。比较函数的类型由`gsl_comparison_fn_t``int (*gsl_comparison_fn_t) (const void * a, const void * b) `定义。一个比较函数应当在第一个参数小于第二个参数时返回一个负整数，在两个参数相等时返回`0`，在第一个元素大于第二个元素时返回一个正整数。例如，下面的函数可以被用于将double类型排序为升序`int compare_doubles (const double * a, const double * b) {   if (*a > *b)     return 1;   else if (*a < *b)     return -1;   else     return 0; } `用于执行排序的合适的函数调用是:`gsl_heapsort (array, count, sizeof(double), compare_doubles); `注意不像`qsort()`，堆排序算法不能通过指针算术来实现稳定排序。对相同元素的指针进行比较的技巧对于堆排序算法不起作用。堆排序算法将对数据执行内部重排会破坏其初始顺序。

- int `gsl_heapsort_index`(size_t * *p*, const void * *array*, size_t *count*, size_t *size*, [gsl_comparison_fn_t](https://www.gnu.org/software/gsl/doc/html/sort.html#c.gsl_comparison_fn_t) *compare*)

  This function indirectly sorts the `count` elements of the array `array`, each of size `size`, into ascending order using the comparison function `compare`. The resulting permutation is stored in `p`, an array of length `n`. The elements of `p` give the index of the array element which would have been stored in that position if the array had been sorted in place. The first element of `p`gives the index of the least element in `array`, and the last element of `p` gives the index of the greatest element in `array`. The array itself is not changed.

  这个函数使用比较函数`compare`对每个元素大小为`size`的数组`array`的前`count`个元素进行升序排序。置换的结果被存储在一个长度为`n`的数组`p`中。`p`的元素给出了数组元素的下标，其表示了如果数组被原位排序后元素应当存储的位置。`p`的第一个元素给出了`array`中最小元素的下标，`p`的最后一个元素给出了`array`中最大元素的下标。数组本身没有改变。

## Sorting vectors对向量排序

The following functions will sort the elements of an array or vector, either directly or indirectly. They are defined for all real and integer types using the normal suffix rules. For example, the `float`versions of the array functions are `gsl_sort_float()` and `gsl_sort_float_index()`. The corresponding vector functions are `gsl_sort_vector_float()` and `gsl_sort_vector_float_index()`. The prototypes are available in the header files `gsl_sort_float.h` `gsl_sort_vector_float.h`. The complete set of prototypes can be included using the header files `gsl_sort.h` and`gsl_sort_vector.h`.

下面的函数将会对一个数组或者向量的元素进行排序，要么是直接地要么是间接地。它们使用标准的后缀规则为所有实数和整数类型而定义。例如数组函数的`float`版本是`gsl_sort_float()`和`gsl_sort_float_index()`。相应的向量函数是`gsl_sort_vector_float()`和`gsl_sort_vector_float_index()`。函数原型可以从头文件`gsl_sort_float.h`和`gsl_sort_vector_float.h`访问。完整的原型可以通过头文件`gsl_sort.h`和`gsl_sort_vector.h`进行包含。

There are no functions for sorting complex arrays or vectors, since the ordering of complex numbers is not uniquely defined. To sort a complex vector by magnitude compute a real vector containing the magnitudes of the complex elements, and sort this vector indirectly. The resulting index gives the appropriate ordering of the original complex vector.



- void `gsl_sort`(double * *data*, const size_t *stride*, size_t *n*)

  This function sorts the `n` elements of the array `data` with stride `stride` into ascending numerical order.

- void `gsl_sort2`(double * *data1*, const size_t *stride1*, double * *data2*, const size_t *stride2*, size_t *n*)

  This function sorts the `n` elements of the array `data1` with stride `stride1` into ascending numerical order, while making the same rearrangement of the array `data2` with stride `stride2`, also of size `n`.

- void `gsl_sort_vector`([gsl_vector](https://www.gnu.org/software/gsl/doc/html/vectors.html#c.gsl_vector) * *v*)

  This function sorts the elements of the vector `v` into ascending numerical order.

- void `gsl_sort_vector2`([gsl_vector](https://www.gnu.org/software/gsl/doc/html/vectors.html#c.gsl_vector) * *v1*, [gsl_vector](https://www.gnu.org/software/gsl/doc/html/vectors.html#c.gsl_vector) * *v2*)

  This function sorts the elements of the vector `v1` into ascending numerical order, while making the same rearrangement of the vector `v2`.



- void `gsl_sort_index`(size_t * *p*, const double * *data*, size_t *stride*, size_t *n*)

  This function indirectly sorts the `n` elements of the array `data` with stride `stride` into ascending order, storing the resulting permutation in `p`. The array `p` must be allocated with a sufficient length to store the `n` elements of the permutation. The elements of `p` give the index of the array element which would have been stored in that position if the array had been sorted in place. The array `data` is not changed.

- int `gsl_sort_vector_index`([gsl_permutation](https://www.gnu.org/software/gsl/doc/html/permutation.html#c.gsl_permutation) * *p*, const [gsl_vector](https://www.gnu.org/software/gsl/doc/html/vectors.html#c.gsl_vector) * *v*)

  This function indirectly sorts the elements of the vector `v` into ascending order, storing the resulting permutation in `p`. The elements of `p` give the index of the vector element which would have been stored in that position if the vector had been sorted in place. The first element of `p` gives the index of the least element in `v`, and the last element of `p` gives the index of the greatest element in `v`. The vector `v` is not changed.

## Selecting the k smallest or largest elements

The functions described in this section select the ![k](https://www.gnu.org/software/gsl/doc/html/_images/math/90645b5c2abc13d7de104fd87abeb2f19406e5f8.png) smallest or largest elements of a data set of size ![N](https://www.gnu.org/software/gsl/doc/html/_images/math/99f8a3c97c19212ca3eee7d00dddc8d68b2d1c27.png). The routines use an ![O(kN)](https://www.gnu.org/software/gsl/doc/html/_images/math/2c137111d9c4dd58d4b1a99e237a0523706c675f.png) direct insertion algorithm which is suited to subsets that are small compared with the total size of the dataset. For example, the routines are useful for selecting the 10 largest values from one million data points, but not for selecting the largest 100,000 values. If the subset is a significant part of the total dataset it may be faster to sort all the elements of the dataset directly with an ![O(N \log N)](https://www.gnu.org/software/gsl/doc/html/_images/math/e0fa68606d706eb7042cf5d3e71f416939b46971.png) algorithm and obtain the smallest or largest values that way.

- int `gsl_sort_smallest`(double * *dest*, size_t *k*, const double * *src*, size_t *stride*, size_t *n*)

  This function copies the `k` smallest elements of the array `src`, of size `n` and stride `stride`, in ascending numerical order into the array `dest`. The size `k` of the subset must be less than or equal to `n`. The data `src` is not modified by this operation.

- int `gsl_sort_largest`(double * *dest*, size_t *k*, const double * *src*, size_t *stride*, size_t *n*)

  This function copies the `k` largest elements of the array `src`, of size `n` and stride `stride`, in descending numerical order into the array `dest`. `k` must be less than or equal to `n`. The data `src` is not modified by this operation.

- int `gsl_sort_vector_smallest`(double * *dest*, size_t *k*, const [gsl_vector](https://www.gnu.org/software/gsl/doc/html/vectors.html#c.gsl_vector) * *v*)

- int `gsl_sort_vector_largest`(double * *dest*, size_t *k*, const [gsl_vector](https://www.gnu.org/software/gsl/doc/html/vectors.html#c.gsl_vector) * *v*)

  These functions copy the `k` smallest or largest elements of the vector `v` into the array `dest`. `k` must be less than or equal to the length of the vector `v`.

The following functions find the indices of the ![k](https://www.gnu.org/software/gsl/doc/html/_images/math/90645b5c2abc13d7de104fd87abeb2f19406e5f8.png) smallest or largest elements of a dataset.

- int `gsl_sort_smallest_index`(size_t * *p*, size_t *k*, const double * *src*, size_t *stride*, size_t *n*)

  This function stores the indices of the `k` smallest elements of the array `src`, of size `n` and stride `stride`, in the array `p`. The indices are chosen so that the corresponding data is in ascending numerical order. `k` must be less than or equal to `n`. The data `src` is not modified by this operation.

- int `gsl_sort_largest_index`(size_t * *p*, size_t *k*, const double * *src*, size_t *stride*, size_t *n*)

  This function stores the indices of the `k` largest elements of the array `src`, of size `n` and stride `stride`, in the array `p`. The indices are chosen so that the corresponding data is in descending numerical order. `k` must be less than or equal to `n`. The data `src` is not modified by this operation.

- int `gsl_sort_vector_smallest_index`(size_t * *p*, size_t *k*, const [gsl_vector](https://www.gnu.org/software/gsl/doc/html/vectors.html#c.gsl_vector) * *v*)

- int `gsl_sort_vector_largest_index`(size_t * *p*, size_t *k*, const [gsl_vector](https://www.gnu.org/software/gsl/doc/html/vectors.html#c.gsl_vector) * *v*)

  These functions store the indices of the `k` smallest or largest elements of the vector `v` in the array `p`. `k` must be less than or equal to the length of the vector `v`.

## Computing the rank

The *rank* of an element is its order in the sorted data. The rank is the inverse of the index permutation, ![p](https://www.gnu.org/software/gsl/doc/html/_images/math/af29d3e8554f748ed3f2d142a3ec9d4e038b9335.png). It can be computed using the following algorithm:

```
for (i = 0; i < p->size; i++)
  {
    size_t pi = p->data[i];
    rank->data[pi] = i;
  }
```

This can be computed directly from the function `gsl_permutation_inverse(rank,p)`.

The following function will print the rank of each element of the vector ![v](https://www.gnu.org/software/gsl/doc/html/_images/math/a3ab041b283bfe2226babbccf0611362d5001323.png):

```
void
print_rank (gsl_vector * v)
{
  size_t i;
  size_t n = v->size;
  gsl_permutation * perm = gsl_permutation_alloc(n);
  gsl_permutation * rank = gsl_permutation_alloc(n);

  gsl_sort_vector_index (perm, v);
  gsl_permutation_inverse (rank, perm);

  for (i = 0; i < n; i++)
    {
      double vi = gsl_vector_get(v, i);
      printf ("element = %d, value = %g, rank = %d\n",
               i, vi, rank->data[i]);
    }

  gsl_permutation_free (perm);
  gsl_permutation_free (rank);
}
```

## Examples

The following example shows how to use the permutation ![p](https://www.gnu.org/software/gsl/doc/html/_images/math/af29d3e8554f748ed3f2d142a3ec9d4e038b9335.png) to print the elements of the vector ![v](https://www.gnu.org/software/gsl/doc/html/_images/math/a3ab041b283bfe2226babbccf0611362d5001323.png) in ascending order:

```
gsl_sort_vector_index (p, v);

for (i = 0; i < v->size; i++)
  {
    double vpi = gsl_vector_get (v, p->data[i]);
    printf ("order = %d, value = %g\n", i, vpi);
  }
```

The next example uses the function [`gsl_sort_smallest()`](https://www.gnu.org/software/gsl/doc/html/sort.html#c.gsl_sort_smallest) to select the 5 smallest numbers from 100000 uniform random variates stored in an array,

```
#include <gsl/gsl_rng.h>
#include <gsl/gsl_sort_double.h>

int
main (void)
{
  const gsl_rng_type * T;
  gsl_rng * r;

  size_t i, k = 5, N = 100000;

  double * x = malloc (N * sizeof(double));
  double * small = malloc (k * sizeof(double));

  gsl_rng_env_setup();

  T = gsl_rng_default;
  r = gsl_rng_alloc (T);

  for (i = 0; i < N; i++)
    {
      x[i] = gsl_rng_uniform(r);
    }

  gsl_sort_smallest (small, k, x, 1, N);

  printf ("%zu smallest values from %zu\n", k, N);

  for (i = 0; i < k; i++)
    {
      printf ("%zu: %.18f\n", i, small[i]);
    }

  free (x);
  free (small);
  gsl_rng_free (r);
  return 0;
}
```

The output lists the 5 smallest values, in ascending order,

```
5 smallest values from 100000
0: 0.000003489200025797
1: 0.000008199829608202
2: 0.000008953968062997
3: 0.000010712770745158
4: 0.000033531803637743
```

## References and Further Reading

The subject of sorting is covered extensively in the following,

- Donald E. Knuth, The Art of Computer Programming: Sorting and Searching (Vol 3, 3rd Ed, 1997), Addison-Wesley, ISBN 0201896850.

The Heapsort algorithm is described in the following book,

- Robert Sedgewick, Algorithms in C, Addison-Wesley, ISBN 0201514257.