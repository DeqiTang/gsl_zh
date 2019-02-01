# 向量和矩阵

The functions described in this chapter provide a simple vector and matrix interface to ordinary C arrays. The memory management of these arrays is implemented using a single underlying type, known as a block. By writing your functions in terms of vectors and matrices you can pass a single structure containing both data and dimensions as an argument without needing additional function parameters. The structures are compatible with the vector and matrix formats used by BLAS routines.

本章描述的函数提供了对于C数组的简单的向量和矩阵接口。这些数组的内存管理是通过使用一个单一的称作block的底层类型来实现的。通过编写你自己的关于向量和矩阵的函数，你可以传递一个单一结构做为实参其包含数据以及维度，而不需要额外的函数参数。这些结构与BLAS程序所使用的关于向量和矩阵的格式是兼容的。

## 数据类型

All the functions are available for each of the standard data-types. The versions for `double` have the prefix `gsl_block`, `gsl_vector` and `gsl_matrix`. Similarly the versions for single-precision `float` arrays have the prefix `gsl_block_float`, `gsl_vector_float` and `gsl_matrix_float`. The full list of available types is given below,

| Prefix                        | Type                |
| ----------------------------- | ------------------- |
| gsl_block                     | double              |
| gsl_block_float               | float               |
| gsl_block_long_double         | long double         |
| gsl_block_int                 | int                 |
| gsl_block_uint                | unsigned int        |
| gsl_block_long                | long                |
| gsl_block_ulong               | unsigned long       |
| gsl_block_short               | short               |
| gsl_block_ushort              | unsigned short      |
| gsl_block_char                | char                |
| gsl_block_uchar               | unsigned char       |
| gsl_block_complex             | complex double      |
| gsl_block_complex_float       | complex float       |
| gsl_block_complex_long_double | complex long double |

所有的函数对于每一个标准数据类型都是可获取的。`double`版本的函数具有前缀`gsl_block`，`gsl_vector`和`gsl_matrix`。类似的单精度`float`数组版本具有前缀`gsl_block_float`，`gsl_vector_float`和`gsl_matrix_float`。可用类型的完整列表给在下面，

| Prefix                        | Type                |
| ----------------------------- | ------------------- |
| gsl_block                     | double              |
| gsl_block_float               | float               |
| gsl_block_long_double         | long double         |
| gsl_block_int                 | int                 |
| gsl_block_uint                | unsigned int        |
| gsl_block_long                | long                |
| gsl_block_ulong               | unsigned long       |
| gsl_block_short               | short               |
| gsl_block_ushort              | unsigned short      |
| gsl_block_char                | char                |
| gsl_block_uchar               | unsigned char       |
| gsl_block_complex             | complex double      |
| gsl_block_complex_float       | complex float       |
| gsl_block_complex_long_double | complex long double |

Corresponding types exist for the `gsl_vector` and `gsl_matrix` functions.

相关类型对于`gsl_vector`和`gsl_matrix`函数是存在的。



## Blocks

For consistency all memory is allocated through a `gsl_block` structure. The structure contains two components, the size of an area of memory and a pointer to the memory. The `gsl_block` structure looks like this,

- `gsl_block`

  `typedef struct {   size_t size;   double * data; } gsl_block; `

为了一致性所有的内存都通过一个`gsl_block`结构进行分配。这个结构包含两个成员，内存区域的大小以及一个指向该内存的指针。`gsl_block`结构看起来像这样，

* `gsl_block`

  `typedef struct {   size_t size;   double * data; } gsl_block; `

Vectors and matrices are made by *slicing* an underlying block. A slice is a set of elements formed from an initial offset and a combination of indices and step-sizes. In the case of a matrix the step-size for the column index represents the row-length. The step-size for a vector is known as the *stride*.

向量和矩阵是通过对下层的block进行切片得到的。一个切片有一个初始偏移量和下标以及步长的组合形成的一个元素的集合。在矩阵情况下列指数的步长代表了行的长度。一个向量的步长被称作跨度。

The functions for allocating and deallocating blocks are defined in `gsl_block.h`.

用于分配和释放blocks的函数声明在`gsl_block.h`中。

### Block 分配

The functions for allocating memory to a block follow the style of `malloc` and `free`. In addition they also perform their own error checking. If there is insufficient memory available to allocate a block then the functions call the GSL error handler (with an error number of [`GSL_ENOMEM`](https://www.gnu.org/software/gsl/doc/html/err.html#c.GSL_ENOMEM)) in addition to returning a null pointer. Thus if you use the library error handler to abort your program then it isn’t necessary to check every `alloc`.

为block分配内存的函数遵循`malloc`和`free`的样式。此外它们也会执行它们自己的错误检查。如果没有足够的内存可以使用来分配一个block，那么函数将会除了返回一个空指针以外调用GSL错误处理器(使用[`GSL_ENOMEM`](https://www.gnu.org/software/gsl/doc/html/err.html#c.GSL_ENOMEM)错误编号)。

- [gsl_block](https://www.gnu.org/software/gsl/doc/html/vectors.html#c.gsl_block) * `gsl_block_alloc`(size_t *n*)

  This function allocates memory for a block of `n` double-precision elements, returning a pointer to the block struct. The block is not initialized and so the values of its elements are undefined. Use the function [`gsl_block_calloc()`](https://www.gnu.org/software/gsl/doc/html/vectors.html#c.gsl_block_calloc) if you want to ensure that all the elements are initialized to zero.Zero-sized requests are valid and return a non-null result. A null pointer is returned if insufficient memory is available to create the block.

  这个为一个具有n个双精度元素的block分配内存，返回一个指向block结构的指针。block不会被初始化并因此其元素的值是未定义的。如果你想要确保所有的元素都被初始化为零，使用函数 [`gsl_block_calloc()`](https://www.gnu.org/software/gsl/doc/html/vectors.html#c.gsl_block_calloc)函数。零尺寸的请求也是有效的并会返回一个非空结果。如果内存不足以创建block将返回一个空指针。

- [gsl_block](https://www.gnu.org/software/gsl/doc/html/vectors.html#c.gsl_block) * `gsl_block_calloc`(size_t *n*)

  This function allocates memory for a block and initializes all the elements of the block to zero.

  这个函数为一个block分配内存并将所有元素初始化为零。

- void `gsl_block_free`([gsl_block](https://www.gnu.org/software/gsl/doc/html/vectors.html#c.gsl_block) * *b*)

  This function frees the memory used by a block `b` previously allocated with [`gsl_block_alloc()`](https://www.gnu.org/software/gsl/doc/html/vectors.html#c.gsl_block_alloc) or [`gsl_block_calloc()`](https://www.gnu.org/software/gsl/doc/html/vectors.html#c.gsl_block_calloc).

  这个函数释放由block`b`使用的之前由[`gsl_block_alloc()`](https://www.gnu.org/software/gsl/doc/html/vectors.html#c.gsl_block_alloc) 或者[`gsl_block_calloc()`](https://www.gnu.org/software/gsl/doc/html/vectors.html#c.gsl_block_calloc)函数分配内存。

  

### Blocks的读入和写出

The library provides functions for reading and writing blocks to a file as binary data or formatted text.

这个库提供函数用于按照二进制数据或者格式化文本来读入写出blocks到一个文件。

- int `gsl_block_fwrite`(FILE * *stream*, const [gsl_block](https://www.gnu.org/software/gsl/doc/html/vectors.html#c.gsl_block) * *b*)

  This function writes the elements of the block `b` to the stream `stream` in binary format. The return value is 0 for success and `GSL_EFAILED` if there was a problem writing to the file. Since the data is written in the native binary format it may not be portable between different architectures.

  这个函数将block`b`的元素以二进制格式写到流`stream`中。成功时返回的值是0，如果写出到文件有错误返回值就是`GSL_EFAILED`。因为数据是按照原声二进制格式来写出，在不同架构上有可能是不可移植的。

- int `gsl_block_fread`(FILE * *stream*, [gsl_block](https://www.gnu.org/software/gsl/doc/html/vectors.html#c.gsl_block) * *b*)

  This function reads into the block `b` from the open stream `stream` in binary format. The block `b` must be preallocated with the correct length since the function uses the size of `b` to determine how many bytes to read. The return value is 0 for success and `GSL_EFAILED` if there was a problem reading from the file. The data is assumed to have been written in the native binary format on the same architecture.

  合格函数从打开的流`stream`中以二进制格式读入block`b`。block`b`必须提前分配为正确的长度，因为函数会使用`b`的尺寸来决定读入多少字节。如果成功返回的值是0，如果从文件读入时发生错误将返回`GSL_EFAILED`。在同样的架构上数据被认为是写入为原生二进制格式。

- int `gsl_block_fprintf`(FILE * *stream*, const [gsl_block](https://www.gnu.org/software/gsl/doc/html/vectors.html#c.gsl_block) * *b*, const char * *format*)

  This function writes the elements of the block `b` line-by-line to the stream `stream` using the format specifier `format`, which should be one of the `%g`, `%e` or `%f` formats for floating point numbers and `%d` for integers. The function returns 0 for success and `GSL_EFAILED` if there was a problem writing to the file.

  这个函数将block`b`的元素按行使用格式声明`format`来写入到流`stream`中，格式声明应当是对于浮点数为`%g`，`%e`或者`%f`格式之一对于整数是`%d`格式。如果成功函数返回0，而当写入文件出现问题时返回`GSL_EFAILED`。

- int `gsl_block_fscanf`(FILE * *stream*, [gsl_block](https://www.gnu.org/software/gsl/doc/html/vectors.html#c.gsl_block) * *b*)

  This function reads formatted data from the stream `stream` into the block `b`. The block `b`must be preallocated with the correct length since the function uses the size of `b` to determine how many numbers to read. The function returns 0 for success and `GSL_EFAILED` if there was a problem reading from the file.

  这个函数从流`stream`中读取格式化数据到block`b`中。block`b`必须提前分配为正确的尺寸，因为函数将使用`b`的尺寸来决定读取多少数字。如果成功函数返回0，而当从文件读取发生错误时返回`GSL_EFAILED`。

  

### 用于block的示例程序

The following program shows how to allocate a block,

下列程序展示如何分配一个block，

```
#include <stdio.h>
#include <gsl/gsl_block.h>

int
main (void)
{
  gsl_block * b = gsl_block_alloc (100);

  printf ("length of block = %zu\n", b->size);
  printf ("block data address = %p\n", b->data);

  gsl_block_free (b);
  return 0;
}
```

Here is the output from the program,

下面是程序的输出，

```
length of block = 100
block data address = 0x804b0d8
```



## 向量

Vectors are defined by a [`gsl_vector`](https://www.gnu.org/software/gsl/doc/html/vectors.html#c.gsl_vector) structure which describes a slice of a block. Different vectors can be created which point to the same block. A vector slice is a set of equally-spaced elements of an area of memory.

向量被结构`gsl_vector`定义，其描述了一个block的切片。可以创建指向统一个block的不同向量。一个向量切片是一个由一块内存区域的等间距元素的集合。

The [`gsl_vector`](https://www.gnu.org/software/gsl/doc/html/vectors.html#c.gsl_vector) structure contains five components, the *size*, the *stride*, a pointer to the memory where the elements are stored, `data`, a pointer to the block owned by the vector, `block`, if any, and an ownership flag, `owner`. The structure is very simple and looks like this,

- `gsl_vector`

  `typedef struct {   size_t size;   size_t stride;   double * data;   gsl_block * block;   int owner; } gsl_vector; `

`gsl_vector`结构包含五个成员，尺寸，步长，一个指向存储有元素的内存区域的指针`data`，一个指向由向量拥有的block的指针`block`，若有的话，还有一个拥有权标识，`owner`。结构非常简单且看起来就像这样，

* `gsl_vector`

  `typedef struct {   size_t size;   size_t stride;   double * data;   gsl_block * block;   int owner; } gsl_vector; `

The `size` is simply the number of vector elements. The range of valid indices runs from 0 to `size-1`. The `stride` is the step-size from one element to the next in physical memory, measured in units of the appropriate datatype. The pointer `data` gives the location of the first element of the vector in memory. The pointer `block` stores the location of the memory block in which the vector elements are located (if any). If the vector owns this block then the `owner` field is set to one and the block will be deallocated when the vector is freed. If the vector points to a block owned by another object then the `owner` field is zero and any underlying block will not be deallocated with the vector.

`size`很简单就是向量元素的数量。有效指数的范围经0到`size-1`.`stride`是是物理内存上从一个元素到下一个元素的步长大小，由合适的数据类型作为单位进行衡量。指针`data`给出了内存中向量的第一个元素的位置。指针`block`存储了向量所位于的内存块的位置(如果有的话)。如果向量拥有该block，那么`owner`域被设为1，并且当向量被释放以后该block也会被释放。如果向量指向一个由另一个对象拥有的block，那么`owner`域被设为0，并且任何低层的block都不会随向量一起被释放。

The functions for allocating and accessing vectors are defined in `gsl_vector.h`.

用于分配和访问向量的函数被声明在`gsl_vector.h`中。

### 向量分配

The functions for allocating memory to a vector follow the style of `malloc` and `free`. In addition they also perform their own error checking. If there is insufficient memory available to allocate a vector then the functions call the GSL error handler (with an error number of [`GSL_ENOMEM`](https://www.gnu.org/software/gsl/doc/html/err.html#c.GSL_ENOMEM)) in addition to returning a null pointer. Thus if you use the library error handler to abort your program then it isn’t necessary to check every `alloc`.

- [gsl_vector](https://www.gnu.org/software/gsl/doc/html/vectors.html#c.gsl_vector) * `gsl_vector_alloc`(size_t *n*)

  This function creates a vector of length n, returning a pointer to a newly initialized vector struct. A new block is allocated for the elements of the vector, and stored in the `block` component of the vector struct. The block is “owned” by the vector, and will be deallocated when the vector is deallocated. Zero-sized requests are valid and return a non-null result.

- [gsl_vector](https://www.gnu.org/software/gsl/doc/html/vectors.html#c.gsl_vector) * `gsl_vector_calloc`(size_t *n*)

  This function allocates memory for a vector of length `n` and initializes all the elements of the vector to zero.

- void `gsl_vector_free`([gsl_vector](https://www.gnu.org/software/gsl/doc/html/vectors.html#c.gsl_vector) * *v*)

  This function frees a previously allocated vector `v`. If the vector was created using [`gsl_vector_alloc()`](https://www.gnu.org/software/gsl/doc/html/vectors.html#c.gsl_vector_alloc) then the block underlying the vector will also be deallocated. If the vector has been created from another object then the memory is still owned by that object and will not be deallocated.



### Accessing vector elements

Unlike Fortran compilers, C compilers do not usually provide support for range checking of vectors and matrices. [[1\]](https://www.gnu.org/software/gsl/doc/html/vectors.html#f1) The functions [`gsl_vector_get()`](https://www.gnu.org/software/gsl/doc/html/vectors.html#c.gsl_vector_get) and [`gsl_vector_set()`](https://www.gnu.org/software/gsl/doc/html/vectors.html#c.gsl_vector_set) can perform portable range checking for you and report an error if you attempt to access elements outside the allowed range.

The functions for accessing the elements of a vector or matrix are defined in `gsl_vector.h` and declared `extern inline` to eliminate function-call overhead. You must compile your program with the preprocessor macro `HAVE_INLINE` defined to use these functions.

- `GSL_RANGE_CHECK_OFF`

  If necessary you can turn off range checking completely without modifying any source files by recompiling your program with the preprocessor definition [`GSL_RANGE_CHECK_OFF`](https://www.gnu.org/software/gsl/doc/html/vectors.html#c.GSL_RANGE_CHECK_OFF). Provided your compiler supports inline functions the effect of turning off range checking is to replace calls to `gsl_vector_get(v,i)` by `v->data[i*v->stride]` and calls to `gsl_vector_set(v,i,x)` by`v->data[i*v->stride]=x`. Thus there should be no performance penalty for using the range checking functions when range checking is turned off.

- `GSL_C99_INLINE`

  If you use a C99 compiler which requires inline functions in header files to be declared `inline`instead of `extern inline`, define the macro [`GSL_C99_INLINE`](https://www.gnu.org/software/gsl/doc/html/vectors.html#c.GSL_C99_INLINE) (see [Inline functions](https://www.gnu.org/software/gsl/doc/html/usage.html#sec-inline-functions)). With GCC this is selected automatically when compiling in C99 mode (`-std=c99`).

- `gsl_check_range`

  If inline functions are not used, calls to the functions [`gsl_vector_get()`](https://www.gnu.org/software/gsl/doc/html/vectors.html#c.gsl_vector_get) and [`gsl_vector_set()`](https://www.gnu.org/software/gsl/doc/html/vectors.html#c.gsl_vector_set)will link to the compiled versions of these functions in the library itself. The range checking in these functions is controlled by the global integer variable `gsl_check_range`. It is enabled by default—to disable it, set `gsl_check_range` to zero. Due to function-call overhead, there is less benefit in disabling range checking here than for inline functions.

- double `gsl_vector_get`(const [gsl_vector](https://www.gnu.org/software/gsl/doc/html/vectors.html#c.gsl_vector) * *v*, const size_t *i*)

  This function returns the `i`-th element of a vector `v`. If `i` lies outside the allowed range of 0 to `size - 1` then the error handler is invoked and 0 is returned. An inline version of this function is used when `HAVE_INLINE` is defined.

- void `gsl_vector_set`([gsl_vector](https://www.gnu.org/software/gsl/doc/html/vectors.html#c.gsl_vector) * *v*, const size_t *i*, double *x*)

  This function sets the value of the `i`-th element of a vector `v` to `x`. If `i` lies outside the allowed range of 0 to `size - 1` then the error handler is invoked. An inline version of this function is used when `HAVE_INLINE` is defined.

- double * `gsl_vector_ptr`([gsl_vector](https://www.gnu.org/software/gsl/doc/html/vectors.html#c.gsl_vector) * *v*, size_t *i*)

- const double * `gsl_vector_const_ptr`(const [gsl_vector](https://www.gnu.org/software/gsl/doc/html/vectors.html#c.gsl_vector) * *v*, size_t *i*)

  These functions return a pointer to the `i`-th element of a vector `v`. If `i` lies outside the allowed range of 0 to `size - 1` then the error handler is invoked and a null pointer is returned. Inline versions of these functions are used when `HAVE_INLINE` is defined.



### Initializing vector elements

- void `gsl_vector_set_all`([gsl_vector](https://www.gnu.org/software/gsl/doc/html/vectors.html#c.gsl_vector) * *v*, double *x*)

  This function sets all the elements of the vector `v` to the value `x`.

- void `gsl_vector_set_zero`([gsl_vector](https://www.gnu.org/software/gsl/doc/html/vectors.html#c.gsl_vector) * *v*)

  This function sets all the elements of the vector `v` to zero.

- int `gsl_vector_set_basis`([gsl_vector](https://www.gnu.org/software/gsl/doc/html/vectors.html#c.gsl_vector) * *v*, size_t *i*)

  This function makes a basis vector by setting all the elements of the vector `v` to zero except for the `i`-th element which is set to one.

### Reading and writing vectors

The library provides functions for reading and writing vectors to a file as binary data or formatted text.

- int `gsl_vector_fwrite`(FILE * *stream*, const [gsl_vector](https://www.gnu.org/software/gsl/doc/html/vectors.html#c.gsl_vector) * *v*)

  This function writes the elements of the vector `v` to the stream `stream` in binary format. The return value is 0 for success and `GSL_EFAILED` if there was a problem writing to the file. Since the data is written in the native binary format it may not be portable between different architectures.

- int `gsl_vector_fread`(FILE * *stream*, [gsl_vector](https://www.gnu.org/software/gsl/doc/html/vectors.html#c.gsl_vector) * *v*)

  This function reads into the vector `v` from the open stream `stream` in binary format. The vector `v` must be preallocated with the correct length since the function uses the size of `v` to determine how many bytes to read. The return value is 0 for success and `GSL_EFAILED` if there was a problem reading from the file. The data is assumed to have been written in the native binary format on the same architecture.

- int `gsl_vector_fprintf`(FILE * *stream*, const [gsl_vector](https://www.gnu.org/software/gsl/doc/html/vectors.html#c.gsl_vector) * *v*, const char * *format*)

  This function writes the elements of the vector `v` line-by-line to the stream `stream` using the format specifier `format`, which should be one of the `%g`, `%e` or `%f` formats for floating point numbers and `%d` for integers. The function returns 0 for success and `GSL_EFAILED` if there was a problem writing to the file.

- int `gsl_vector_fscanf`(FILE * *stream*, [gsl_vector](https://www.gnu.org/software/gsl/doc/html/vectors.html#c.gsl_vector) * *v*)

  This function reads formatted data from the stream `stream` into the vector `v`. The vector `v`must be preallocated with the correct length since the function uses the size of `v` to determine how many numbers to read. The function returns 0 for success and `GSL_EFAILED` if there was a problem reading from the file.

### Vector views

In addition to creating vectors from slices of blocks it is also possible to slice vectors and create vector views. For example, a subvector of another vector can be described with a view, or two views can be made which provide access to the even and odd elements of a vector.

- `gsl_vector_view`

- `gsl_vector_const_view`

  A vector view is a temporary object, stored on the stack, which can be used to operate on a subset of vector elements. Vector views can be defined for both constant and non-constant vectors, using separate types that preserve constness. A vector view has the type[`gsl_vector_view`](https://www.gnu.org/software/gsl/doc/html/vectors.html#c.gsl_vector_view) and a constant vector view has the type [`gsl_vector_const_view`](https://www.gnu.org/software/gsl/doc/html/vectors.html#c.gsl_vector_const_view). In both cases the elements of the view can be accessed as a [`gsl_vector`](https://www.gnu.org/software/gsl/doc/html/vectors.html#c.gsl_vector) using the `vector` component of the view object. A pointer to a vector of type `gsl_vector *` or `const gsl_vector *` can be obtained by taking the address of this component with the `&` operator.When using this pointer it is important to ensure that the view itself remains in scope—the simplest way to do so is by always writing the pointer as `&view.vector`, and never storing this value in another variable.

- [gsl_vector_view](https://www.gnu.org/software/gsl/doc/html/vectors.html#c.gsl_vector_view) `gsl_vector_subvector`([gsl_vector](https://www.gnu.org/software/gsl/doc/html/vectors.html#c.gsl_vector) * *v*, size_t *offset*, size_t *n*)

- [gsl_vector_const_view](https://www.gnu.org/software/gsl/doc/html/vectors.html#c.gsl_vector_const_view) `gsl_vector_const_subvector`(const [gsl_vector](https://www.gnu.org/software/gsl/doc/html/vectors.html#c.gsl_vector) * *v*, size_t *offset*, size_t *n*)

  These functions return a vector view of a subvector of another vector `v`. The start of the new vector is offset by `offset` elements from the start of the original vector. The new vector has `n`elements. Mathematically, the `i`-th element of the new vector `v'` is given by:`v'(i) = v->data[(offset + i)*v->stride] `where the index `i` runs from 0 to `n - 1`.The `data` pointer of the returned vector struct is set to null if the combined parameters (`offset`, `n`) overrun the end of the original vector.The new vector is only a view of the block underlying the original vector, `v`. The block containing the elements of `v` is not owned by the new vector. When the view goes out of scope the original vector `v` and its block will continue to exist. The original memory can only be deallocated by freeing the original vector. Of course, the original vector should not be deallocated while the view is still in use.The function [`gsl_vector_const_subvector()`](https://www.gnu.org/software/gsl/doc/html/vectors.html#c.gsl_vector_const_subvector) is equivalent to [`gsl_vector_subvector()`](https://www.gnu.org/software/gsl/doc/html/vectors.html#c.gsl_vector_subvector) but can be used for vectors which are declared `const`.

- [gsl_vector_view](https://www.gnu.org/software/gsl/doc/html/vectors.html#c.gsl_vector_view) `gsl_vector_subvector_with_stride`([gsl_vector](https://www.gnu.org/software/gsl/doc/html/vectors.html#c.gsl_vector) * *v*, size_t *offset*, size_t *stride*, size_t *n*)

- [gsl_vector_const_view](https://www.gnu.org/software/gsl/doc/html/vectors.html#c.gsl_vector_const_view) `gsl_vector_const_subvector_with_stride`(const [gsl_vector](https://www.gnu.org/software/gsl/doc/html/vectors.html#c.gsl_vector) * *v*, size_t *offset*, size_t *stride*, size_t *n*)

  These functions return a vector view of a subvector of another vector `v` with an additional stride argument. The subvector is formed in the same way as for [`gsl_vector_subvector()`](https://www.gnu.org/software/gsl/doc/html/vectors.html#c.gsl_vector_subvector) but the new vector has `n` elements with a step-size of `stride` from one element to the next in the original vector. Mathematically, the `i`-th element of the new vector `v'` is given by:`v'(i) = v->data[(offset + i*stride)*v->stride] `where the index `i` runs from 0 to `n - 1`.Note that subvector views give direct access to the underlying elements of the original vector. For example, the following code will zero the even elements of the vector `v` of length `n`, while leaving the odd elements untouched:`gsl_vector_view v_even = gsl_vector_subvector_with_stride (v, 0, 2, n/2); gsl_vector_set_zero (&v_even.vector); `A vector view can be passed to any subroutine which takes a vector argument just as a directly allocated vector would be, using `&view.vector`. For example, the following code computes the norm of the odd elements of `v` using the BLAS routine `dnrm2`:`gsl_vector_view v_odd = gsl_vector_subvector_with_stride (v, 1, 2, n/2); double r = gsl_blas_dnrm2 (&v_odd.vector); `The function [`gsl_vector_const_subvector_with_stride()`](https://www.gnu.org/software/gsl/doc/html/vectors.html#c.gsl_vector_const_subvector_with_stride) is equivalent to [`gsl_vector_subvector_with_stride()`](https://www.gnu.org/software/gsl/doc/html/vectors.html#c.gsl_vector_subvector_with_stride) but can be used for vectors which are declared `const`.

- [gsl_vector_view](https://www.gnu.org/software/gsl/doc/html/vectors.html#c.gsl_vector_view) `gsl_vector_complex_real`(gsl_vector_complex * *v*)

- [gsl_vector_const_view](https://www.gnu.org/software/gsl/doc/html/vectors.html#c.gsl_vector_const_view) `gsl_vector_complex_const_real`(const gsl_vector_complex * *v*)

  These functions return a vector view of the real parts of the complex vector `v`.The function [`gsl_vector_complex_const_real()`](https://www.gnu.org/software/gsl/doc/html/vectors.html#c.gsl_vector_complex_const_real) is equivalent to [`gsl_vector_complex_real()`](https://www.gnu.org/software/gsl/doc/html/vectors.html#c.gsl_vector_complex_real) but can be used for vectors which are declared `const`.

- [gsl_vector_view](https://www.gnu.org/software/gsl/doc/html/vectors.html#c.gsl_vector_view) `gsl_vector_complex_imag`(gsl_vector_complex * *v*)

- [gsl_vector_const_view](https://www.gnu.org/software/gsl/doc/html/vectors.html#c.gsl_vector_const_view) `gsl_vector_complex_const_imag`(const gsl_vector_complex * *v*)

  These functions return a vector view of the imaginary parts of the complex vector `v`.The function [`gsl_vector_complex_const_imag()`](https://www.gnu.org/software/gsl/doc/html/vectors.html#c.gsl_vector_complex_const_imag) is equivalent to [`gsl_vector_complex_imag()`](https://www.gnu.org/software/gsl/doc/html/vectors.html#c.gsl_vector_complex_imag) but can be used for vectors which are declared `const`.

- [gsl_vector_view](https://www.gnu.org/software/gsl/doc/html/vectors.html#c.gsl_vector_view) `gsl_vector_view_array`(double * *base*, size_t *n*)

- [gsl_vector_const_view](https://www.gnu.org/software/gsl/doc/html/vectors.html#c.gsl_vector_const_view) `gsl_vector_const_view_array`(const double * *base*, size_t *n*)

  These functions return a vector view of an array. The start of the new vector is given by `base`and has `n` elements. Mathematically, the `i`-th element of the new vector `v'` is given by:`v'(i) = base[i] `where the index `i` runs from 0 to `n - 1`.The array containing the elements of `v` is not owned by the new vector view. When the view goes out of scope the original array will continue to exist. The original memory can only be deallocated by freeing the original pointer `base`. Of course, the original array should not be deallocated while the view is still in use.The function [`gsl_vector_const_view_array()`](https://www.gnu.org/software/gsl/doc/html/vectors.html#c.gsl_vector_const_view_array) is equivalent to [`gsl_vector_view_array()`](https://www.gnu.org/software/gsl/doc/html/vectors.html#c.gsl_vector_view_array) but can be used for arrays which are declared `const`.

- [gsl_vector_view](https://www.gnu.org/software/gsl/doc/html/vectors.html#c.gsl_vector_view) `gsl_vector_view_array_with_stride`(double * *base*, size_t *stride*, size_t *n*)

- [gsl_vector_const_view](https://www.gnu.org/software/gsl/doc/html/vectors.html#c.gsl_vector_const_view) `gsl_vector_const_view_array_with_stride`(const double * *base*, size_t *stride*, size_t *n*)

  These functions return a vector view of an array `base` with an additional stride argument. The subvector is formed in the same way as for [`gsl_vector_view_array()`](https://www.gnu.org/software/gsl/doc/html/vectors.html#c.gsl_vector_view_array) but the new vector has `n`elements with a step-size of `stride` from one element to the next in the original array. Mathematically, the `i`-th element of the new vector `v'` is given by:`v'(i) = base[i*stride] `where the index `i` runs from 0 to `n - 1`.Note that the view gives direct access to the underlying elements of the original array. A vector view can be passed to any subroutine which takes a vector argument just as a directly allocated vector would be, using `&view.vector`.The function [`gsl_vector_const_view_array_with_stride()`](https://www.gnu.org/software/gsl/doc/html/vectors.html#c.gsl_vector_const_view_array_with_stride) is equivalent to [`gsl_vector_view_array_with_stride()`](https://www.gnu.org/software/gsl/doc/html/vectors.html#c.gsl_vector_view_array_with_stride) but can be used for arrays which are declared `const`.

### Copying vectors

Common operations on vectors such as addition and multiplication are available in the BLAS part of the library (see [BLAS Support](https://www.gnu.org/software/gsl/doc/html/blas.html#chap-blas-support)). However, it is useful to have a small number of utility functions which do not require the full BLAS code. The following functions fall into this category.

- int `gsl_vector_memcpy`([gsl_vector](https://www.gnu.org/software/gsl/doc/html/vectors.html#c.gsl_vector) * *dest*, const [gsl_vector](https://www.gnu.org/software/gsl/doc/html/vectors.html#c.gsl_vector) * *src*)

  This function copies the elements of the vector `src` into the vector `dest`. The two vectors must have the same length.

- int `gsl_vector_swap`([gsl_vector](https://www.gnu.org/software/gsl/doc/html/vectors.html#c.gsl_vector) * *v*, [gsl_vector](https://www.gnu.org/software/gsl/doc/html/vectors.html#c.gsl_vector) * *w*)

  This function exchanges the elements of the vectors `v` and `w` by copying. The two vectors must have the same length.

### Exchanging elements

The following functions can be used to exchange, or permute, the elements of a vector.

- int `gsl_vector_swap_elements`([gsl_vector](https://www.gnu.org/software/gsl/doc/html/vectors.html#c.gsl_vector) * *v*, size_t *i*, size_t *j*)

  This function exchanges the `i`-th and `j`-th elements of the vector `v` in-place.

- int `gsl_vector_reverse`([gsl_vector](https://www.gnu.org/software/gsl/doc/html/vectors.html#c.gsl_vector) * *v*)

  This function reverses the order of the elements of the vector `v`.

### Vector operations

- int `gsl_vector_add`([gsl_vector](https://www.gnu.org/software/gsl/doc/html/vectors.html#c.gsl_vector) * *a*, const [gsl_vector](https://www.gnu.org/software/gsl/doc/html/vectors.html#c.gsl_vector) * *b*)

  This function adds the elements of vector `b` to the elements of vector `a`. The result ![a_i \leftarrow a_i + b_i](https://www.gnu.org/software/gsl/doc/html/_images/math/352e2b12fe9de5e7ae5391b12865465ebe6f7af2.png) is stored in `a` and `b` remains unchanged. The two vectors must have the same length.

- int `gsl_vector_sub`([gsl_vector](https://www.gnu.org/software/gsl/doc/html/vectors.html#c.gsl_vector) * *a*, const [gsl_vector](https://www.gnu.org/software/gsl/doc/html/vectors.html#c.gsl_vector) * *b*)

  This function subtracts the elements of vector `b` from the elements of vector `a`. The result ![a_i \leftarrow a_i - b_i](https://www.gnu.org/software/gsl/doc/html/_images/math/4638c1964aac587c42c01bc55746cb5fe84099ed.png) is stored in `a` and `b` remains unchanged. The two vectors must have the same length.

- int `gsl_vector_mul`([gsl_vector](https://www.gnu.org/software/gsl/doc/html/vectors.html#c.gsl_vector) * *a*, const [gsl_vector](https://www.gnu.org/software/gsl/doc/html/vectors.html#c.gsl_vector) * *b*)

  This function multiplies the elements of vector `a` by the elements of vector `b`. The result ![a_i \leftarrow a_i * b_i](https://www.gnu.org/software/gsl/doc/html/_images/math/2b7200dfe676821322d687223cf70261f84c6c24.png) is stored in `a` and `b` remains unchanged. The two vectors must have the same length.

- int `gsl_vector_div`([gsl_vector](https://www.gnu.org/software/gsl/doc/html/vectors.html#c.gsl_vector) * *a*, const [gsl_vector](https://www.gnu.org/software/gsl/doc/html/vectors.html#c.gsl_vector) * *b*)

  This function divides the elements of vector `a` by the elements of vector `b`. The result ![a_i \leftarrow a_i / b_i](https://www.gnu.org/software/gsl/doc/html/_images/math/2a44d2abda16f4aa3e28eb9a7716ceff0cea8e0b.png) is stored in `a` and `b` remains unchanged. The two vectors must have the same length.

- int `gsl_vector_scale`([gsl_vector](https://www.gnu.org/software/gsl/doc/html/vectors.html#c.gsl_vector) * *a*, const double *x*)

  This function multiplies the elements of vector `a` by the constant factor `x`. The result ![a_i \leftarrow x a_i](https://www.gnu.org/software/gsl/doc/html/_images/math/f8372f9fee22038c59e7b5788b461b067159dc1b.png) is stored in `a`.

- int `gsl_vector_add_constant`([gsl_vector](https://www.gnu.org/software/gsl/doc/html/vectors.html#c.gsl_vector) * *a*, const double *x*)

  This function adds the constant value `x` to the elements of the vector `a`. The result ![a_i \leftarrow a_i + x](https://www.gnu.org/software/gsl/doc/html/_images/math/66fafd4da5646bb95447ebd368cc645a268dafa9.png) is stored in `a`.

### Finding maximum and minimum elements of vectors

The following operations are only defined for real vectors.

- double `gsl_vector_max`(const [gsl_vector](https://www.gnu.org/software/gsl/doc/html/vectors.html#c.gsl_vector) * *v*)

  This function returns the maximum value in the vector `v`.

- double `gsl_vector_min`(const [gsl_vector](https://www.gnu.org/software/gsl/doc/html/vectors.html#c.gsl_vector) * *v*)

  This function returns the minimum value in the vector `v`.

- void `gsl_vector_minmax`(const [gsl_vector](https://www.gnu.org/software/gsl/doc/html/vectors.html#c.gsl_vector) * *v*, double * *min_out*, double * *max_out*)

  This function returns the minimum and maximum values in the vector `v`, storing them in `min_out` and `max_out`.

- size_t `gsl_vector_max_index`(const [gsl_vector](https://www.gnu.org/software/gsl/doc/html/vectors.html#c.gsl_vector) * *v*)

  This function returns the index of the maximum value in the vector `v`. When there are several equal maximum elements then the lowest index is returned.

- size_t `gsl_vector_min_index`(const [gsl_vector](https://www.gnu.org/software/gsl/doc/html/vectors.html#c.gsl_vector) * *v*)

  This function returns the index of the minimum value in the vector `v`. When there are several equal minimum elements then the lowest index is returned.

- void `gsl_vector_minmax_index`(const [gsl_vector](https://www.gnu.org/software/gsl/doc/html/vectors.html#c.gsl_vector) * *v*, size_t * *imin*, size_t * *imax*)

  This function returns the indices of the minimum and maximum values in the vector `v`, storing them in `imin` and `imax`. When there are several equal minimum or maximum elements then the lowest indices are returned.

### Vector properties

The following functions are defined for real and complex vectors. For complex vectors both the real and imaginary parts must satisfy the conditions.

- int `gsl_vector_isnull`(const [gsl_vector](https://www.gnu.org/software/gsl/doc/html/vectors.html#c.gsl_vector) * *v*)

- int `gsl_vector_ispos`(const [gsl_vector](https://www.gnu.org/software/gsl/doc/html/vectors.html#c.gsl_vector) * *v*)

- int `gsl_vector_isneg`(const [gsl_vector](https://www.gnu.org/software/gsl/doc/html/vectors.html#c.gsl_vector) * *v*)

- int `gsl_vector_isnonneg`(const [gsl_vector](https://www.gnu.org/software/gsl/doc/html/vectors.html#c.gsl_vector) * *v*)

  These functions return 1 if all the elements of the vector `v` are zero, strictly positive, strictly negative, or non-negative respectively, and 0 otherwise.

- int `gsl_vector_equal`(const [gsl_vector](https://www.gnu.org/software/gsl/doc/html/vectors.html#c.gsl_vector) * *u*, const [gsl_vector](https://www.gnu.org/software/gsl/doc/html/vectors.html#c.gsl_vector) * *v*)

  This function returns 1 if the vectors `u` and `v` are equal (by comparison of element values) and 0 otherwise.

### Example programs for vectors

This program shows how to allocate, initialize and read from a vector using the functions [`gsl_vector_alloc()`](https://www.gnu.org/software/gsl/doc/html/vectors.html#c.gsl_vector_alloc), [`gsl_vector_set()`](https://www.gnu.org/software/gsl/doc/html/vectors.html#c.gsl_vector_set) and [`gsl_vector_get()`](https://www.gnu.org/software/gsl/doc/html/vectors.html#c.gsl_vector_get).

```
#include <stdio.h>
#include <gsl/gsl_vector.h>

int
main (void)
{
  int i;
  gsl_vector * v = gsl_vector_alloc (3);

  for (i = 0; i < 3; i++)
    {
      gsl_vector_set (v, i, 1.23 + i);
    }

  for (i = 0; i < 100; i++) /* OUT OF RANGE ERROR */
    {
      printf ("v_%d = %g\n", i, gsl_vector_get (v, i));
    }

  gsl_vector_free (v);
  return 0;
}
```

Here is the output from the program. The final loop attempts to read outside the range of the vector `v`, and the error is trapped by the range-checking code in [`gsl_vector_get()`](https://www.gnu.org/software/gsl/doc/html/vectors.html#c.gsl_vector_get).

```
$ ./a.out
v_0 = 1.23
v_1 = 2.23
v_2 = 3.23
gsl: vector_source.c:12: ERROR: index out of range
Default GSL error handler invoked.
Aborted (core dumped)
```

The next program shows how to write a vector to a file.

```
#include <stdio.h>
#include <gsl/gsl_vector.h>

int
main (void)
{
  int i;
  gsl_vector * v = gsl_vector_alloc (100);

  for (i = 0; i < 100; i++)
    {
      gsl_vector_set (v, i, 1.23 + i);
    }

  {
     FILE * f = fopen ("test.dat", "w");
     gsl_vector_fprintf (f, v, "%.5g");
     fclose (f);
  }

  gsl_vector_free (v);
  return 0;
}
```

After running this program the file `test.dat` should contain the elements of `v`, written using the format specifier `%.5g`. The vector could then be read back in using the function`gsl_vector_fscanf (f, v)` as follows:

```
#include <stdio.h>
#include <gsl/gsl_vector.h>

int
main (void)
{
  int i;
  gsl_vector * v = gsl_vector_alloc (10);

  {
     FILE * f = fopen ("test.dat", "r");
     gsl_vector_fscanf (f, v);
     fclose (f);
  }

  for (i = 0; i < 10; i++)
    {
      printf ("%g\n", gsl_vector_get(v, i));
    }

  gsl_vector_free (v);
  return 0;
}
```



## Matrices

Matrices are defined by a [`gsl_matrix`](https://www.gnu.org/software/gsl/doc/html/vectors.html#c.gsl_matrix) structure which describes a generalized slice of a block. Like a vector it represents a set of elements in an area of memory, but uses two indices instead of one.

- `gsl_matrix`

  The [`gsl_matrix`](https://www.gnu.org/software/gsl/doc/html/vectors.html#c.gsl_matrix) structure contains six components, the two dimensions of the matrix, a physical dimension, a pointer to the memory where the elements of the matrix are stored, `data`, a pointer to the block owned by the matrix `block`, if any, and an ownership flag, `owner`. The physical dimension determines the memory layout and can differ from the matrix dimension to allow the use of submatrices. The [`gsl_matrix`](https://www.gnu.org/software/gsl/doc/html/vectors.html#c.gsl_matrix) structure is very simple and looks like this:`typedef struct {   size_t size1;   size_t size2;   size_t tda;   double * data;   gsl_block * block;   int owner; } gsl_matrix; `

Matrices are stored in row-major order, meaning that each row of elements forms a contiguous block in memory. This is the standard “C-language ordering” of two-dimensional arrays. Note that Fortran stores arrays in column-major order. The number of rows is `size1`. The range of valid row indices runs from 0 to `size1 - 1`. Similarly `size2` is the number of columns. The range of valid column indices runs from 0 to `size2 - 1`. The physical row dimension `tda`, or *trailing dimension*, specifies the size of a row of the matrix as laid out in memory.

For example, in the following matrix `size1` is 3, `size2` is 4, and `tda` is 8. The physical memory layout of the matrix begins in the top left hand-corner and proceeds from left to right along each row in turn.

```
00 01 02 03 XX XX XX XX
10 11 12 13 XX XX XX XX
20 21 22 23 XX XX XX XX
```

Each unused memory location is represented by “`XX`”. The pointer `data` gives the location of the first element of the matrix in memory. The pointer `block` stores the location of the memory block in which the elements of the matrix are located (if any). If the matrix owns this block then the `owner` field is set to one and the block will be deallocated when the matrix is freed. If the matrix is only a slice of a block owned by another object then the `owner` field is zero and any underlying block will not be freed.

The functions for allocating and accessing matrices are defined in `gsl_matrix.h`.

### Matrix allocation

The functions for allocating memory to a matrix follow the style of `malloc` and `free`. They also perform their own error checking. If there is insufficient memory available to allocate a matrix then the functions call the GSL error handler (with an error number of [`GSL_ENOMEM`](https://www.gnu.org/software/gsl/doc/html/err.html#c.GSL_ENOMEM)) in addition to returning a null pointer. Thus if you use the library error handler to abort your program then it isn’t necessary to check every `alloc`.

- [gsl_matrix](https://www.gnu.org/software/gsl/doc/html/vectors.html#c.gsl_matrix) * `gsl_matrix_alloc`(size_t *n1*, size_t *n2*)

  This function creates a matrix of size `n1` rows by `n2` columns, returning a pointer to a newly initialized matrix struct. A new block is allocated for the elements of the matrix, and stored in the `block` component of the matrix struct. The block is “owned” by the matrix, and will be deallocated when the matrix is deallocated. Requesting zero for `n1` or `n2` is valid and returns a non-null result.

- [gsl_matrix](https://www.gnu.org/software/gsl/doc/html/vectors.html#c.gsl_matrix) * `gsl_matrix_calloc`(size_t *n1*, size_t *n2*)

  This function allocates memory for a matrix of size `n1` rows by `n2` columns and initializes all the elements of the matrix to zero.

- void `gsl_matrix_free`([gsl_matrix](https://www.gnu.org/software/gsl/doc/html/vectors.html#c.gsl_matrix) * *m*)

  This function frees a previously allocated matrix `m`. If the matrix was created using [`gsl_matrix_alloc()`](https://www.gnu.org/software/gsl/doc/html/vectors.html#c.gsl_matrix_alloc) then the block underlying the matrix will also be deallocated. If the matrix has been created from another object then the memory is still owned by that object and will not be deallocated.



### Accessing matrix elements

The functions for accessing the elements of a matrix use the same range checking system as vectors. You can turn off range checking by recompiling your program with the preprocessor definition [`GSL_RANGE_CHECK_OFF`](https://www.gnu.org/software/gsl/doc/html/vectors.html#c.GSL_RANGE_CHECK_OFF).

The elements of the matrix are stored in “C-order”, where the second index moves continuously through memory. More precisely, the element accessed by the function `gsl_matrix_get(m,i,j)` and`gsl_matrix_set(m,i,j,x)` is:

```
m->data[i * m->tda + j]
```

where `tda` is the physical row-length of the matrix.

- double `gsl_matrix_get`(const [gsl_matrix](https://www.gnu.org/software/gsl/doc/html/vectors.html#c.gsl_matrix) * *m*, const size_t *i*, const size_t *j*)

  This function returns the ![(i,j)](https://www.gnu.org/software/gsl/doc/html/_images/math/4b982cbeb3e7a1b695da41fa331848138a69542a.png)-th element of a matrix `m`. If `i` or `j` lie outside the allowed range of 0 to `n1 - 1` and 0 to `n2 - 1` then the error handler is invoked and 0 is returned. An inline version of this function is used when `HAVE_INLINE` is defined.

- void `gsl_matrix_set`([gsl_matrix](https://www.gnu.org/software/gsl/doc/html/vectors.html#c.gsl_matrix) * *m*, const size_t *i*, const size_t *j*, double *x*)

  This function sets the value of the ![(i,j)](https://www.gnu.org/software/gsl/doc/html/_images/math/4b982cbeb3e7a1b695da41fa331848138a69542a.png)-th element of a matrix `m` to `x`. If `i` or `j` lies outside the allowed range of 0 to `n1 - 1` and 0 to `n2 - 1` then the error handler is invoked. An inline version of this function is used when `HAVE_INLINE` is defined.

- double * `gsl_matrix_ptr`([gsl_matrix](https://www.gnu.org/software/gsl/doc/html/vectors.html#c.gsl_matrix) * *m*, size_t *i*, size_t *j*)

- const double * `gsl_matrix_const_ptr`(const [gsl_matrix](https://www.gnu.org/software/gsl/doc/html/vectors.html#c.gsl_matrix) * *m*, size_t *i*, size_t *j*)

  These functions return a pointer to the ![(i,j)](https://www.gnu.org/software/gsl/doc/html/_images/math/4b982cbeb3e7a1b695da41fa331848138a69542a.png)-th element of a matrix `m`. If `i` or `j` lie outside the allowed range of 0 to `n1 - 1` and 0 to `n2 - 1` then the error handler is invoked and a null pointer is returned. Inline versions of these functions are used when `HAVE_INLINE` is defined.



### Initializing matrix elements

- void `gsl_matrix_set_all`([gsl_matrix](https://www.gnu.org/software/gsl/doc/html/vectors.html#c.gsl_matrix) * *m*, double *x*)

  This function sets all the elements of the matrix `m` to the value `x`.

- void `gsl_matrix_set_zero`([gsl_matrix](https://www.gnu.org/software/gsl/doc/html/vectors.html#c.gsl_matrix) * *m*)

  This function sets all the elements of the matrix `m` to zero.

- void `gsl_matrix_set_identity`([gsl_matrix](https://www.gnu.org/software/gsl/doc/html/vectors.html#c.gsl_matrix) * *m*)

  This function sets the elements of the matrix `m` to the corresponding elements of the identity matrix, ![m(i,j) = \delta(i,j)](https://www.gnu.org/software/gsl/doc/html/_images/math/e0722a7fd1ef92adb03968e5ccf3f83ecee7f9b1.png), i.e. a unit diagonal with all off-diagonal elements zero. This applies to both square and rectangular matrices.

### Reading and writing matrices

The library provides functions for reading and writing matrices to a file as binary data or formatted text.

- int `gsl_matrix_fwrite`(FILE * *stream*, const [gsl_matrix](https://www.gnu.org/software/gsl/doc/html/vectors.html#c.gsl_matrix) * *m*)

  This function writes the elements of the matrix `m` to the stream `stream` in binary format. The return value is 0 for success and `GSL_EFAILED` if there was a problem writing to the file. Since the data is written in the native binary format it may not be portable between different architectures.

- int `gsl_matrix_fread`(FILE * *stream*, [gsl_matrix](https://www.gnu.org/software/gsl/doc/html/vectors.html#c.gsl_matrix) * *m*)

  This function reads into the matrix `m` from the open stream `stream` in binary format. The matrix `m` must be preallocated with the correct dimensions since the function uses the size of `m` to determine how many bytes to read. The return value is 0 for success and `GSL_EFAILED` if there was a problem reading from the file. The data is assumed to have been written in the native binary format on the same architecture.

- int `gsl_matrix_fprintf`(FILE * *stream*, const [gsl_matrix](https://www.gnu.org/software/gsl/doc/html/vectors.html#c.gsl_matrix) * *m*, const char * *format*)

  This function writes the elements of the matrix `m` line-by-line to the stream `stream` using the format specifier `format`, which should be one of the `%g`, `%e` or `%f` formats for floating point numbers and `%d` for integers. The function returns 0 for success and `GSL_EFAILED` if there was a problem writing to the file.

- int `gsl_matrix_fscanf`(FILE * *stream*, [gsl_matrix](https://www.gnu.org/software/gsl/doc/html/vectors.html#c.gsl_matrix) * *m*)

  This function reads formatted data from the stream `stream` into the matrix `m`. The matrix `m`must be preallocated with the correct dimensions since the function uses the size of `m` to determine how many numbers to read. The function returns 0 for success and `GSL_EFAILED` if there was a problem reading from the file.

### Matrix views

- `gsl_matrix_view`

- `gsl_matrix_const_view`

  A matrix view is a temporary object, stored on the stack, which can be used to operate on a subset of matrix elements. Matrix views can be defined for both constant and non-constant matrices using separate types that preserve constness. A matrix view has the type[`gsl_matrix_view`](https://www.gnu.org/software/gsl/doc/html/vectors.html#c.gsl_matrix_view) and a constant matrix view has the type [`gsl_matrix_const_view`](https://www.gnu.org/software/gsl/doc/html/vectors.html#c.gsl_matrix_const_view). In both cases the elements of the view can by accessed using the `matrix` component of the view object. A pointer `gsl_matrix *` or `const gsl_matrix *` can be obtained by taking the address of the `matrix` component with the `&` operator. In addition to matrix views it is also possible to create vector views of a matrix, such as row or column views.

- [gsl_matrix_view](https://www.gnu.org/software/gsl/doc/html/vectors.html#c.gsl_matrix_view) `gsl_matrix_submatrix`([gsl_matrix](https://www.gnu.org/software/gsl/doc/html/vectors.html#c.gsl_matrix) * *m*, size_t *k1*, size_t *k2*, size_t *n1*, size_t *n2*)

- [gsl_matrix_const_view](https://www.gnu.org/software/gsl/doc/html/vectors.html#c.gsl_matrix_const_view) `gsl_matrix_const_submatrix`(const [gsl_matrix](https://www.gnu.org/software/gsl/doc/html/vectors.html#c.gsl_matrix) * *m*, size_t *k1*, size_t *k2*, size_t *n1*, size_t *n2*)

  These functions return a matrix view of a submatrix of the matrix `m`. The upper-left element of the submatrix is the element (`k1`, `k2`) of the original matrix. The submatrix has `n1` rows and `n2` columns. The physical number of columns in memory given by `tda` is unchanged. Mathematically, the ![(i,j)](https://www.gnu.org/software/gsl/doc/html/_images/math/4b982cbeb3e7a1b695da41fa331848138a69542a.png)-th element of the new matrix is given by:`m'(i,j) = m->data[(k1*m->tda + k2) + i*m->tda + j] `where the index `i` runs from 0 to `n1 - 1` and the index `j` runs from 0 to `n2 - 1`.The `data` pointer of the returned matrix struct is set to null if the combined parameters (`i`, `j`, `n1`, `n2`, `tda`) overrun the ends of the original matrix.The new matrix view is only a view of the block underlying the existing matrix, `m`. The block containing the elements of `m` is not owned by the new matrix view. When the view goes out of scope the original matrix `m` and its block will continue to exist. The original memory can only be deallocated by freeing the original matrix. Of course, the original matrix should not be deallocated while the view is still in use.The function [`gsl_matrix_const_submatrix()`](https://www.gnu.org/software/gsl/doc/html/vectors.html#c.gsl_matrix_const_submatrix) is equivalent to [`gsl_matrix_submatrix()`](https://www.gnu.org/software/gsl/doc/html/vectors.html#c.gsl_matrix_submatrix) but can be used for matrices which are declared `const`.

- [gsl_matrix_view](https://www.gnu.org/software/gsl/doc/html/vectors.html#c.gsl_matrix_view) `gsl_matrix_view_array`(double * *base*, size_t *n1*, size_t *n2*)

- [gsl_matrix_const_view](https://www.gnu.org/software/gsl/doc/html/vectors.html#c.gsl_matrix_const_view) `gsl_matrix_const_view_array`(const double * *base*, size_t *n1*, size_t *n2*)

  These functions return a matrix view of the array `base`. The matrix has `n1` rows and `n2`columns. The physical number of columns in memory is also given by `n2`. Mathematically, the ![(i,j)](https://www.gnu.org/software/gsl/doc/html/_images/math/4b982cbeb3e7a1b695da41fa331848138a69542a.png)-th element of the new matrix is given by:`m'(i,j) = base[i*n2 + j] `where the index `i` runs from 0 to `n1 - 1` and the index `j` runs from 0 to `n2 - 1`.The new matrix is only a view of the array `base`. When the view goes out of scope the original array `base` will continue to exist. The original memory can only be deallocated by freeing the original array. Of course, the original array should not be deallocated while the view is still in use.The function [`gsl_matrix_const_view_array()`](https://www.gnu.org/software/gsl/doc/html/vectors.html#c.gsl_matrix_const_view_array) is equivalent to [`gsl_matrix_view_array()`](https://www.gnu.org/software/gsl/doc/html/vectors.html#c.gsl_matrix_view_array) but can be used for matrices which are declared `const`.

- [gsl_matrix_view](https://www.gnu.org/software/gsl/doc/html/vectors.html#c.gsl_matrix_view) `gsl_matrix_view_array_with_tda`(double * *base*, size_t *n1*, size_t *n2*, size_t *tda*)

- [gsl_matrix_const_view](https://www.gnu.org/software/gsl/doc/html/vectors.html#c.gsl_matrix_const_view) `gsl_matrix_const_view_array_with_tda`(const double * *base*, size_t *n1*, size_t *n2*, size_t *tda*)

  These functions return a matrix view of the array `base` with a physical number of columns `tda`which may differ from the corresponding dimension of the matrix. The matrix has `n1` rows and `n2` columns, and the physical number of columns in memory is given by `tda`. Mathematically, the ![(i,j)](https://www.gnu.org/software/gsl/doc/html/_images/math/4b982cbeb3e7a1b695da41fa331848138a69542a.png)-th element of the new matrix is given by:`m'(i,j) = base[i*tda + j] `where the index `i` runs from 0 to `n1 - 1` and the index `j` runs from 0 to `n2 - 1`.The new matrix is only a view of the array `base`. When the view goes out of scope the original array `base` will continue to exist. The original memory can only be deallocated by freeing the original array. Of course, the original array should not be deallocated while the view is still in use.The function [`gsl_matrix_const_view_array_with_tda()`](https://www.gnu.org/software/gsl/doc/html/vectors.html#c.gsl_matrix_const_view_array_with_tda) is equivalent to [`gsl_matrix_view_array_with_tda()`](https://www.gnu.org/software/gsl/doc/html/vectors.html#c.gsl_matrix_view_array_with_tda) but can be used for matrices which are declared `const`.

- [gsl_matrix_view](https://www.gnu.org/software/gsl/doc/html/vectors.html#c.gsl_matrix_view) `gsl_matrix_view_vector`([gsl_vector](https://www.gnu.org/software/gsl/doc/html/vectors.html#c.gsl_vector) * *v*, size_t *n1*, size_t *n2*)

- [gsl_matrix_const_view](https://www.gnu.org/software/gsl/doc/html/vectors.html#c.gsl_matrix_const_view) `gsl_matrix_const_view_vector`(const [gsl_vector](https://www.gnu.org/software/gsl/doc/html/vectors.html#c.gsl_vector) * *v*, size_t *n1*, size_t *n2*)

  These functions return a matrix view of the vector `v`. The matrix has `n1` rows and `n2`columns. The vector must have unit stride. The physical number of columns in memory is also given by `n2`. Mathematically, the ![(i,j)](https://www.gnu.org/software/gsl/doc/html/_images/math/4b982cbeb3e7a1b695da41fa331848138a69542a.png)-th element of the new matrix is given by:`m'(i,j) = v->data[i*n2 + j] `where the index `i` runs from 0 to `n1 - 1` and the index `j` runs from 0 to `n2 - 1`.The new matrix is only a view of the vector `v`. When the view goes out of scope the original vector `v` will continue to exist. The original memory can only be deallocated by freeing the original vector. Of course, the original vector should not be deallocated while the view is still in use.The function [`gsl_matrix_const_view_vector()`](https://www.gnu.org/software/gsl/doc/html/vectors.html#c.gsl_matrix_const_view_vector) is equivalent to [`gsl_matrix_view_vector()`](https://www.gnu.org/software/gsl/doc/html/vectors.html#c.gsl_matrix_view_vector) but can be used for matrices which are declared `const`.

- [gsl_matrix_view](https://www.gnu.org/software/gsl/doc/html/vectors.html#c.gsl_matrix_view) `gsl_matrix_view_vector_with_tda`([gsl_vector](https://www.gnu.org/software/gsl/doc/html/vectors.html#c.gsl_vector) * *v*, size_t *n1*, size_t *n2*, size_t *tda*)

- [gsl_matrix_const_view](https://www.gnu.org/software/gsl/doc/html/vectors.html#c.gsl_matrix_const_view) `gsl_matrix_const_view_vector_with_tda`(const [gsl_vector](https://www.gnu.org/software/gsl/doc/html/vectors.html#c.gsl_vector) * *v*, size_t *n1*, size_t *n2*, size_t *tda*)

  These functions return a matrix view of the vector `v` with a physical number of columns `tda`which may differ from the corresponding matrix dimension. The vector must have unit stride. The matrix has `n1` rows and `n2` columns, and the physical number of columns in memory is given by `tda`. Mathematically, the ![(i,j)](https://www.gnu.org/software/gsl/doc/html/_images/math/4b982cbeb3e7a1b695da41fa331848138a69542a.png)-th element of the new matrix is given by:`m'(i,j) = v->data[i*tda + j] `where the index `i` runs from 0 to `n1 - 1` and the index `j` runs from 0 to `n2 - 1`.The new matrix is only a view of the vector `v`. When the view goes out of scope the original vector `v` will continue to exist. The original memory can only be deallocated by freeing the original vector. Of course, the original vector should not be deallocated while the view is still in use.The function [`gsl_matrix_const_view_vector_with_tda()`](https://www.gnu.org/software/gsl/doc/html/vectors.html#c.gsl_matrix_const_view_vector_with_tda) is equivalent to [`gsl_matrix_view_vector_with_tda()`](https://www.gnu.org/software/gsl/doc/html/vectors.html#c.gsl_matrix_view_vector_with_tda) but can be used for matrices which are declared `const`.

### Creating row and column views

In general there are two ways to access an object, by reference or by copying. The functions described in this section create vector views which allow access to a row or column of a matrix by reference. Modifying elements of the view is equivalent to modifying the matrix, since both the vector view and the matrix point to the same memory block.

- [gsl_vector_view](https://www.gnu.org/software/gsl/doc/html/vectors.html#c.gsl_vector_view) `gsl_matrix_row`([gsl_matrix](https://www.gnu.org/software/gsl/doc/html/vectors.html#c.gsl_matrix) * *m*, size_t *i*)

- [gsl_vector_const_view](https://www.gnu.org/software/gsl/doc/html/vectors.html#c.gsl_vector_const_view) `gsl_matrix_const_row`(const [gsl_matrix](https://www.gnu.org/software/gsl/doc/html/vectors.html#c.gsl_matrix) * *m*, size_t *i*)

  These functions return a vector view of the `i`-th row of the matrix `m`. The `data` pointer of the new vector is set to null if `i` is out of range.The function [`gsl_matrix_const_row()`](https://www.gnu.org/software/gsl/doc/html/vectors.html#c.gsl_matrix_const_row) is equivalent to [`gsl_matrix_row()`](https://www.gnu.org/software/gsl/doc/html/vectors.html#c.gsl_matrix_row) but can be used for matrices which are declared `const`.

- [gsl_vector_view](https://www.gnu.org/software/gsl/doc/html/vectors.html#c.gsl_vector_view) `gsl_matrix_column`([gsl_matrix](https://www.gnu.org/software/gsl/doc/html/vectors.html#c.gsl_matrix) * *m*, size_t *j*)

- [gsl_vector_const_view](https://www.gnu.org/software/gsl/doc/html/vectors.html#c.gsl_vector_const_view) `gsl_matrix_const_column`(const [gsl_matrix](https://www.gnu.org/software/gsl/doc/html/vectors.html#c.gsl_matrix) * *m*, size_t *j*)

  These functions return a vector view of the `j`-th column of the matrix `m`. The `data` pointer of the new vector is set to null if `j` is out of range.The function [`gsl_matrix_const_column()`](https://www.gnu.org/software/gsl/doc/html/vectors.html#c.gsl_matrix_const_column) is equivalent to [`gsl_matrix_column()`](https://www.gnu.org/software/gsl/doc/html/vectors.html#c.gsl_matrix_column) but can be used for matrices which are declared `const`.

- [gsl_vector_view](https://www.gnu.org/software/gsl/doc/html/vectors.html#c.gsl_vector_view) `gsl_matrix_subrow`([gsl_matrix](https://www.gnu.org/software/gsl/doc/html/vectors.html#c.gsl_matrix) * *m*, size_t *i*, size_t *offset*, size_t *n*)

- [gsl_vector_const_view](https://www.gnu.org/software/gsl/doc/html/vectors.html#c.gsl_vector_const_view) `gsl_matrix_const_subrow`(const [gsl_matrix](https://www.gnu.org/software/gsl/doc/html/vectors.html#c.gsl_matrix) * *m*, size_t *i*, size_t *offset*, size_t *n*)

  These functions return a vector view of the `i`-th row of the matrix `m` beginning at `offset`elements past the first column and containing `n` elements. The `data` pointer of the new vector is set to null if `i`, `offset`, or `n` are out of range.The function [`gsl_matrix_const_subrow()`](https://www.gnu.org/software/gsl/doc/html/vectors.html#c.gsl_matrix_const_subrow) is equivalent to [`gsl_matrix_subrow()`](https://www.gnu.org/software/gsl/doc/html/vectors.html#c.gsl_matrix_subrow) but can be used for matrices which are declared `const`.

- [gsl_vector_view](https://www.gnu.org/software/gsl/doc/html/vectors.html#c.gsl_vector_view) `gsl_matrix_subcolumn`([gsl_matrix](https://www.gnu.org/software/gsl/doc/html/vectors.html#c.gsl_matrix) * *m*, size_t *j*, size_t *offset*, size_t *n*)

- [gsl_vector_const_view](https://www.gnu.org/software/gsl/doc/html/vectors.html#c.gsl_vector_const_view) `gsl_matrix_const_subcolumn`(const [gsl_matrix](https://www.gnu.org/software/gsl/doc/html/vectors.html#c.gsl_matrix) * *m*, size_t *j*, size_t *offset*, size_t *n*)

  These functions return a vector view of the `j`-th column of the matrix `m` beginning at `offset`elements past the first row and containing `n` elements. The `data` pointer of the new vector is set to null if `j`, `offset`, or `n` are out of range.The function [`gsl_matrix_const_subcolumn()`](https://www.gnu.org/software/gsl/doc/html/vectors.html#c.gsl_matrix_const_subcolumn) is equivalent to [`gsl_matrix_subcolumn()`](https://www.gnu.org/software/gsl/doc/html/vectors.html#c.gsl_matrix_subcolumn) but can be used for matrices which are declared `const`.



- [gsl_vector_view](https://www.gnu.org/software/gsl/doc/html/vectors.html#c.gsl_vector_view) `gsl_matrix_diagonal`([gsl_matrix](https://www.gnu.org/software/gsl/doc/html/vectors.html#c.gsl_matrix) * *m*)

- [gsl_vector_const_view](https://www.gnu.org/software/gsl/doc/html/vectors.html#c.gsl_vector_const_view) `gsl_matrix_const_diagonal`(const [gsl_matrix](https://www.gnu.org/software/gsl/doc/html/vectors.html#c.gsl_matrix) * *m*)

  These functions return a vector view of the diagonal of the matrix `m`. The matrix `m` is not required to be square. For a rectangular matrix the length of the diagonal is the same as the smaller dimension of the matrix.The function [`gsl_matrix_const_diagonal()`](https://www.gnu.org/software/gsl/doc/html/vectors.html#c.gsl_matrix_const_diagonal) is equivalent to [`gsl_matrix_diagonal()`](https://www.gnu.org/software/gsl/doc/html/vectors.html#c.gsl_matrix_diagonal) but can be used for matrices which are declared `const`.



- [gsl_vector_view](https://www.gnu.org/software/gsl/doc/html/vectors.html#c.gsl_vector_view) `gsl_matrix_subdiagonal`([gsl_matrix](https://www.gnu.org/software/gsl/doc/html/vectors.html#c.gsl_matrix) * *m*, size_t *k*)

- [gsl_vector_const_view](https://www.gnu.org/software/gsl/doc/html/vectors.html#c.gsl_vector_const_view) `gsl_matrix_const_subdiagonal`(const [gsl_matrix](https://www.gnu.org/software/gsl/doc/html/vectors.html#c.gsl_matrix) * *m*, size_t *k*)

  These functions return a vector view of the `k`-th subdiagonal of the matrix `m`. The matrix `m`is not required to be square. The diagonal of the matrix corresponds to ![k = 0](https://www.gnu.org/software/gsl/doc/html/_images/math/7fe0629d90b771b38027bce3cadc1a36afe514bb.png).The function [`gsl_matrix_const_subdiagonal()`](https://www.gnu.org/software/gsl/doc/html/vectors.html#c.gsl_matrix_const_subdiagonal) is equivalent to [`gsl_matrix_subdiagonal()`](https://www.gnu.org/software/gsl/doc/html/vectors.html#c.gsl_matrix_subdiagonal) but can be used for matrices which are declared `const`.



- [gsl_vector_view](https://www.gnu.org/software/gsl/doc/html/vectors.html#c.gsl_vector_view) `gsl_matrix_superdiagonal`([gsl_matrix](https://www.gnu.org/software/gsl/doc/html/vectors.html#c.gsl_matrix) * *m*, size_t *k*)

- [gsl_vector_const_view](https://www.gnu.org/software/gsl/doc/html/vectors.html#c.gsl_vector_const_view) `gsl_matrix_const_superdiagonal`(const [gsl_matrix](https://www.gnu.org/software/gsl/doc/html/vectors.html#c.gsl_matrix) * *m*, size_t *k*)

  These functions return a vector view of the `k`-th superdiagonal of the matrix `m`. The matrix `m` is not required to be square. The diagonal of the matrix corresponds to ![k = 0](https://www.gnu.org/software/gsl/doc/html/_images/math/7fe0629d90b771b38027bce3cadc1a36afe514bb.png).The function [`gsl_matrix_const_superdiagonal()`](https://www.gnu.org/software/gsl/doc/html/vectors.html#c.gsl_matrix_const_superdiagonal) is equivalent to [`gsl_matrix_superdiagonal()`](https://www.gnu.org/software/gsl/doc/html/vectors.html#c.gsl_matrix_superdiagonal)but can be used for matrices which are declared `const`.

### Copying matrices

- int `gsl_matrix_memcpy`([gsl_matrix](https://www.gnu.org/software/gsl/doc/html/vectors.html#c.gsl_matrix) * *dest*, const [gsl_matrix](https://www.gnu.org/software/gsl/doc/html/vectors.html#c.gsl_matrix) * *src*)

  This function copies the elements of the matrix `src` into the matrix `dest`. The two matrices must have the same size.

- int `gsl_matrix_swap`([gsl_matrix](https://www.gnu.org/software/gsl/doc/html/vectors.html#c.gsl_matrix) * *m1*, [gsl_matrix](https://www.gnu.org/software/gsl/doc/html/vectors.html#c.gsl_matrix) * *m2*)

  This function exchanges the elements of the matrices `m1` and `m2` by copying. The two matrices must have the same size.

### Copying rows and columns

The functions described in this section copy a row or column of a matrix into a vector. This allows the elements of the vector and the matrix to be modified independently. Note that if the matrix and the vector point to overlapping regions of memory then the result will be undefined. The same effect can be achieved with more generality using [`gsl_vector_memcpy()`](https://www.gnu.org/software/gsl/doc/html/vectors.html#c.gsl_vector_memcpy) with vector views of rows and columns.

- int `gsl_matrix_get_row`([gsl_vector](https://www.gnu.org/software/gsl/doc/html/vectors.html#c.gsl_vector) * *v*, const [gsl_matrix](https://www.gnu.org/software/gsl/doc/html/vectors.html#c.gsl_matrix) * *m*, size_t *i*)

  This function copies the elements of the `i`-th row of the matrix `m` into the vector `v`. The length of the vector must be the same as the length of the row.

- int `gsl_matrix_get_col`([gsl_vector](https://www.gnu.org/software/gsl/doc/html/vectors.html#c.gsl_vector) * *v*, const [gsl_matrix](https://www.gnu.org/software/gsl/doc/html/vectors.html#c.gsl_matrix) * *m*, size_t *j*)

  This function copies the elements of the `j`-th column of the matrix `m` into the vector `v`. The length of the vector must be the same as the length of the column.

- int `gsl_matrix_set_row`([gsl_matrix](https://www.gnu.org/software/gsl/doc/html/vectors.html#c.gsl_matrix) * *m*, size_t *i*, const [gsl_vector](https://www.gnu.org/software/gsl/doc/html/vectors.html#c.gsl_vector) * *v*)

  This function copies the elements of the vector `v` into the `i`-th row of the matrix `m`. The length of the vector must be the same as the length of the row.

- int `gsl_matrix_set_col`([gsl_matrix](https://www.gnu.org/software/gsl/doc/html/vectors.html#c.gsl_matrix) * *m*, size_t *j*, const [gsl_vector](https://www.gnu.org/software/gsl/doc/html/vectors.html#c.gsl_vector) * *v*)

  This function copies the elements of the vector `v` into the `j`-th column of the matrix `m`. The length of the vector must be the same as the length of the column.

### Exchanging rows and columns

The following functions can be used to exchange the rows and columns of a matrix.

- int `gsl_matrix_swap_rows`([gsl_matrix](https://www.gnu.org/software/gsl/doc/html/vectors.html#c.gsl_matrix) * *m*, size_t *i*, size_t *j*)

  This function exchanges the `i`-th and `j`-th rows of the matrix `m` in-place.

- int `gsl_matrix_swap_columns`([gsl_matrix](https://www.gnu.org/software/gsl/doc/html/vectors.html#c.gsl_matrix) * *m*, size_t *i*, size_t *j*)

  This function exchanges the `i`-th and `j`-th columns of the matrix `m` in-place.

- int `gsl_matrix_swap_rowcol`([gsl_matrix](https://www.gnu.org/software/gsl/doc/html/vectors.html#c.gsl_matrix) * *m*, size_t *i*, size_t *j*)

  This function exchanges the `i`-th row and `j`-th column of the matrix `m` in-place. The matrix must be square for this operation to be possible.

- int `gsl_matrix_transpose_memcpy`([gsl_matrix](https://www.gnu.org/software/gsl/doc/html/vectors.html#c.gsl_matrix) * *dest*, const [gsl_matrix](https://www.gnu.org/software/gsl/doc/html/vectors.html#c.gsl_matrix) * *src*)

  This function makes the matrix `dest` the transpose of the matrix `src` by copying the elements of `src` into `dest`. This function works for all matrices provided that the dimensions of the matrix `dest` match the transposed dimensions of the matrix `src`.

- int `gsl_matrix_transpose`([gsl_matrix](https://www.gnu.org/software/gsl/doc/html/vectors.html#c.gsl_matrix) * *m*)

  This function replaces the matrix `m` by its transpose by copying the elements of the matrix in-place. The matrix must be square for this operation to be possible.

### Matrix operations

The following operations are defined for real and complex matrices.

- int `gsl_matrix_add`([gsl_matrix](https://www.gnu.org/software/gsl/doc/html/vectors.html#c.gsl_matrix) * *a*, const [gsl_matrix](https://www.gnu.org/software/gsl/doc/html/vectors.html#c.gsl_matrix) * *b*)

  This function adds the elements of matrix `b` to the elements of matrix `a`. The result ![a(i,j) \leftarrow a(i,j) + b(i,j)](https://www.gnu.org/software/gsl/doc/html/_images/math/0bc02b7c27b12d89af55c6542b82dc6a9fff1341.png) is stored in `a` and `b` remains unchanged. The two matrices must have the same dimensions.

- int `gsl_matrix_sub`([gsl_matrix](https://www.gnu.org/software/gsl/doc/html/vectors.html#c.gsl_matrix) * *a*, const [gsl_matrix](https://www.gnu.org/software/gsl/doc/html/vectors.html#c.gsl_matrix) * *b*)

  This function subtracts the elements of matrix `b` from the elements of matrix `a`. The result ![a(i,j) \leftarrow a(i,j) - b(i,j)](https://www.gnu.org/software/gsl/doc/html/_images/math/46aaa0f18b648e3943f34e301abfe7fc84f031db.png) is stored in `a` and `b` remains unchanged. The two matrices must have the same dimensions.

- int `gsl_matrix_mul_elements`([gsl_matrix](https://www.gnu.org/software/gsl/doc/html/vectors.html#c.gsl_matrix) * *a*, const [gsl_matrix](https://www.gnu.org/software/gsl/doc/html/vectors.html#c.gsl_matrix) * *b*)

  This function multiplies the elements of matrix `a` by the elements of matrix `b`. The result ![a(i,j) \leftarrow a(i,j) * b(i,j)](https://www.gnu.org/software/gsl/doc/html/_images/math/9ac89ff82a8fd8539af10562d5267550fb8405f6.png) is stored in `a` and `b` remains unchanged. The two matrices must have the same dimensions.

- int `gsl_matrix_div_elements`([gsl_matrix](https://www.gnu.org/software/gsl/doc/html/vectors.html#c.gsl_matrix) * *a*, const [gsl_matrix](https://www.gnu.org/software/gsl/doc/html/vectors.html#c.gsl_matrix) * *b*)

  This function divides the elements of matrix `a` by the elements of matrix `b`. The result ![a(i,j) \leftarrow a(i,j) / b(i,j)](https://www.gnu.org/software/gsl/doc/html/_images/math/6db76719f396f5103700e2e921bab5ba7afc584d.png) is stored in `a` and `b` remains unchanged. The two matrices must have the same dimensions.

- int `gsl_matrix_scale`([gsl_matrix](https://www.gnu.org/software/gsl/doc/html/vectors.html#c.gsl_matrix) * *a*, const double *x*)

  This function multiplies the elements of matrix `a` by the constant factor `x`. The result ![a(i,j) \leftarrow x a(i,j)](https://www.gnu.org/software/gsl/doc/html/_images/math/dfec17a539ac118278e38e5e98a8182edd587a25.png) is stored in `a`.

- int `gsl_matrix_add_constant`([gsl_matrix](https://www.gnu.org/software/gsl/doc/html/vectors.html#c.gsl_matrix) * *a*, const double *x*)

  This function adds the constant value `x` to the elements of the matrix `a`. The result ![a(i,j) \leftarrow a(i,j) + x](https://www.gnu.org/software/gsl/doc/html/_images/math/8ffa3ca8730fd4230aafdba16ec24c5f54140099.png) is stored in `a`.

### Finding maximum and minimum elements of matrices

The following operations are only defined for real matrices.

- double `gsl_matrix_max`(const [gsl_matrix](https://www.gnu.org/software/gsl/doc/html/vectors.html#c.gsl_matrix) * *m*)

  This function returns the maximum value in the matrix `m`.

- double `gsl_matrix_min`(const [gsl_matrix](https://www.gnu.org/software/gsl/doc/html/vectors.html#c.gsl_matrix) * *m*)

  This function returns the minimum value in the matrix `m`.

- void `gsl_matrix_minmax`(const [gsl_matrix](https://www.gnu.org/software/gsl/doc/html/vectors.html#c.gsl_matrix) * *m*, double * *min_out*, double * *max_out*)

  This function returns the minimum and maximum values in the matrix `m`, storing them in `min_out` and `max_out`.

- void `gsl_matrix_max_index`(const [gsl_matrix](https://www.gnu.org/software/gsl/doc/html/vectors.html#c.gsl_matrix) * *m*, size_t * *imax*, size_t * *jmax*)

  This function returns the indices of the maximum value in the matrix `m`, storing them in `imax`and `jmax`. When there are several equal maximum elements then the first element found is returned, searching in row-major order.

- void `gsl_matrix_min_index`(const [gsl_matrix](https://www.gnu.org/software/gsl/doc/html/vectors.html#c.gsl_matrix) * *m*, size_t * *imin*, size_t * *jmin*)

  This function returns the indices of the minimum value in the matrix `m`, storing them in `imin`and `jmin`. When there are several equal minimum elements then the first element found is returned, searching in row-major order.

- void `gsl_matrix_minmax_index`(const [gsl_matrix](https://www.gnu.org/software/gsl/doc/html/vectors.html#c.gsl_matrix) * *m*, size_t * *imin*, size_t * *jmin*, size_t * *imax*, size_t * *jmax*)

  This function returns the indices of the minimum and maximum values in the matrix `m`, storing them in (`imin`, `jmin`) and (`imax`, `jmax`). When there are several equal minimum or maximum elements then the first elements found are returned, searching in row-major order.

### Matrix properties

The following functions are defined for real and complex matrices. For complex matrices both the real and imaginary parts must satisfy the conditions.

- int `gsl_matrix_isnull`(const [gsl_matrix](https://www.gnu.org/software/gsl/doc/html/vectors.html#c.gsl_matrix) * *m*)

- int `gsl_matrix_ispos`(const [gsl_matrix](https://www.gnu.org/software/gsl/doc/html/vectors.html#c.gsl_matrix) * *m*)

- int `gsl_matrix_isneg`(const [gsl_matrix](https://www.gnu.org/software/gsl/doc/html/vectors.html#c.gsl_matrix) * *m*)

- int `gsl_matrix_isnonneg`(const [gsl_matrix](https://www.gnu.org/software/gsl/doc/html/vectors.html#c.gsl_matrix) * *m*)

  These functions return 1 if all the elements of the matrix `m` are zero, strictly positive, strictly negative, or non-negative respectively, and 0 otherwise. To test whether a matrix is positive-definite, use the [Cholesky decomposition](https://www.gnu.org/software/gsl/doc/html/linalg.html#sec-cholesky-decomposition).

- int `gsl_matrix_equal`(const [gsl_matrix](https://www.gnu.org/software/gsl/doc/html/vectors.html#c.gsl_matrix) * *a*, const [gsl_matrix](https://www.gnu.org/software/gsl/doc/html/vectors.html#c.gsl_matrix) * *b*)

  This function returns 1 if the matrices `a` and `b` are equal (by comparison of element values) and 0 otherwise.

### Example programs for matrices

The program below shows how to allocate, initialize and read from a matrix using the functions [`gsl_matrix_alloc()`](https://www.gnu.org/software/gsl/doc/html/vectors.html#c.gsl_matrix_alloc), [`gsl_matrix_set()`](https://www.gnu.org/software/gsl/doc/html/vectors.html#c.gsl_matrix_set) and [`gsl_matrix_get()`](https://www.gnu.org/software/gsl/doc/html/vectors.html#c.gsl_matrix_get).

```
#include <stdio.h>
#include <gsl/gsl_matrix.h>

int
main (void)
{
  int i, j;
  gsl_matrix * m = gsl_matrix_alloc (10, 3);

  for (i = 0; i < 10; i++)
    for (j = 0; j < 3; j++)
      gsl_matrix_set (m, i, j, 0.23 + 100*i + j);

  for (i = 0; i < 100; i++)  /* OUT OF RANGE ERROR */
    for (j = 0; j < 3; j++)
      printf ("m(%d,%d) = %g\n", i, j,
              gsl_matrix_get (m, i, j));

  gsl_matrix_free (m);

  return 0;
}
```

Here is the output from the program. The final loop attempts to read outside the range of the matrix `m`, and the error is trapped by the range-checking code in [`gsl_matrix_get()`](https://www.gnu.org/software/gsl/doc/html/vectors.html#c.gsl_matrix_get).

```
$ ./a.out
m(0,0) = 0.23
m(0,1) = 1.23
m(0,2) = 2.23
m(1,0) = 100.23
m(1,1) = 101.23
m(1,2) = 102.23
...
m(9,2) = 902.23
gsl: matrix_source.c:13: ERROR: first index out of range
Default GSL error handler invoked.
Aborted (core dumped)
```

The next program shows how to write a matrix to a file.

```
#include <stdio.h>
#include <gsl/gsl_matrix.h>

int
main (void)
{
  int i, j, k = 0;
  gsl_matrix * m = gsl_matrix_alloc (100, 100);
  gsl_matrix * a = gsl_matrix_alloc (100, 100);

  for (i = 0; i < 100; i++)
    for (j = 0; j < 100; j++)
      gsl_matrix_set (m, i, j, 0.23 + i + j);

  {
     FILE * f = fopen ("test.dat", "wb");
     gsl_matrix_fwrite (f, m);
     fclose (f);
  }

  {
     FILE * f = fopen ("test.dat", "rb");
     gsl_matrix_fread (f, a);
     fclose (f);
  }

  for (i = 0; i < 100; i++)
    for (j = 0; j < 100; j++)
      {
        double mij = gsl_matrix_get (m, i, j);
        double aij = gsl_matrix_get (a, i, j);
        if (mij != aij) k++;
      }

  gsl_matrix_free (m);
  gsl_matrix_free (a);

  printf ("differences = %d (should be zero)\n", k);
  return (k > 0);
}
```

After running this program the file `test.dat` should contain the elements of `m`, written in binary format. The matrix which is read back in using the function [`gsl_matrix_fread()`](https://www.gnu.org/software/gsl/doc/html/vectors.html#c.gsl_matrix_fread) should be exactly equal to the original matrix.

The following program demonstrates the use of vector views. The program computes the column norms of a matrix.

```
#include <math.h>
#include <stdio.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_blas.h>

int
main (void)
{
  size_t i,j;

  gsl_matrix *m = gsl_matrix_alloc (10, 10);

  for (i = 0; i < 10; i++)
    for (j = 0; j < 10; j++)
      gsl_matrix_set (m, i, j, sin (i) + cos (j));

  for (j = 0; j < 10; j++)
    {
      gsl_vector_view column = gsl_matrix_column (m, j);
      double d;

      d = gsl_blas_dnrm2 (&column.vector);

      printf ("matrix column %zu, norm = %g\n", j, d);
    }

  gsl_matrix_free (m);

  return 0;
}
```

Here is the output of the program,

```
matrix column 0, norm = 4.31461
matrix column 1, norm = 3.1205
matrix column 2, norm = 2.19316
matrix column 3, norm = 3.26114
matrix column 4, norm = 2.53416
matrix column 5, norm = 2.57281
matrix column 6, norm = 4.20469
matrix column 7, norm = 3.65202
matrix column 8, norm = 2.08524
matrix column 9, norm = 3.07313
```

The results can be confirmed using GNU octave:

```
$ octave
GNU Octave, version 2.0.16.92
octave> m = sin(0:9)' * ones(1,10)
               + ones(10,1) * cos(0:9);
octave> sqrt(sum(m.^2))
ans =
  4.3146  3.1205  2.1932  3.2611  2.5342  2.5728
  4.2047  3.6520  2.0852  3.0731
```

### References and Further Reading

The block, vector and matrix objects in GSL follow the `valarray` model of C++. A description of this model can be found in the following reference,

- B. Stroustrup, The C++ Programming Language (3rd Ed), Section 22.4 Vector Arithmetic. Addison-Wesley 1997, ISBN 0-201-88954-4.

### Footnotes

[[1]](https://www.gnu.org/software/gsl/doc/html/vectors.html#id1)Range checking is available in the GNU C Compiler bounds-checking extension, but it is not part of the default installation of GCC. Memory accesses can also be checked with Valgrind or the `gcc -fmudflap`memory protection option.

