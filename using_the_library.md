## 使用本库

This chapter describes how to compile programs that use GSL, and introduces its conventions.

本章描述了如何编译使用了GSL的程序，以及介绍其惯用手法。

## 一个示例程序

The following short program demonstrates the use of the library by computing the value of the Bessel function ![J_0(x)](https://www.gnu.org/software/gsl/doc/html/_images/math/436b1799db46b1ecb38cc2fde801e882814c7119.png) for ![x=5](https://www.gnu.org/software/gsl/doc/html/_images/math/15d744c081cd9f683f1a4d673569418803e0a99b.png):

下面的小段程序通过计算Bessel函数$J_0(x)​$在$x=5​$的值来展示本库的使用:

下面的小段程序通过计算Bessel函数![J_0(x)](https://www.gnu.org/software/gsl/doc/html/_images/math/436b1799db46b1ecb38cc2fde801e882814c7119.png)在![x=5](https://www.gnu.org/software/gsl/doc/html/_images/math/15d744c081cd9f683f1a4d673569418803e0a99b.png)的值来展示本库的使用:

```
#include <stdio.h>
#include <gsl/gsl_sf_bessel.h>

int
main (void)
{
  double x = 5.0;
  double y = gsl_sf_bessel_J0 (x);
  printf ("J0(%g) = %.18e\n", x, y);
  return 0;
}
```

The output is shown below, and should be correct to double-precision accuracy [[1]](https://www.gnu.org/software/gsl/doc/html/usage.html#f1),

输出如下所示，并且对于双精度情况是正确的[[1]](https://www.gnu.org/software/gsl/doc/html/usage.html#f1)，

```
J0(5) = -1.775967713143382642e-01
```

The steps needed to compile this program are described in the following sections.

需要用来编译本程序的步骤在下面部分中得到描述。



## 编译和链接

The library header files are installed in their own `gsl` directory. You should write any preprocessor include statements with a `gsl/` directory prefix thus:

```
#include <gsl/gsl_math.h>
```

本库的头文件安装在其自身的`gsl`目录下。因此你应该编写任何具有`gsl/`目录前缀的预处理器include声明:

```
#include <gsl/gsl_math.h>
```

If the directory is not installed on the standard search path of your compiler you will also need to provide its location to the preprocessor as a command line flag. The default location of the `gsl`directory is `/usr/local/include/gsl`. A typical compilation command for a source file `example.c`with the GNU C compiler `gcc` is:

```
$ gcc -Wall -I/usr/local/include -c example.c
```

This results in an object file `example.o`. The default include path for `gcc` searches `/usr/local/include` automatically so the `-I` option can actually be omitted when GSL is installed in its default location.

如果头文件目录并非安装在你的编译器的标准搜索路径，你将需要通过命令行标记向预处理器提供它的位置。默认的`gsl`目录位置是`/usr/local/include/gsl`。一个典型的为一个源文件`example.c`所用的使用GNU C编译器`gcc`的编译命令是:

```
$ gcc -Wall -I/usr/local/include -c example.c
```

结果是一个目标文件`example.o`。供`gcc`使用的默认的include路径会自动搜索`/usr/local/include`，因此`-I` 选项实际上能够被忽略掉如果GSL是安装在默认位置。

### 将程序与库链接在一起

The library is installed as a single file, `libgsl.a`. A shared version of the library `libgsl.so` is also installed on systems that support shared libraries. The default location of these files is`/usr/local/lib`. If this directory is not on the standard search path of your linker you will also need to provide its location as a command line flag.

库是作为一个单一文件进行安装，`libgsl.a`。一个共享库版本`libgsl.so`也被安装在支持共享库的系统上。这些文件的默认位置是`/usr/local/lib`。如果这个目录没有在你的链接器的标准搜索路径上，你也应该需要以一个命令行参数来提供其位置。

To link against the library you need to specify both the main library and a supporting CBLAS library, which provides standard basic linear algebra subroutines. A suitable CBLAS implementation is provided in the library `libgslcblas.a` if your system does not provide one. The following example shows how to link an application with the library:

```
$ gcc -L/usr/local/lib example.o -lgsl -lgslcblas -lm
```

为了链接库你需要明确主库以及提供标准线性代数子程序支持的CBLAS库。如果你的系统没有提供一个合适的CBLAS实现，在库`libgslclbas.a`中进行了提供。以下例子展示了如何将库链接到应用:

```
$ gcc -L/usr/local/lib example.o -lgsl -lgslcblas -lm
```

The default library path for `gcc` searches `/usr/local/lib` automatically so the `-L` option can be omitted when GSL is installed in its default location.

The option `-lm` links with the system math library. On some systems it is not needed. [[2]](https://www.gnu.org/software/gsl/doc/html/usage.html#f2)

For a tutorial introduction to the GNU C Compiler and related programs, see “An Introduction to GCC” (ISBN 0954161793). [[3]](https://www.gnu.org/software/gsl/doc/html/usage.html#f3)

`gcc`的默认自动库搜索路径是`/usr/local/lib`，因此当GSL安装在默认路径时`-L`选项可以被忽略掉。

选项`-lm`与系统数学库进行链接。在一些系统上不需要指定。[[2]](https://www.gnu.org/software/gsl/doc/html/usage.html#f2)



### 链接到一个可选的BLAS库

The following command line shows how you would link the same application with an alternative CBLAS library `libcblas.a`:

```
$ gcc example.o -lgsl -lcblas -lm
```

以下命令展示了你应当如何将同一个应用链接到一个可替代的CBLAS库`libcblas.a`:

```
$ gcc example.o -lgsl -lcblas -lm
```

For the best performance an optimized platform-specific CBLAS library should be used for `-lcblas`. The library must conform to the CBLAS standard. The ATLAS package provides a portable high-performance BLAS library with a CBLAS interface. It is free software and should be installed for any work requiring fast vector and matrix operations. The following command line will link with the ATLAS library and its CBLAS interface:

```
$ gcc example.o -lgsl -lcblas -latlas -lm
```

为了最好的性能一个优化的平台特定的CBLAS库应该被通过`-lcblas`使用。库必须要符合CBLAS标准。ATLAS包提供了一个具有CBLAS接口的可移植高性能BLAS库。它是免费软件并应该被任何需要快速向量和矩阵操作的工作所安装。以下命令行将会链接到ATLAS库以及其CBLAS接口:

```
$ gcc example.o -lgsl -lcblas -latlas -lm
```

If the ATLAS library is installed in a non-standard directory use the `-L` option to add it to the search path, as described above.

For more information about BLAS functions see [BLAS Support](https://www.gnu.org/software/gsl/doc/html/blas.html#chap-blas-support).

如果ATLAS库被安装在一个非标准目录，使用`-L`选项来将其添加到搜索路径，如上所示。

了解关于BLAS函数的更多信息，参见[BLAS Support](https://www.gnu.org/software/gsl/doc/html/blas.html#chap-blas-support)。



## 共享库

To run a program linked with the shared version of the library the operating system must be able to locate the corresponding `.so` file at runtime. If the library cannot be found, the following error will occur:

```
$ ./a.out
./a.out: error while loading shared libraries:
libgsl.so.0: cannot open shared object file: No such file or directory
```

为了运行一个链接到共享库版本的程序，操作系统需要能够在运行时定位相应的`.so`文件。如果库无法被找到，下面的错误将会发生:

```
$ ./a.out
./a.out: error while loading shared libraries:
libgsl.so.0: cannot open shared object file: No such file or directory
```

To avoid this error, either modify the system dynamic linker configuration [[4]](https://www.gnu.org/software/gsl/doc/html/usage.html#f4) or define the shell variable `LD_LIBRARY_PATH` to include the directory where the library is installed.

For example, in the Bourne shell (`/bin/sh` or `/bin/bash`), the library search path can be set with the following commands:

```
$ LD_LIBRARY_PATH=/usr/local/lib
$ export LD_LIBRARY_PATH
$ ./example
```

为了避免这个错误，要么更改系统动态链接器配置[[4]](https://www.gnu.org/software/gsl/doc/html/usage.html#f4) 或者定义shell环境变量`LD_LIBRARY_PATH`以包含本库安装的目录。

例如，在Bourne shell(`/bin/sh`或者`/bin/bash`)中，库搜索路径能通过下列命令设置:

```
$ LD_LIBRARY_PATH=/usr/local/lib
$ export LD_LIBRARY_PATH
$ ./example
```

In the C-shell (`/bin/csh` or `/bin/tcsh`) the equivalent command is:

```
% setenv LD_LIBRARY_PATH /usr/local/lib
```

在C-shell(`/bin/csh`或者`/bin/tcsh`)中等价的命令是:

```
% setenv LD_LIBRARY_PATH /usr/local/lib
```

The standard prompt for the C-shell in the example above is the percent character %, and should not be typed as part of the command.

To save retyping these commands each session they can be placed in an individual or system-wide login file.

在上面例子中提到的C-shell的标准提示符是百分符号%，并不应当被当做命令一起输入。

To compile a statically linked version of the program, use the `-static` flag in `gcc`:

```
$ gcc -static example.o -lgsl -lgslcblas -lm
```

为了编译一个静态链接版本的程序，在`gcc`中使用`-static`标签:

## 符合ANSI C

The library is written in ANSI C and is intended to conform to the ANSI C standard (C89). It should be portable to any system with a working ANSI C compiler.

本库是以ANSI C进行编写的并被打算来遵守ANSI C标准(C89)。其应当能够一直到任何具有可工作的ANSI C编译器的系统上。

The library does not rely on any non-ANSI extensions in the interface it exports to the user. Programs you write using GSL can be ANSI compliant. Extensions which can be used in a way compatible with pure ANSI C are supported, however, via conditional compilation. This allows the library to take advantage of compiler extensions on those platforms which support them.

本库所导出给用户的接口中不会依赖于任何非ANSI扩展。你使用GSL编写的程序会是遵从ANSI的。能够被以与纯ANSI C兼容的方式使用的扩展也是得到了支持，然而，需要经过条件编译。这允许本库去利用那些支持它们的平台上的编译器扩展的优势。

When an ANSI C feature is known to be broken on a particular system the library will exclude any related functions at compile-time. This should make it impossible to link a program that would use these functions and give incorrect results.

当一个ANSI C特征被发现在一个特定系统上崩溃是，本库将会在编译时排除任何相关的函数。这应当能使得链接一个将要使用这些函数并给出错误结果的程序是不可能的。

To avoid namespace conflicts all exported function names and variables have the prefix `gsl_`, while exported macros have the prefix `GSL_`.

为了避免名字空间冲突，所有导出的函数名字和变量具有前缀`gsl_`，而导出的宏具有前缀`GSL_`。



## 内联函数

The `inline` keyword is not part of the original ANSI C standard (C89) so the library does not export any inline function definitions by default. Inline functions were introduced officially in the newer C99 standard but most C89 compilers have also included `inline` as an extension for a long time.

To allow the use of inline functions, the library provides optional inline versions of performance-critical routines by conditional compilation in the exported header files. The inline versions of these functions can be included by defining the macro `HAVE_INLINE` when compiling an application:

```
$ gcc -Wall -c -DHAVE_INLINE example.c
```

If you use `autoconf` this macro can be defined automatically. If you do not define the macro `HAVE_INLINE` then the slower non-inlined versions of the functions will be used instead.

By default, the actual form of the inline keyword is `extern inline`, which is a `gcc` extension that eliminates unnecessary function definitions. If the form `extern inline` causes problems with other compilers a stricter autoconf test can be used, see [Autoconf Macros](https://www.gnu.org/software/gsl/doc/html/autoconf.html#chap-autoconf-macros).

When compiling with `gcc` in C99 mode (`gcc -std=c99`) the header files automatically switch to C99-compatible inline function declarations instead of `extern inline`. With other C99 compilers, define the macro `GSL_C99_INLINE` to use these declarations.

## Long double

In general, the algorithms in the library are written for double precision only. The `long double` type is not supported for actual computation.

One reason for this choice is that the precision of `long double` is platform dependent. The IEEE standard only specifies the minimum precision of extended precision numbers, while the precision of `double` is the same on all platforms.

However, it is sometimes necessary to interact with external data in long-double format, so the vector and matrix datatypes include long-double versions.

It should be noted that in some system libraries the `stdio.h` formatted input/output functions `printf` and `scanf` are not implemented correctly for `long double`. Undefined or incorrect results are avoided by testing these functions during the `configure` stage of library compilation and eliminating certain GSL functions which depend on them if necessary. The corresponding line in the `configure` output looks like this:

```
checking whether printf works with long double... no
```

Consequently when `long double` formatted input/output does not work on a given system it should be impossible to link a program which uses GSL functions dependent on this.

If it is necessary to work on a system which does not support formatted `long double` input/output then the options are to use binary formats or to convert `long double` results into `double` for reading and writing.



## Portability functions

To help in writing portable applications GSL provides some implementations of functions that are found in other libraries, such as the BSD math library. You can write your application to use the native versions of these functions, and substitute the GSL versions via a preprocessor macro if they are unavailable on another platform.

For example, after determining whether the BSD function `hypot()` is available you can include the following macro definitions in a file `config.h` with your application:

```
/* Substitute gsl_hypot for missing system hypot */

#ifndef HAVE_HYPOT
#define hypot gsl_hypot
#endif
```

The application source files can then use the include command `#include <config.h>` to replace each occurrence of `hypot()` by [`gsl_hypot()`](https://www.gnu.org/software/gsl/doc/html/math.html#c.gsl_hypot) when `hypot()` is not available. This substitution can be made automatically if you use `autoconf`, see [Autoconf Macros](https://www.gnu.org/software/gsl/doc/html/autoconf.html#chap-autoconf-macros).

In most circumstances the best strategy is to use the native versions of these functions when available, and fall back to GSL versions otherwise, since this allows your application to take advantage of any platform-specific optimizations in the system library. This is the strategy used within GSL itself.



## Alternative optimized functions

The main implementation of some functions in the library will not be optimal on all architectures. For example, there are several ways to compute a Gaussian random variate and their relative speeds are platform-dependent. In cases like this the library provides alternative implementations of these functions with the same interface. If you write your application using calls to the standard implementation you can select an alternative version later via a preprocessor definition. It is also possible to introduce your own optimized functions this way while retaining portability. The following lines demonstrate the use of a platform-dependent choice of methods for sampling from the Gaussian distribution:

```
#ifdef SPARC
#define gsl_ran_gaussian gsl_ran_gaussian_ratio_method
#endif
#ifdef INTEL
#define gsl_ran_gaussian my_gaussian
#endif
```

These lines would be placed in the configuration header file `config.h` of the application, which should then be included by all the source files. Note that the alternative implementations will not produce bit-for-bit identical results, and in the case of random number distributions will produce an entirely different stream of random variates.

## Support for different numeric types

Many functions in the library are defined for different numeric types. This feature is implemented by varying the name of the function with a type-related modifier—a primitive form of C++ templates. The modifier is inserted into the function name after the initial module prefix. The following table shows the function names defined for all the numeric types of an imaginary module `gsl_foo` with function `fn()`:

```
gsl_foo_fn               double
gsl_foo_long_double_fn   long double
gsl_foo_float_fn         float
gsl_foo_long_fn          long
gsl_foo_ulong_fn         unsigned long
gsl_foo_int_fn           int
gsl_foo_uint_fn          unsigned int
gsl_foo_short_fn         short
gsl_foo_ushort_fn        unsigned short
gsl_foo_char_fn          char
gsl_foo_uchar_fn         unsigned char
```

The normal numeric precision `double` is considered the default and does not require a suffix. For example, the function [`gsl_stats_mean()`](https://www.gnu.org/software/gsl/doc/html/statistics.html#c.gsl_stats_mean) computes the mean of double precision numbers, while the function `gsl_stats_int_mean()` computes the mean of integers.

A corresponding scheme is used for library defined types, such as `gsl_vector` and `gsl_matrix`. In this case the modifier is appended to the type name. For example, if a module defines a new type-dependent struct or typedef `gsl_foo` it is modified for other types in the following way:

```
gsl_foo                  double
gsl_foo_long_double      long double
gsl_foo_float            float
gsl_foo_long             long
gsl_foo_ulong            unsigned long
gsl_foo_int              int
gsl_foo_uint             unsigned int
gsl_foo_short            short
gsl_foo_ushort           unsigned short
gsl_foo_char             char
gsl_foo_uchar            unsigned char
```

When a module contains type-dependent definitions the library provides individual header files for each type. The filenames are modified as shown in the below. For convenience the default header includes the definitions for all the types. To include only the double precision header file, or any other specific type, use its individual filename:

```
#include <gsl/gsl_foo.h>               All types
#include <gsl/gsl_foo_double.h>        double
#include <gsl/gsl_foo_long_double.h>   long double
#include <gsl/gsl_foo_float.h>         float
#include <gsl/gsl_foo_long.h>          long
#include <gsl/gsl_foo_ulong.h>         unsigned long
#include <gsl/gsl_foo_int.h>           int
#include <gsl/gsl_foo_uint.h>          unsigned int
#include <gsl/gsl_foo_short.h>         short
#include <gsl/gsl_foo_ushort.h>        unsigned short
#include <gsl/gsl_foo_char.h>          char
#include <gsl/gsl_foo_uchar.h>         unsigned char
```



## Compatibility with C++

The library header files automatically define functions to have `extern "C"` linkage when included in C++ programs. This allows the functions to be called directly from C++.

To use C++ exception handling within user-defined functions passed to the library as parameters, the library must be built with the additional `CFLAGS` compilation option `-fexceptions`.





## 数组别名

The library assumes that arrays, vectors and matrices passed as modifiable arguments are not aliased and do not overlap with each other. This removes the need for the library to handle overlapping memory regions as a special case, and allows additional optimizations to be used. If overlapping memory regions are passed as modifiable arguments then the results of such functions will be undefined. If the arguments will not be modified (for example, if a function prototype declares them as `const` arguments) then overlapping or aliased memory regions can be safely used.

## 线程安全

The library can be used in multi-threaded programs. All the functions are thread-safe, in the sense that they do not use static variables. Memory is always associated with objects and not with functions. For functions which use *workspace* objects as temporary storage the workspaces should be allocated on a per-thread basis. For functions which use *table* objects as read-only memory the tables can be used by multiple threads simultaneously. Table arguments are always declared `const`in function prototypes, to indicate that they may be safely accessed by different threads.

本库可以在多线程程序中进行使用。所有函数都是线程安全的，只要它们不使用静态变量。内存永远与对象而不是函数相一致。对于使用*workspace*对象作为临时存储的函数工作空间应该在每个线程的基础上分配。列表变量在函数原型中总是被声明为`const`，以表明其在不同线程中是可以安全访问的。

There are a small number of static global variables which are used to control the overall behavior of the library (e.g. whether to use range-checking, the function to call on fatal error, etc). These variables are set directly by the user, so they should be initialized once at program startup and not modified by different threads.

有少量静态全局变量被用来控制库的整体行为(e.g. 是否使用范围检查，出现致命错误时调用的函数，等)。这些变阿玲由用户直接设置，因此它们应该在程序启东市被初始化病不能被不同线程更改。



## 弃用的函数

From time to time, it may be necessary for the definitions of some functions to be altered or removed from the library. In these circumstances the functions will first be declared *deprecated* and then removed from subsequent versions of the library. Functions that are deprecated can be disabled in the current release by setting the preprocessor definition `GSL_DISABLE_DEPRECATED`. This allows existing code to be tested for forwards compatibility.

有些情况下，也许需要从库中更改或移除一些函数的定义。在这些情况下函数首先会被声明为*deprecated(译注: 弃用)*并在随后的库版本中被移除。被弃用的函数可以在当前发布版本中通过设置预处理器定义`GSL_DISABLE_DEPRECATED`来禁用。这允许已存代码来进行前向兼容性测试。

## 代码复用

Where possible the routines in the library have been written to avoid dependencies between modules and files. This should make it possible to extract individual functions for use in your own applications, without needing to have the whole library installed. You may need to define certain macros such as `GSL_ERROR` and remove some `#include` statements in order to compile the files as standalone units. Reuse of the library code in this way is encouraged, subject to the terms of the GNU General Public License.

在任何有可能的地方本库中的程序都被编写来可以避免模块和文件之间的依赖。这使得不需要安装正割库而将独立的函数抽取出来用在自己的应用中是可行的。你也许需要定义特定的宏比如`GSL_ERROR`和移除一些`#include`声明来将文件能作为独立单元编译。以这个方式重用库代码是鼓励的，遵从GNU通用公共许可证的条款。

### 脚注

[[1]](https://www.gnu.org/software/gsl/doc/html/usage.html#id1)最后几位数字有可能会轻微不同取决于使用的编译器和平台—这是正常的

[[2]](https://www.gnu.org/software/gsl/doc/html/usage.html#id2)在MacOS X系统上不需要

[[3]](https://www.gnu.org/software/gsl/doc/html/usage.html#id3)<http://www.network-theory.co.uk/gcc/intro/>

[[4]](https://www.gnu.org/software/gsl/doc/html/usage.html#id4)在GNU/Linux系统上是`/etc/ld.so.conf`