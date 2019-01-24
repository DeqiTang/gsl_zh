## 简介

GNU科学库(GSL)是一个用于数值计算的程序集。这些程序使用C语言从零开始编写，并为C编程人员提供应用程序编程接口，也允许为更高级语言编写包装器。源码在GNU通用公共许可证下发布。

## GSL中可使用的程序

这个库涵盖数值计算中众多的主题。程序对于一下领域可以使用，

| Complex Numbers            | Roots of Polynomials      | Special Functions      |
| -------------------------- | ------------------------- | ---------------------- |
| Vectors and Matrices       | Permutations              | Combinations           |
| Sorting                    | BLAS Support              | Linear Algebra         |
| CBLAS Library              | Fast Fourier Transforms   | Eigensystems           |
| Random Numbers             | Quadrature                | Random Distributions   |
| Quasi-Random Sequences     | Histograms                | Statistics             |
| Monte Carlo Integration    | N-Tuples                  | Differential Equations |
| Simulated Annealing        | Numerical Differentiation | Interpolation          |
| Series Acceleration        | Chebyshev Approximations  | Root-Finding           |
| Discrete Hankel Transforms | Least-Squares Fitting     | Minimization           |
| IEEE Floating-Point        | Physical Constants        | Basis Splines          |
| Wavelets                   | Sparse BLAS Support       | Sparse Linear Algebra  |

这些程序的使用在本手册中得以描述。每章提供了函数的详细的定义，紧随其后的是示例程序以及所使用的算法基于的参考文章。

在一些情况下，程序是基于可靠的公用程序包比如FFTPACK和QUADPACK，由GSL的开发者使用C按照现代编程约定进行重新实现。

## GSL是自由软件

GNU科学计算库中的子程序属于“自由软件”；这意味着每个人可以免费使用以及在其它自由软件中重新分发该程序库。这个库不属于公用软件；它是受版权保护的并且是按照条件进行发布。这些条件被设计来用于允许一个好的合作者可以做任何其想做的事情。不被允许的是试图防止他人进一步分享他们可能从你这儿获取的本软件的任何版本。

Specifically, we want to make sure that you have the right to share copies of programs that you are given which use the GNU Scientific Library, that you receive their source code or else can get it if you want it, that you can change these programs or use pieces of them in new free programs, and that you know you can do these things.

特别地，我们想要确保你有权分享你得到的使用了GNU科学计算库的程序的拷贝，以及你能收到程序的源码或者在想要的情况下能够获取，以及你能够改变这些程序或者在新程序中使用它们的片段，以及你知道你能够做这些事情。

To make sure that everyone has such rights, we have to forbid you to deprive anyone else of these rights. For example, if you distribute copies of any code which uses the GNU Scientific Library, you must give the recipients all the rights that you have received. You must make sure that they, too, receive or can get the source code, both to the library and the code which uses it. And you must tell them their rights. This means that the library should not be redistributed in proprietary programs.



Also, for our own protection, we must make certain that everyone finds out that there is no warranty for the GNU Scientific Library. If these programs are modified by someone else and passed on, we want their recipients to know that what they have is not what we distributed, so that any problems introduced by others will not reflect on our reputation.

The precise conditions for the distribution of software related to the GNU Scientific Library are found in the [GNU General Public License](https://www.gnu.org/software/gsl/manual/html_node/GNU-General-Public-License.html#GNU-General-Public-License). Further information about this license is available from the GNU Project webpage [Frequently Asked Questions about the GNU GPL](http://www.gnu.org/copyleft/gpl-faq.html).

The Free Software Foundation also operates a license consulting service for commercial users (contact details available from [http://www.fsf.org](http://www.fsf.org/).

## Obtaining GSL

The source code for the library can be obtained in different ways, by copying it from a friend, purchasing it on CDROM or downloading it from the internet. A list of public ftp servers which carry the source code can be found on the GNU website, <http://www.gnu.org/software/gsl/>.

The preferred platform for the library is a GNU system, which allows it to take advantage of additional features in the GNU C compiler and GNU C library. However, the library is fully portable and should compile on most systems with a C compiler.

Announcements of new releases, updates and other relevant events are made on the [info-gsl@gnu.org](mailto:info-gsl%40gnu.org) mailing list. To subscribe to this low-volume list, send an email of the following form:

```
To: info-gsl-request@gnu.org
Subject: subscribe
```

You will receive a response asking you to reply in order to confirm your subscription.



## No Warranty

The software described in this manual has no warranty, it is provided “as is”. It is your responsibility to validate the behavior of the routines and their accuracy using the source code provided, or to purchase support and warranties from commercial redistributors. Consult the [GNU General Public License](https://www.gnu.org/software/gsl/manual/html_node/GNU-General-Public-License.html#GNU-General-Public-License) for further details.



## Reporting Bugs

A list of known bugs can be found in the `BUGS` file included in the GSL distribution or online in the GSL bug tracker. [[1\]](https://www.gnu.org/software/gsl/doc/html/intro.html#f1) Details of compilation problems can be found in the `INSTALL` file.

If you find a bug which is not listed in these files, please report it to [bug-gsl@gnu.org](mailto:bug-gsl%40gnu.org).

All bug reports should include:

- The version number of GSL
- The hardware and operating system
- The compiler used, including version number and compilation options
- A description of the bug behavior
- A short program which exercises the bug

It is useful if you can check whether the same problem occurs when the library is compiled without optimization. Thank you.

Any errors or omissions in this manual can also be reported to the same address.



## Further Information

Additional information, including online copies of this manual, links to related projects, and mailing list archives are available from the website mentioned above.

Any questions about the use and installation of the library can be asked on the mailing list [help-gsl@gnu.org](mailto:help-gsl%40gnu.org). To subscribe to this list, send an email of the following form:

```
To: help-gsl-request@gnu.org
Subject: subscribe
```

This mailing list can be used to ask questions not covered by this manual, and to contact the developers of the library.

If you would like to refer to the GNU Scientific Library in a journal article, the recommended way is to cite this reference manual, e.g.:

```
M. Galassi et al, GNU Scientific Library Reference Manual (3rd Ed.), ISBN 0954612078.
```

If you want to give a url, use “<http://www.gnu.org/software/gsl/>”.



## Conventions used in this manual

This manual contains many examples which can be typed at the keyboard. A command entered at the terminal is shown like this:

```
$ command
```

The first character on the line is the terminal prompt, and should not be typed. The dollar sign $ is used as the standard prompt in this manual, although some systems may use a different character.

The examples assume the use of the GNU operating system. There may be minor differences in the output on other systems. The commands for setting environment variables use the Bourne shell syntax of the standard GNU shell (`bash`).

### Footnotes

[1] <http://savannah.gnu.org/bugs/?group=gsl>