# Multisets

This chapter describes functions for creating and manipulating multisets. A multiset ![c](https://www.gnu.org/software/gsl/doc/html/_images/math/a2ff0b909e008fecae54f397fc3bd6ef4ae96b5c.png) is represented by an array of ![k](https://www.gnu.org/software/gsl/doc/html/_images/math/90645b5c2abc13d7de104fd87abeb2f19406e5f8.png) integers in the range 0 to ![n - 1](https://www.gnu.org/software/gsl/doc/html/_images/math/f3e382739f6fe75af097ea5bc46a5d816872827d.png), where each value ![c_i](https://www.gnu.org/software/gsl/doc/html/_images/math/4cb4e830e3ae02f4a0ad77ce68a58f27915935e1.png) may occur more than once. The multiset ![c](https://www.gnu.org/software/gsl/doc/html/_images/math/a2ff0b909e008fecae54f397fc3bd6ef4ae96b5c.png) corresponds to indices of ![k](https://www.gnu.org/software/gsl/doc/html/_images/math/90645b5c2abc13d7de104fd87abeb2f19406e5f8.png) elements chosen from an ![n](https://www.gnu.org/software/gsl/doc/html/_images/math/a24554f1502cf204e7aa24a6c064962c3504de48.png) element vector with replacement. In mathematical terms, ![n](https://www.gnu.org/software/gsl/doc/html/_images/math/a24554f1502cf204e7aa24a6c064962c3504de48.png) is the cardinality of the multiset while ![k](https://www.gnu.org/software/gsl/doc/html/_images/math/90645b5c2abc13d7de104fd87abeb2f19406e5f8.png) is the maximum multiplicity of any value. Multisets are useful, for example, when iterating over the indices of a ![k](https://www.gnu.org/software/gsl/doc/html/_images/math/90645b5c2abc13d7de104fd87abeb2f19406e5f8.png)-th order symmetric tensor in ![n](https://www.gnu.org/software/gsl/doc/html/_images/math/a24554f1502cf204e7aa24a6c064962c3504de48.png)-space.

The functions described in this chapter are defined in the header file `gsl_multiset.h`.

## The Multiset struct

- `gsl_multiset`

  A multiset is defined by a structure containing three components, the values of ![n](https://www.gnu.org/software/gsl/doc/html/_images/math/a24554f1502cf204e7aa24a6c064962c3504de48.png) and ![k](https://www.gnu.org/software/gsl/doc/html/_images/math/90645b5c2abc13d7de104fd87abeb2f19406e5f8.png), and a pointer to the multiset array. The elements of the multiset array are all of type `size_t`, and are stored in increasing order. The [`gsl_multiset`](https://www.gnu.org/software/gsl/doc/html/multiset.html#c.gsl_multiset) structure looks like this:`typedef struct {   size_t n;   size_t k;   size_t *data; } gsl_multiset; `

## Multiset allocation

- [gsl_multiset](https://www.gnu.org/software/gsl/doc/html/multiset.html#c.gsl_multiset) * `gsl_multiset_alloc`(size_t *n*, size_t *k*)

  This function allocates memory for a new multiset with parameters `n`, `k`. The multiset is not initialized and its elements are undefined. Use the function [`gsl_multiset_calloc()`](https://www.gnu.org/software/gsl/doc/html/multiset.html#c.gsl_multiset_calloc) if you want to create a multiset which is initialized to the lexicographically first multiset element. A null pointer is returned if insufficient memory is available to create the multiset.

- [gsl_multiset](https://www.gnu.org/software/gsl/doc/html/multiset.html#c.gsl_multiset) * `gsl_multiset_calloc`(size_t *n*, size_t *k*)

  This function allocates memory for a new multiset with parameters `n`, `k` and initializes it to the lexicographically first multiset element. A null pointer is returned if insufficient memory is available to create the multiset.

- void `gsl_multiset_init_first`([gsl_multiset](https://www.gnu.org/software/gsl/doc/html/multiset.html#c.gsl_multiset) * *c*)

  This function initializes the multiset `c` to the lexicographically first multiset element, i.e. ![0](https://www.gnu.org/software/gsl/doc/html/_images/math/3b9bed1b0ccefe4f2e5e9be7ff710ff847c91a5c.png)repeated ![k](https://www.gnu.org/software/gsl/doc/html/_images/math/90645b5c2abc13d7de104fd87abeb2f19406e5f8.png) times.

- void `gsl_multiset_init_last`([gsl_multiset](https://www.gnu.org/software/gsl/doc/html/multiset.html#c.gsl_multiset) * *c*)

  This function initializes the multiset `c` to the lexicographically last multiset element, i.e. ![n-1](https://www.gnu.org/software/gsl/doc/html/_images/math/0acbb742031c41a269215d223bd0d699a0cc0522.png)repeated ![k](https://www.gnu.org/software/gsl/doc/html/_images/math/90645b5c2abc13d7de104fd87abeb2f19406e5f8.png) times.

- void `gsl_multiset_free`([gsl_multiset](https://www.gnu.org/software/gsl/doc/html/multiset.html#c.gsl_multiset) * *c*)

  This function frees all the memory used by the multiset `c`.

- int `gsl_multiset_memcpy`([gsl_multiset](https://www.gnu.org/software/gsl/doc/html/multiset.html#c.gsl_multiset) * *dest*, const [gsl_multiset](https://www.gnu.org/software/gsl/doc/html/multiset.html#c.gsl_multiset) * *src*)

  This function copies the elements of the multiset `src` into the multiset `dest`. The two multisets must have the same size.

## Accessing multiset elements

The following function can be used to access the elements of a multiset.

- size_t `gsl_multiset_get`(const [gsl_multiset](https://www.gnu.org/software/gsl/doc/html/multiset.html#c.gsl_multiset) * *c*, const size_t *i*)

  This function returns the value of the `i`-th element of the multiset `c`. If `i` lies outside the allowed range of 0 to ![k - 1](https://www.gnu.org/software/gsl/doc/html/_images/math/4730fcd0c2f6acdafd1ec39d67415edab5c38bbb.png) then the error handler is invoked and 0 is returned. An inline version of this function is used when `HAVE_INLINE` is defined.

## Multiset properties

- size_t `gsl_multiset_n`(const [gsl_multiset](https://www.gnu.org/software/gsl/doc/html/multiset.html#c.gsl_multiset) * *c*)

  This function returns the range (![n](https://www.gnu.org/software/gsl/doc/html/_images/math/a24554f1502cf204e7aa24a6c064962c3504de48.png)) of the multiset `c`.

- size_t `gsl_multiset_k`(const [gsl_multiset](https://www.gnu.org/software/gsl/doc/html/multiset.html#c.gsl_multiset) * *c*)

  This function returns the number of elements (![k](https://www.gnu.org/software/gsl/doc/html/_images/math/90645b5c2abc13d7de104fd87abeb2f19406e5f8.png)) in the multiset `c`.

- size_t * `gsl_multiset_data`(const [gsl_multiset](https://www.gnu.org/software/gsl/doc/html/multiset.html#c.gsl_multiset) * *c*)

  This function returns a pointer to the array of elements in the multiset `c`.



- int `gsl_multiset_valid`([gsl_multiset](https://www.gnu.org/software/gsl/doc/html/multiset.html#c.gsl_multiset) * *c*)

  This function checks that the multiset `c` is valid. The `k` elements should lie in the range 0 to ![n - 1](https://www.gnu.org/software/gsl/doc/html/_images/math/f3e382739f6fe75af097ea5bc46a5d816872827d.png), with each value occurring in nondecreasing order.

## Multiset functions



- int `gsl_multiset_next`([gsl_multiset](https://www.gnu.org/software/gsl/doc/html/multiset.html#c.gsl_multiset) * *c*)

  This function advances the multiset `c` to the next multiset element in lexicographic order and returns `GSL_SUCCESS`. If no further multisets elements are available it returns `GSL_FAILURE` and leaves `c` unmodified. Starting with the first multiset and repeatedly applying this function will iterate through all possible multisets of a given order.

- int `gsl_multiset_prev`([gsl_multiset](https://www.gnu.org/software/gsl/doc/html/multiset.html#c.gsl_multiset) * *c*)

  This function steps backwards from the multiset `c` to the previous multiset element in lexicographic order, returning `GSL_SUCCESS`. If no previous multiset is available it returns `GSL_FAILURE` and leaves `c` unmodified.

## Reading and writing multisets

The library provides functions for reading and writing multisets to a file as binary data or formatted text.

- int `gsl_multiset_fwrite`(FILE * *stream*, const [gsl_multiset](https://www.gnu.org/software/gsl/doc/html/multiset.html#c.gsl_multiset) * *c*)

  This function writes the elements of the multiset `c` to the stream `stream` in binary format. The function returns `GSL_EFAILED` if there was a problem writing to the file. Since the data is written in the native binary format it may not be portable between different architectures.

- int `gsl_multiset_fread`(FILE * *stream*, [gsl_multiset](https://www.gnu.org/software/gsl/doc/html/multiset.html#c.gsl_multiset) * *c*)

  This function reads elements from the open stream `stream` into the multiset `c` in binary format. The multiset `c` must be preallocated with correct values of ![n](https://www.gnu.org/software/gsl/doc/html/_images/math/a24554f1502cf204e7aa24a6c064962c3504de48.png) and ![k](https://www.gnu.org/software/gsl/doc/html/_images/math/90645b5c2abc13d7de104fd87abeb2f19406e5f8.png) since the function uses the size of `c` to determine how many bytes to read. The function returns `GSL_EFAILED` if there was a problem reading from the file. The data is assumed to have been written in the native binary format on the same architecture.

- int `gsl_multiset_fprintf`(FILE * *stream*, const [gsl_multiset](https://www.gnu.org/software/gsl/doc/html/multiset.html#c.gsl_multiset) * *c*, const char * *format*)

  This function writes the elements of the multiset `c` line-by-line to the stream `stream` using the format specifier `format`, which should be suitable for a type of `size_t`. In ISO C99 the type modifier `z` represents `size_t`, so `"%zu\n"` is a suitable format [[1\]](https://www.gnu.org/software/gsl/doc/html/multiset.html#f1). The function returns `GSL_EFAILED` if there was a problem writing to the file.

- int `gsl_multiset_fscanf`(FILE * *stream*, [gsl_multiset](https://www.gnu.org/software/gsl/doc/html/multiset.html#c.gsl_multiset) * *c*)

  This function reads formatted data from the stream `stream` into the multiset `c`. The multiset `c` must be preallocated with correct values of ![n](https://www.gnu.org/software/gsl/doc/html/_images/math/a24554f1502cf204e7aa24a6c064962c3504de48.png) and ![k](https://www.gnu.org/software/gsl/doc/html/_images/math/90645b5c2abc13d7de104fd87abeb2f19406e5f8.png) since the function uses the size of `c` to determine how many numbers to read. The function returns `GSL_EFAILED` if there was a problem reading from the file.

## Examples

The example program below prints all multisets elements containing the values ![{0,1,2,3}](https://www.gnu.org/software/gsl/doc/html/_images/math/b59c77ffc893eda6c7ca7f0622c41e7815a64037.png) ordered by size. Multiset elements of the same size are ordered lexicographically.

```
#include <stdio.h>
#include <gsl/gsl_multiset.h>

int
main (void)
{
  gsl_multiset * c;
  size_t i;

  printf ("All multisets of {0,1,2,3} by size:\n") ;
  for (i = 0; i <= 4; i++)
    {
      c = gsl_multiset_calloc (4, i);
      do
        {
          printf ("{");
          gsl_multiset_fprintf (stdout, c, " %u");
          printf (" }\n");
        }
      while (gsl_multiset_next (c) == GSL_SUCCESS);
      gsl_multiset_free (c);
    }

  return 0;
}
```

Here is the output from the program,

```
All multisets of {0,1,2,3} by size:
{ }
{ 0 }
{ 1 }
{ 2 }
{ 3 }
{ 0 0 }
{ 0 1 }
{ 0 2 }
{ 0 3 }
{ 1 1 }
{ 1 2 }
{ 1 3 }
{ 2 2 }
{ 2 3 }
{ 3 3 }
{ 0 0 0 }
{ 0 0 1 }
{ 0 0 2 }
{ 0 0 3 }
{ 0 1 1 }
{ 0 1 2 }
{ 0 1 3 }
{ 0 2 2 }
{ 0 2 3 }
{ 0 3 3 }
{ 1 1 1 }
{ 1 1 2 }
{ 1 1 3 }
{ 1 2 2 }
{ 1 2 3 }
{ 1 3 3 }
{ 2 2 2 }
{ 2 2 3 }
{ 2 3 3 }
{ 3 3 3 }
{ 0 0 0 0 }
{ 0 0 0 1 }
{ 0 0 0 2 }
{ 0 0 0 3 }
{ 0 0 1 1 }
{ 0 0 1 2 }
{ 0 0 1 3 }
{ 0 0 2 2 }
{ 0 0 2 3 }
{ 0 0 3 3 }
{ 0 1 1 1 }
{ 0 1 1 2 }
{ 0 1 1 3 }
{ 0 1 2 2 }
{ 0 1 2 3 }
{ 0 1 3 3 }
{ 0 2 2 2 }
{ 0 2 2 3 }
{ 0 2 3 3 }
{ 0 3 3 3 }
{ 1 1 1 1 }
{ 1 1 1 2 }
{ 1 1 1 3 }
{ 1 1 2 2 }
{ 1 1 2 3 }
{ 1 1 3 3 }
{ 1 2 2 2 }
{ 1 2 2 3 }
{ 1 2 3 3 }
{ 1 3 3 3 }
{ 2 2 2 2 }
{ 2 2 2 3 }
{ 2 2 3 3 }
{ 2 3 3 3 }
{ 3 3 3 3 }
```

All 70 multisets are generated and sorted lexicographically.

### Footnotes

[[1]](https://www.gnu.org/software/gsl/doc/html/multiset.html#id1)In versions of the GNU C library prior to the ISO C99 standard, the type modifier `Z` was used instead.