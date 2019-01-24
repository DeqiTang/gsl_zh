# Combinations

This chapter describes functions for creating and manipulating combinations. A combination ![c](https://www.gnu.org/software/gsl/doc/html/_images/math/a2ff0b909e008fecae54f397fc3bd6ef4ae96b5c.png) is represented by an array of ![k](https://www.gnu.org/software/gsl/doc/html/_images/math/90645b5c2abc13d7de104fd87abeb2f19406e5f8.png) integers in the range 0 to ![n - 1](https://www.gnu.org/software/gsl/doc/html/_images/math/f3e382739f6fe75af097ea5bc46a5d816872827d.png), where each value ![c_i](https://www.gnu.org/software/gsl/doc/html/_images/math/4cb4e830e3ae02f4a0ad77ce68a58f27915935e1.png) occurs at most once. The combination ![c](https://www.gnu.org/software/gsl/doc/html/_images/math/a2ff0b909e008fecae54f397fc3bd6ef4ae96b5c.png) corresponds to indices of ![k](https://www.gnu.org/software/gsl/doc/html/_images/math/90645b5c2abc13d7de104fd87abeb2f19406e5f8.png) elements chosen from an ![n](https://www.gnu.org/software/gsl/doc/html/_images/math/a24554f1502cf204e7aa24a6c064962c3504de48.png) element vector. Combinations are useful for iterating over all ![k](https://www.gnu.org/software/gsl/doc/html/_images/math/90645b5c2abc13d7de104fd87abeb2f19406e5f8.png)-element subsets of a set.

The functions described in this chapter are defined in the header file `gsl_combination.h`.

## The Combination struct

- `gsl_combination`

  A combination is defined by a structure containing three components, the values of ![n](https://www.gnu.org/software/gsl/doc/html/_images/math/a24554f1502cf204e7aa24a6c064962c3504de48.png) and ![k](https://www.gnu.org/software/gsl/doc/html/_images/math/90645b5c2abc13d7de104fd87abeb2f19406e5f8.png), and a pointer to the combination array. The elements of the combination array are all of type `size_t`, and are stored in increasing order. The [`gsl_combination`](https://www.gnu.org/software/gsl/doc/html/combination.html#c.gsl_combination) structure looks like this:`typedef struct {   size_t n;   size_t k;   size_t *data; } gsl_combination; `

## Combination allocation

- [gsl_combination](https://www.gnu.org/software/gsl/doc/html/combination.html#c.gsl_combination) * `gsl_combination_alloc`(size_t *n*, size_t *k*)

  This function allocates memory for a new combination with parameters `n`, `k`. The combination is not initialized and its elements are undefined. Use the function [`gsl_combination_calloc()`](https://www.gnu.org/software/gsl/doc/html/combination.html#c.gsl_combination_calloc) if you want to create a combination which is initialized to the lexicographically first combination. A null pointer is returned if insufficient memory is available to create the combination.

- [gsl_combination](https://www.gnu.org/software/gsl/doc/html/combination.html#c.gsl_combination) * `gsl_combination_calloc`(size_t *n*, size_t *k*)

  This function allocates memory for a new combination with parameters `n`, `k` and initializes it to the lexicographically first combination. A null pointer is returned if insufficient memory is available to create the combination.

- void `gsl_combination_init_first`([gsl_combination](https://www.gnu.org/software/gsl/doc/html/combination.html#c.gsl_combination) * *c*)

  This function initializes the combination `c` to the lexicographically first combination, i.e.![(0, 1, 2, \dots, k - 1)](https://www.gnu.org/software/gsl/doc/html/_images/math/416a638690b0a142c23639c08a34107507afa580.png).

- void `gsl_combination_init_last`([gsl_combination](https://www.gnu.org/software/gsl/doc/html/combination.html#c.gsl_combination) * *c*)

  This function initializes the combination `c` to the lexicographically last combination, i.e.![(n - k, n - k + 1, \dots, n - 1)](https://www.gnu.org/software/gsl/doc/html/_images/math/23309716d8e1d9bb3cbc7019c28e5c23834cc75b.png).

- void `gsl_combination_free`([gsl_combination](https://www.gnu.org/software/gsl/doc/html/combination.html#c.gsl_combination) * *c*)

  This function frees all the memory used by the combination `c`.

- int `gsl_combination_memcpy`([gsl_combination](https://www.gnu.org/software/gsl/doc/html/combination.html#c.gsl_combination) * *dest*, const [gsl_combination](https://www.gnu.org/software/gsl/doc/html/combination.html#c.gsl_combination) * *src*)

  This function copies the elements of the combination `src` into the combination `dest`. The two combinations must have the same size.

## Accessing combination elements

The following function can be used to access the elements of a combination.

- size_t `gsl_combination_get`(const [gsl_combination](https://www.gnu.org/software/gsl/doc/html/combination.html#c.gsl_combination) * *c*, const size_t *i*)

  This function returns the value of the `i`-th element of the combination `c`. If `i` lies outside the allowed range of 0 to ![k - 1](https://www.gnu.org/software/gsl/doc/html/_images/math/4730fcd0c2f6acdafd1ec39d67415edab5c38bbb.png) then the error handler is invoked and 0 is returned. An inline version of this function is used when `HAVE_INLINE` is defined.

## Combination properties

- size_t `gsl_combination_n`(const [gsl_combination](https://www.gnu.org/software/gsl/doc/html/combination.html#c.gsl_combination) * *c*)

  This function returns the range (![n](https://www.gnu.org/software/gsl/doc/html/_images/math/a24554f1502cf204e7aa24a6c064962c3504de48.png)) of the combination c.

- size_t `gsl_combination_k`(const [gsl_combination](https://www.gnu.org/software/gsl/doc/html/combination.html#c.gsl_combination) * *c*)

  This function returns the number of elements (![k](https://www.gnu.org/software/gsl/doc/html/_images/math/90645b5c2abc13d7de104fd87abeb2f19406e5f8.png)) in the combination `c`.

- size_t * `gsl_combination_data`(const [gsl_combination](https://www.gnu.org/software/gsl/doc/html/combination.html#c.gsl_combination) * *c*)

  This function returns a pointer to the array of elements in the combination `c`.



- int `gsl_combination_valid`([gsl_combination](https://www.gnu.org/software/gsl/doc/html/combination.html#c.gsl_combination) * *c*)

  This function checks that the combination `c` is valid. The `k` elements should lie in the range 0 to ![n - 1](https://www.gnu.org/software/gsl/doc/html/_images/math/f3e382739f6fe75af097ea5bc46a5d816872827d.png), with each value occurring once at most and in increasing order.

## Combination functions



- int `gsl_combination_next`([gsl_combination](https://www.gnu.org/software/gsl/doc/html/combination.html#c.gsl_combination) * *c*)

  This function advances the combination `c` to the next combination in lexicographic order and returns `GSL_SUCCESS`. If no further combinations are available it returns `GSL_FAILURE` and leaves`c` unmodified. Starting with the first combination and repeatedly applying this function will iterate through all possible combinations of a given order.

- int `gsl_combination_prev`([gsl_combination](https://www.gnu.org/software/gsl/doc/html/combination.html#c.gsl_combination) * *c*)

  This function steps backwards from the combination `c` to the previous combination in lexicographic order, returning `GSL_SUCCESS`. If no previous combination is available it returns`GSL_FAILURE` and leaves `c` unmodified.

## Reading and writing combinations

The library provides functions for reading and writing combinations to a file as binary data or formatted text.

- int `gsl_combination_fwrite`(FILE * *stream*, const [gsl_combination](https://www.gnu.org/software/gsl/doc/html/combination.html#c.gsl_combination) * *c*)

  This function writes the elements of the combination `c` to the stream `stream` in binary format. The function returns `GSL_EFAILED` if there was a problem writing to the file. Since the data is written in the native binary format it may not be portable between different architectures.

- int `gsl_combination_fread`(FILE * *stream*, [gsl_combination](https://www.gnu.org/software/gsl/doc/html/combination.html#c.gsl_combination) * *c*)

  This function reads elements from the open stream `stream` into the combination `c` in binary format. The combination `c` must be preallocated with correct values of ![n](https://www.gnu.org/software/gsl/doc/html/_images/math/a24554f1502cf204e7aa24a6c064962c3504de48.png) and ![k](https://www.gnu.org/software/gsl/doc/html/_images/math/90645b5c2abc13d7de104fd87abeb2f19406e5f8.png) since the function uses the size of `c` to determine how many bytes to read. The function returns `GSL_EFAILED` if there was a problem reading from the file. The data is assumed to have been written in the native binary format on the same architecture.

- int `gsl_combination_fprintf`(FILE * *stream*, const [gsl_combination](https://www.gnu.org/software/gsl/doc/html/combination.html#c.gsl_combination) * *c*, const char * *format*)

  This function writes the elements of the combination `c` line-by-line to the stream `stream`using the format specifier `format`, which should be suitable for a type of `size_t`. In ISO C99 the type modifier `z` represents `size_t`, so `"%zu\n"` is a suitable format [[1\]](https://www.gnu.org/software/gsl/doc/html/combination.html#f1). The function returns `GSL_EFAILED` if there was a problem writing to the file.

- int `gsl_combination_fscanf`(FILE * *stream*, [gsl_combination](https://www.gnu.org/software/gsl/doc/html/combination.html#c.gsl_combination) * *c*)

  This function reads formatted data from the stream `stream` into the combination `c`. The combination `c` must be preallocated with correct values of ![n](https://www.gnu.org/software/gsl/doc/html/_images/math/a24554f1502cf204e7aa24a6c064962c3504de48.png) and ![k](https://www.gnu.org/software/gsl/doc/html/_images/math/90645b5c2abc13d7de104fd87abeb2f19406e5f8.png) since the function uses the size of `c` to determine how many numbers to read. The function returns `GSL_EFAILED` if there was a problem reading from the file.

## Examples

The example program below prints all subsets of the set ![{0,1,2,3}](https://www.gnu.org/software/gsl/doc/html/_images/math/b59c77ffc893eda6c7ca7f0622c41e7815a64037.png) ordered by size. Subsets of the same size are ordered lexicographically.

```
#include <stdio.h>
#include <gsl/gsl_combination.h>

int
main (void)
{
  gsl_combination * c;
  size_t i;

  printf ("All subsets of {0,1,2,3} by size:\n") ;
  for (i = 0; i <= 4; i++)
    {
      c = gsl_combination_calloc (4, i);
      do
        {
          printf ("{");
          gsl_combination_fprintf (stdout, c, " %u");
          printf (" }\n");
        }
      while (gsl_combination_next (c) == GSL_SUCCESS);
      gsl_combination_free (c);
    }

  return 0;
}
```

Here is the output from the program,

```
All subsets of {0,1,2,3} by size:
{ }
{ 0 }
{ 1 }
{ 2 }
{ 3 }
{ 0 1 }
{ 0 2 }
{ 0 3 }
{ 1 2 }
{ 1 3 }
{ 2 3 }
{ 0 1 2 }
{ 0 1 3 }
{ 0 2 3 }
{ 1 2 3 }
{ 0 1 2 3 }
```

All 16 subsets are generated, and the subsets of each size are sorted lexicographically.

## References and Further Reading

Further information on combinations can be found in,

- Donald L. Kreher, Douglas R. Stinson, Combinatorial Algorithms: Generation, Enumeration and Search, 1998, CRC Press LLC, ISBN 084933988X

### Footnotes

[[1]](https://www.gnu.org/software/gsl/doc/html/combination.html#id1)In versions of the GNU C library prior to the ISO C99 standard, the type modifier `Z` was used instead.