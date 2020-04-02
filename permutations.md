# Permutations

This chapter describes functions for creating and manipulating permutations. A permutation ![p](https://www.gnu.org/software/gsl/doc/html/_images/math/af29d3e8554f748ed3f2d142a3ec9d4e038b9335.png) is represented by an array of ![n](https://www.gnu.org/software/gsl/doc/html/_images/math/a24554f1502cf204e7aa24a6c064962c3504de48.png) integers in the range 0 to ![n-1](https://www.gnu.org/software/gsl/doc/html/_images/math/0acbb742031c41a269215d223bd0d699a0cc0522.png), where each value ![p_i](https://www.gnu.org/software/gsl/doc/html/_images/math/1467ba8f9c83e8a7880d5878f586cee0d1753acf.png) occurs once and only once. The application of a permutation ![p](https://www.gnu.org/software/gsl/doc/html/_images/math/af29d3e8554f748ed3f2d142a3ec9d4e038b9335.png) to a vector ![v](https://www.gnu.org/software/gsl/doc/html/_images/math/a3ab041b283bfe2226babbccf0611362d5001323.png) yields a new vector ![v'](https://www.gnu.org/software/gsl/doc/html/_images/math/d407b116fdbda2ab5626258afc5c52ad270d74ea.png) where ![v'_i = v_{p_i}](https://www.gnu.org/software/gsl/doc/html/_images/math/285187811e4885ec945de5ecc9d0e680d4be0db9.png). For example, the array ![(0,1,3,2)](https://www.gnu.org/software/gsl/doc/html/_images/math/a617066f98aedbf95e7baee6fe3a6fbe7ddaa372.png) represents a permutation which exchanges the last two elements of a four element vector. The corresponding identity permutation is ![(0,1,2,3)](https://www.gnu.org/software/gsl/doc/html/_images/math/89e10001a8a6b262d1b1a238091e8253746a2d3d.png).

Note that the permutations produced by the linear algebra routines correspond to the exchange of matrix columns, and so should be considered as applying to row-vectors in the form ![v' = v P](https://www.gnu.org/software/gsl/doc/html/_images/math/7e55bca38cd932ddb5a672050104d4ff2f7a4917.png) rather than column-vectors, when permuting the elements of a vector.

The functions described in this chapter are defined in the header file `gsl_permutation.h`.

## The Permutation struct

- `gsl_permutation`

  A permutation is defined by a structure containing two components, the size of the permutation and a pointer to the permutation array. The elements of the permutation array are all of type `size_t`. The [`gsl_permutation`](https://www.gnu.org/software/gsl/doc/html/permutation.html#c.gsl_permutation) structure looks like this:`typedef struct {   size_t size;   size_t * data; } gsl_permutation; `

## Permutation allocation

- [gsl_permutation](https://www.gnu.org/software/gsl/doc/html/permutation.html#c.gsl_permutation) * `gsl_permutation_alloc`(size_t *n*)

  This function allocates memory for a new permutation of size `n`. The permutation is not initialized and its elements are undefined. Use the function [`gsl_permutation_calloc()`](https://www.gnu.org/software/gsl/doc/html/permutation.html#c.gsl_permutation_calloc) if you want to create a permutation which is initialized to the identity. A null pointer is returned if insufficient memory is available to create the permutation.

- [gsl_permutation](https://www.gnu.org/software/gsl/doc/html/permutation.html#c.gsl_permutation) * `gsl_permutation_calloc`(size_t *n*)

  This function allocates memory for a new permutation of size `n` and initializes it to the identity. A null pointer is returned if insufficient memory is available to create the permutation.



- void `gsl_permutation_init`([gsl_permutation](https://www.gnu.org/software/gsl/doc/html/permutation.html#c.gsl_permutation) * *p*)

  This function initializes the permutation `p` to the identity, i.e. ![(0, 1, 2, \dots, n - 1)](https://www.gnu.org/software/gsl/doc/html/_images/math/fb1a661db5da03ccf9323547595d314a2154d258.png).

- void `gsl_permutation_free`([gsl_permutation](https://www.gnu.org/software/gsl/doc/html/permutation.html#c.gsl_permutation) * *p*)

  This function frees all the memory used by the permutation `p`.

- int `gsl_permutation_memcpy`([gsl_permutation](https://www.gnu.org/software/gsl/doc/html/permutation.html#c.gsl_permutation) * *dest*, const [gsl_permutation](https://www.gnu.org/software/gsl/doc/html/permutation.html#c.gsl_permutation) * *src*)

  This function copies the elements of the permutation `src` into the permutation `dest`. The two permutations must have the same size.

## Accessing permutation elements

The following functions can be used to access and manipulate permutations.

- size_t `gsl_permutation_get`(const [gsl_permutation](https://www.gnu.org/software/gsl/doc/html/permutation.html#c.gsl_permutation) * *p*, const size_t *i*)

  This function returns the value of the `i`-th element of the permutation `p`. If `i` lies outside the allowed range of 0 to ![n - 1](https://www.gnu.org/software/gsl/doc/html/_images/math/f3e382739f6fe75af097ea5bc46a5d816872827d.png) then the error handler is invoked and 0 is returned. An inline version of this function is used when `HAVE_INLINE` is defined.



- int `gsl_permutation_swap`([gsl_permutation](https://www.gnu.org/software/gsl/doc/html/permutation.html#c.gsl_permutation) * *p*, const size_t *i*, const size_t *j*)

  This function exchanges the `i`-th and `j`-th elements of the permutation `p`.

## Permutation properties

- size_t `gsl_permutation_size`(const [gsl_permutation](https://www.gnu.org/software/gsl/doc/html/permutation.html#c.gsl_permutation) * *p*)

  This function returns the size of the permutation `p`.

- size_t * `gsl_permutation_data`(const [gsl_permutation](https://www.gnu.org/software/gsl/doc/html/permutation.html#c.gsl_permutation) * *p*)

  This function returns a pointer to the array of elements in the permutation `p`.



- int `gsl_permutation_valid`(const [gsl_permutation](https://www.gnu.org/software/gsl/doc/html/permutation.html#c.gsl_permutation) * *p*)

  This function checks that the permutation `p` is valid. The `n` elements should contain each of the numbers 0 to `n - 1` once and only once.

## Permutation functions



- void `gsl_permutation_reverse`([gsl_permutation](https://www.gnu.org/software/gsl/doc/html/permutation.html#c.gsl_permutation) * *p*)

  This function reverses the elements of the permutation `p`.



- int `gsl_permutation_inverse`([gsl_permutation](https://www.gnu.org/software/gsl/doc/html/permutation.html#c.gsl_permutation) * *inv*, const [gsl_permutation](https://www.gnu.org/software/gsl/doc/html/permutation.html#c.gsl_permutation) * *p*)

  This function computes the inverse of the permutation `p`, storing the result in `inv`.



- int `gsl_permutation_next`([gsl_permutation](https://www.gnu.org/software/gsl/doc/html/permutation.html#c.gsl_permutation) * *p*)

  This function advances the permutation `p` to the next permutation in lexicographic order and returns `GSL_SUCCESS`. If no further permutations are available it returns `GSL_FAILURE` and leaves`p` unmodified. Starting with the identity permutation and repeatedly applying this function will iterate through all possible permutations of a given order.

- int `gsl_permutation_prev`([gsl_permutation](https://www.gnu.org/software/gsl/doc/html/permutation.html#c.gsl_permutation) * *p*)

  This function steps backwards from the permutation `p` to the previous permutation in lexicographic order, returning `GSL_SUCCESS`. If no previous permutation is available it returns`GSL_FAILURE` and leaves `p` unmodified.

## Applying Permutations

- int `gsl_permute`(const size_t * *p*, double * *data*, size_t *stride*, size_t *n*)

  This function applies the permutation `p` to the array `data` of size `n` with stride `stride`.

- int `gsl_permute_inverse`(const size_t * *p*, double * *data*, size_t *stride*, size_t *n*)

  This function applies the inverse of the permutation `p` to the array `data` of size `n` with stride `stride`.

- int `gsl_permute_vector`(const [gsl_permutation](https://www.gnu.org/software/gsl/doc/html/permutation.html#c.gsl_permutation) * *p*, [gsl_vector](https://www.gnu.org/software/gsl/doc/html/vectors.html#c.gsl_vector) * *v*)

  This function applies the permutation `p` to the elements of the vector `v`, considered as a row-vector acted on by a permutation matrix from the right, ![v' = v P](https://www.gnu.org/software/gsl/doc/html/_images/math/7e55bca38cd932ddb5a672050104d4ff2f7a4917.png). The ![j](https://www.gnu.org/software/gsl/doc/html/_images/math/7bc2ad053edeecad39b867d60495ac2caa4555dc.png)-th column of the permutation matrix ![P](https://www.gnu.org/software/gsl/doc/html/_images/math/d1a0a8c0064fe114c19c8fd71ffec7a98480f707.png) is given by the ![p_j](https://www.gnu.org/software/gsl/doc/html/_images/math/81d175bb9fa9d51aa8125101ed2c69af5f0c6bae.png)-th column of the identity matrix. The permutation `p`and the vector `v` must have the same length.

- int `gsl_permute_vector_inverse`(const [gsl_permutation](https://www.gnu.org/software/gsl/doc/html/permutation.html#c.gsl_permutation) * *p*, [gsl_vector](https://www.gnu.org/software/gsl/doc/html/vectors.html#c.gsl_vector) * *v*)

  This function applies the inverse of the permutation `p` to the elements of the vector `v`, considered as a row-vector acted on by an inverse permutation matrix from the right, ![v' = v P^T](https://www.gnu.org/software/gsl/doc/html/_images/math/636fd314735f8be20409a7690968c53723172ee8.png). Note that for permutation matrices the inverse is the same as the transpose. The ![j](https://www.gnu.org/software/gsl/doc/html/_images/math/7bc2ad053edeecad39b867d60495ac2caa4555dc.png)-th column of the permutation matrix ![P](https://www.gnu.org/software/gsl/doc/html/_images/math/d1a0a8c0064fe114c19c8fd71ffec7a98480f707.png) is given by the ![:data:`p`_j](https://www.gnu.org/software/gsl/doc/html/_images/math/ebcdc5ce26c84f26052f066c67a88340c9c50f6b.png)-th column of the identity matrix. The permutation `p` and the vector `v` must have the same length.

- int `gsl_permute_matrix`(const [gsl_permutation](https://www.gnu.org/software/gsl/doc/html/permutation.html#c.gsl_permutation) * *p*, [gsl_matrix](https://www.gnu.org/software/gsl/doc/html/vectors.html#c.gsl_matrix) * *A*)

  This function applies the permutation `p` to the matrix `A` from the right, ![A' = A P](https://www.gnu.org/software/gsl/doc/html/_images/math/df1af50258a65ee2773538e13153509c6fc394b5.png). The ![j](https://www.gnu.org/software/gsl/doc/html/_images/math/7bc2ad053edeecad39b867d60495ac2caa4555dc.png)-th column of the permutation matrix ![P](https://www.gnu.org/software/gsl/doc/html/_images/math/d1a0a8c0064fe114c19c8fd71ffec7a98480f707.png) is given by the ![p_j](https://www.gnu.org/software/gsl/doc/html/_images/math/81d175bb9fa9d51aa8125101ed2c69af5f0c6bae.png)-th column of the identity matrix. This effectively permutes the columns of `A` according to the permutation `p`, and so the number of columns of `A` must equal the size of the permutation `p`.

- int `gsl_permutation_mul`([gsl_permutation](https://www.gnu.org/software/gsl/doc/html/permutation.html#c.gsl_permutation) * *p*, const [gsl_permutation](https://www.gnu.org/software/gsl/doc/html/permutation.html#c.gsl_permutation) * *pa*, const [gsl_permutation](https://www.gnu.org/software/gsl/doc/html/permutation.html#c.gsl_permutation) * *pb*)

  This function combines the two permutations `pa` and `pb` into a single permutation `p`, where ![p = pa * pb](https://www.gnu.org/software/gsl/doc/html/_images/math/4ff77ea13cd913ffddc8cd297770b3bbd3dc9209.png) The permutation `p` is equivalent to applying `pb` first and then `pa`.

## Reading and writing permutations

The library provides functions for reading and writing permutations to a file as binary data or formatted text.

- int `gsl_permutation_fwrite`(FILE * *stream*, const [gsl_permutation](https://www.gnu.org/software/gsl/doc/html/permutation.html#c.gsl_permutation) * *p*)

  This function writes the elements of the permutation `p` to the stream `stream` in binary format. The function returns `GSL_EFAILED` if there was a problem writing to the file. Since the data is written in the native binary format it may not be portable between different architectures.

- int `gsl_permutation_fread`(FILE * *stream*, [gsl_permutation](https://www.gnu.org/software/gsl/doc/html/permutation.html#c.gsl_permutation) * *p*)

  This function reads into the permutation `p` from the open stream `stream` in binary format. The permutation `p` must be preallocated with the correct length since the function uses the size of `p` to determine how many bytes to read. The function returns `GSL_EFAILED` if there was a problem reading from the file. The data is assumed to have been written in the native binary format on the same architecture.

- int `gsl_permutation_fprintf`(FILE * *stream*, const [gsl_permutation](https://www.gnu.org/software/gsl/doc/html/permutation.html#c.gsl_permutation) * *p*, const char * *format*)

  This function writes the elements of the permutation `p` line-by-line to the stream `stream`using the format specifier `format`, which should be suitable for a type of `size_t`. In ISO C99 the type modifier `z` represents `size_t`, so `"%zu\n"` is a suitable format [[1\]](https://www.gnu.org/software/gsl/doc/html/permutation.html#f1). The function returns `GSL_EFAILED` if there was a problem writing to the file.

- int `gsl_permutation_fscanf`(FILE * *stream*, [gsl_permutation](https://www.gnu.org/software/gsl/doc/html/permutation.html#c.gsl_permutation) * *p*)

  This function reads formatted data from the stream `stream` into the permutation `p`. The permutation `p` must be preallocated with the correct length since the function uses the size of `p` to determine how many numbers to read. The function returns `GSL_EFAILED` if there was a problem reading from the file.

## Permutations in cyclic form

A permutation can be represented in both *linear* and *cyclic* notations. The functions described in this section convert between the two forms. The linear notation is an index mapping, and has already been described above. The cyclic notation expresses a permutation as a series of circular rearrangements of groups of elements, or *cycles*.

For example, under the cycle (1 2 3), 1 is replaced by 2, 2 is replaced by 3 and 3 is replaced by 1 in a circular fashion. Cycles of different sets of elements can be combined independently, for example (1 2 3) (4 5) combines the cycle (1 2 3) with the cycle (4 5), which is an exchange of elements 4 and 5. A cycle of length one represents an element which is unchanged by the permutation and is referred to as a *singleton*.

It can be shown that every permutation can be decomposed into combinations of cycles. The decomposition is not unique, but can always be rearranged into a standard *canonical form* by a reordering of elements. The library uses the canonical form defined in Knuth’s *Art of Computer Programming* (Vol 1, 3rd Ed, 1997) Section 1.3.3, p.178.

The procedure for obtaining the canonical form given by Knuth is,

1. Write all singleton cycles explicitly
2. Within each cycle, put the smallest number first
3. Order the cycles in decreasing order of the first number in the cycle.

For example, the linear representation (2 4 3 0 1) is represented as (1 4) (0 2 3) in canonical form. The permutation corresponds to an exchange of elements 1 and 4, and rotation of elements 0, 2 and 3.

The important property of the canonical form is that it can be reconstructed from the contents of each cycle without the brackets. In addition, by removing the brackets it can be considered as a linear representation of a different permutation. In the example given above the permutation (2 4 3 0 1) would become (1 4 0 2 3). This mapping has many applications in the theory of permutations.

- int `gsl_permutation_linear_to_canonical`([gsl_permutation](https://www.gnu.org/software/gsl/doc/html/permutation.html#c.gsl_permutation) * *q*, const [gsl_permutation](https://www.gnu.org/software/gsl/doc/html/permutation.html#c.gsl_permutation) * *p*)

  This function computes the canonical form of the permutation `p` and stores it in the output argument `q`.

- int `gsl_permutation_canonical_to_linear`([gsl_permutation](https://www.gnu.org/software/gsl/doc/html/permutation.html#c.gsl_permutation) * *p*, const [gsl_permutation](https://www.gnu.org/software/gsl/doc/html/permutation.html#c.gsl_permutation) * *q*)

  This function converts a permutation `q` in canonical form back into linear form storing it in the output argument `p`.

- size_t `gsl_permutation_inversions`(const [gsl_permutation](https://www.gnu.org/software/gsl/doc/html/permutation.html#c.gsl_permutation) * *p*)

  This function counts the number of inversions in the permutation `p`. An inversion is any pair of elements that are not in order. For example, the permutation 2031 has three inversions, corresponding to the pairs (2,0) (2,1) and (3,1). The identity permutation has no inversions.

- size_t `gsl_permutation_linear_cycles`(const [gsl_permutation](https://www.gnu.org/software/gsl/doc/html/permutation.html#c.gsl_permutation) * *p*)

  This function counts the number of cycles in the permutation `p`, given in linear form.

- size_t `gsl_permutation_canonical_cycles`(const [gsl_permutation](https://www.gnu.org/software/gsl/doc/html/permutation.html#c.gsl_permutation) * *q*)

  This function counts the number of cycles in the permutation `q`, given in canonical form.

## Examples

The example program below creates a random permutation (by shuffling the elements of the identity) and finds its inverse.

```
#include <stdio.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>
#include <gsl/gsl_permutation.h>

int
main (void)
{
  const size_t N = 10;
  const gsl_rng_type * T;
  gsl_rng * r;

  gsl_permutation * p = gsl_permutation_alloc (N);
  gsl_permutation * q = gsl_permutation_alloc (N);

  gsl_rng_env_setup();
  T = gsl_rng_default;
  r = gsl_rng_alloc (T);

  printf ("initial permutation:");
  gsl_permutation_init (p);
  gsl_permutation_fprintf (stdout, p, " %u");
  printf ("\n");

  printf (" random permutation:");
  gsl_ran_shuffle (r, p->data, N, sizeof(size_t));
  gsl_permutation_fprintf (stdout, p, " %u");
  printf ("\n");

  printf ("inverse permutation:");
  gsl_permutation_inverse (q, p);
  gsl_permutation_fprintf (stdout, q, " %u");
  printf ("\n");

  gsl_permutation_free (p);
  gsl_permutation_free (q);
  gsl_rng_free (r);

  return 0;
}
```

Here is the output from the program:

```
$ ./a.out
initial permutation: 0 1 2 3 4 5 6 7 8 9
 random permutation: 1 3 5 2 7 6 0 4 9 8
inverse permutation: 6 0 3 1 7 2 5 4 9 8
```

The random permutation `p[i]` and its inverse `q[i]` are related through the identity `p[q[i]] = i`, which can be verified from the output.

The next example program steps forwards through all possible third order permutations, starting from the identity,

```
#include <stdio.h>
#include <gsl/gsl_permutation.h>

int
main (void)
{
  gsl_permutation * p = gsl_permutation_alloc (3);

  gsl_permutation_init (p);

  do
   {
      gsl_permutation_fprintf (stdout, p, " %u");
      printf ("\n");
   }
  while (gsl_permutation_next(p) == GSL_SUCCESS);

  gsl_permutation_free (p);

  return 0;
}
```

Here is the output from the program:

```
$ ./a.out
 0 1 2
 0 2 1
 1 0 2
 1 2 0
 2 0 1
 2 1 0
```

The permutations are generated in lexicographic order. To reverse the sequence, begin with the final permutation (which is the reverse of the identity) and replace [`gsl_permutation_next()`](https://www.gnu.org/software/gsl/doc/html/permutation.html#c.gsl_permutation_next) with[`gsl_permutation_prev()`](https://www.gnu.org/software/gsl/doc/html/permutation.html#c.gsl_permutation_prev).

## References and Further Reading

The subject of permutations is covered extensively in the following,

- Donald E. Knuth, The Art of Computer Programming: Sorting and Searching (Vol 3, 3rd Ed, 1997), Addison-Wesley, ISBN 0201896850.

For the definition of the *canonical form* see,

- Donald E. Knuth, The Art of Computer Programming: Fundamental Algorithms (Vol 1, 3rd Ed, 1997), Addison-Wesley, ISBN 0201896850. Section 1.3.3, An Unusual Correspondence, p.178–179.

### Footnotes

[[1]](https://www.gnu.org/software/gsl/doc/html/permutation.html#id1)In versions of the GNU C library prior to the ISO C99 standard, the type modifier `Z` was used instead.