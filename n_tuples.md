# N-tuplesN重元组

This chapter describes functions for creating and manipulating *ntuples*, sets of values associated with events. The ntuples are stored in files. Their values can be extracted in any combination and *booked* in a histogram using a selection function.

本章描述用于创造和操控*ntuples*的函数，与事件相关的值的集合。N重元组是存储在文件中。它们的值可以以任何组合提取以及使用一个选择函数在直方图中进行登记。

The values to be stored are held in a user-defined data structure, and an ntuple is created associating this data structure with a file. The values are then written to the file (normally inside a loop) using the ntuple functions described below.

要存储的值被保存在一个用户定义数据结构中，并且一个N重元组是被创造以将该数据结构与一个文件相关联。这些值进而会被通过下面描述的N重元组函数写入文件(通常是在一个循环中)。

A histogram can be created from ntuple data by providing a selection function and a value function. The selection function specifies whether an event should be included in the subset to be analyzed or not. The value function computes the entry to be added to the histogram for each event.

通过提供一个选择函数和一个值函数可以从N重元组数据中产生一个直方图。选择函数确定了是否一个时间应该被包含在将要被分析的子集中。值函数计算每个实践中将要加入直方图的条目的值。

All the ntuple functions are defined in the header file `gsl_ntuple.h`.

所有的N重元组函数被声明在头文件`gsl_ntuple.h`中。

## N重元组结构

- `gsl_ntuple`

  Ntuples are manipulated using the [`gsl_ntuple`](https://www.gnu.org/software/gsl/doc/html/ntuple.html#c.gsl_ntuple) struct. This struct contains information on the file where the ntuple data is stored, a pointer to the current ntuple data row and the size of the user-defined ntuple data struct:`typedef struct   {     FILE * file;     void * ntuple_data;     size_t size;   } gsl_ntuple; `

  N重元组通过`gsl_ntuple`结构进行操作。这个包含存储有N重元组数据的文件以及指向当前N重元组数据行和用户定义的N冲远足数据结构尺寸信息:

  ```
  typedef struct {
      FILE *file;
      void *ntuple_data;
      size_t size;
  } gsl_ntuple;
  ```

  

## 生成N重元组

- [gsl_ntuple](https://www.gnu.org/software/gsl/doc/html/ntuple.html#c.gsl_ntuple) * `gsl_ntuple_create`(char * *filename*, void * *ntuple_data*, size_t *size*)

  This function creates a new write-only ntuple file `filename` for ntuples of size `size` and returns a pointer to the newly created ntuple struct. Any existing file with the same name is truncated to zero length and overwritten. A pointer to memory for the current ntuple row `ntuple_data` must be supplied—this is used to copy ntuples in and out of the file.

  这个函数为一个具有尺寸为`size`的N重元组生成一个新的只写N重元组文件`filename`并返回一个指向新产生的N重元组结构的指针。任何已存的同名文件将被删节为零长度并被重写。一个指向当前N重元组行内存的指针`ntuple_data`应当被提供—这是被用来将N重元组拷贝进以及出一个文件。

## 打开已存的N重元组文件

- [gsl_ntuple](https://www.gnu.org/software/gsl/doc/html/ntuple.html#c.gsl_ntuple) * `gsl_ntuple_open`(char * *filename*, void * *ntuple_data*, size_t *size*)

  This function opens an existing ntuple file `filename` for reading and returns a pointer to a corresponding ntuple struct. The ntuples in the file must have size `size`. A pointer to memory for the current ntuple row `ntuple_data` must be supplied—this is used to copy ntuples in and out of the file.

  这个函数打开一个已存的N重元组文件`filename`以进行读取病返回一个指向相应N重元组结构的指针。文件中的N重元组必须要具有尺寸`size`。一个指向当前N重元组行内存的指针`ntuple_data`应当被提供—这是被用来将N重元组拷贝进以及出一个文件。

## 写入N重元组

- int `gsl_ntuple_write`([gsl_ntuple](https://www.gnu.org/software/gsl/doc/html/ntuple.html#c.gsl_ntuple) * *ntuple*)

  This function writes the current ntuple `ntuple->ntuple_data` of size `ntuple->size` to the corresponding file.

- int `gsl_ntuple_bookdata`([gsl_ntuple](https://www.gnu.org/software/gsl/doc/html/ntuple.html#c.gsl_ntuple) * *ntuple*)

  This function is a synonym for [`gsl_ntuple_write()`](https://www.gnu.org/software/gsl/doc/html/ntuple.html#c.gsl_ntuple_write).

## Reading ntuples

- int `gsl_ntuple_read`([gsl_ntuple](https://www.gnu.org/software/gsl/doc/html/ntuple.html#c.gsl_ntuple) * *ntuple*)

  This function reads the current row of the ntuple file for `ntuple` and stores the values in `ntuple->data`.

## Closing an ntuple file

- int `gsl_ntuple_close`([gsl_ntuple](https://www.gnu.org/software/gsl/doc/html/ntuple.html#c.gsl_ntuple) * *ntuple*)

  This function closes the ntuple file `ntuple` and frees its associated allocated memory.

## Histogramming ntuple values

Once an ntuple has been created its contents can be histogrammed in various ways using the function [`gsl_ntuple_project()`](https://www.gnu.org/software/gsl/doc/html/ntuple.html#c.gsl_ntuple_project). Two user-defined functions must be provided, a function to select events and a function to compute scalar values. The selection function and the value function both accept the ntuple row as a first argument and other parameters as a second argument.



- `gsl_ntuple_select_fn`

  The *selection function* determines which ntuple rows are selected for histogramming. It is defined by the following struct:`typedef struct   {     int (* function) (void * ntuple_data, void * params);     void * params;   } gsl_ntuple_select_fn; `The struct component `function` should return a non-zero value for each ntuple row that is to be included in the histogram.



- `gsl_ntuple_value_fn`

  The *value function* computes scalar values for those ntuple rows selected by the selection function:`typedef struct   {     double (* function) (void * ntuple_data, void * params);     void * params;   } gsl_ntuple_value_fn; `In this case the struct component `function` should return the value to be added to the histogram for the ntuple row.



- int `gsl_ntuple_project`([gsl_histogram](https://www.gnu.org/software/gsl/doc/html/histogram.html#c.gsl_histogram) * *h*, [gsl_ntuple](https://www.gnu.org/software/gsl/doc/html/ntuple.html#c.gsl_ntuple) * *ntuple*, [gsl_ntuple_value_fn](https://www.gnu.org/software/gsl/doc/html/ntuple.html#c.gsl_ntuple_value_fn) * *value_func*, [gsl_ntuple_select_fn](https://www.gnu.org/software/gsl/doc/html/ntuple.html#c.gsl_ntuple_select_fn) * *select_func*)

  This function updates the histogram `h` from the ntuple `ntuple` using the functions `value_func` and `select_func`. For each ntuple row where the selection function `select_func`is non-zero the corresponding value of that row is computed using the function `value_func`and added to the histogram. Those ntuple rows where `select_func` returns zero are ignored. New entries are added to the histogram, so subsequent calls can be used to accumulate further data in the same histogram.

## Examples

The following example programs demonstrate the use of ntuples in managing a large dataset. The first program creates a set of 10,000 simulated “events”, each with 3 associated values ![(x,y,z)](https://www.gnu.org/software/gsl/doc/html/_images/math/62c54178349a8df2b5ba6d028eadd0cb60c837a6.png). These are generated from a Gaussian distribution with unit variance, for demonstration purposes, and written to the ntuple file `test.dat`.

```
#include <gsl/gsl_ntuple.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>

struct data
{
  double x;
  double y;
  double z;
};

int
main (void)
{
  const gsl_rng_type * T;
  gsl_rng * r;

  struct data ntuple_row;
  int i;

  gsl_ntuple *ntuple
    = gsl_ntuple_create ("test.dat", &ntuple_row,
                         sizeof (ntuple_row));

  gsl_rng_env_setup ();

  T = gsl_rng_default;
  r = gsl_rng_alloc (T);

  for (i = 0; i < 10000; i++)
    {
      ntuple_row.x = gsl_ran_ugaussian (r);
      ntuple_row.y = gsl_ran_ugaussian (r);
      ntuple_row.z = gsl_ran_ugaussian (r);

      gsl_ntuple_write (ntuple);
    }

  gsl_ntuple_close (ntuple);
  gsl_rng_free (r);

  return 0;
}
```

The next program analyses the ntuple data in the file `test.dat`. The analysis procedure is to compute the squared-magnitude of each event, ![E^2=x^2+y^2+z^2](https://www.gnu.org/software/gsl/doc/html/_images/math/b0a8513b799891d29d2f494dc5a5a8d5b42c2aaf.png), and select only those which exceed a lower limit of 1.5. The selected events are then histogrammed using their ![E^2](https://www.gnu.org/software/gsl/doc/html/_images/math/78972277721e8a672f1493e92a5e8f7037af44b7.png) values.

```
#include <math.h>
#include <gsl/gsl_ntuple.h>
#include <gsl/gsl_histogram.h>

struct data
{
  double x;
  double y;
  double z;
};

int sel_func (void *ntuple_data, void *params);
double val_func (void *ntuple_data, void *params);

int
main (void)
{
  struct data ntuple_row;

  gsl_ntuple *ntuple
    = gsl_ntuple_open ("test.dat", &ntuple_row,
                       sizeof (ntuple_row));
  double lower = 1.5;

  gsl_ntuple_select_fn S;
  gsl_ntuple_value_fn V;

  gsl_histogram *h = gsl_histogram_alloc (100);
  gsl_histogram_set_ranges_uniform(h, 0.0, 10.0);

  S.function = &sel_func;
  S.params = &lower;

  V.function = &val_func;
  V.params = 0;

  gsl_ntuple_project (h, ntuple, &V, &S);
  gsl_histogram_fprintf (stdout, h, "%f", "%f");
  gsl_histogram_free (h);
  gsl_ntuple_close (ntuple);

  return 0;
}

int
sel_func (void *ntuple_data, void *params)
{
  struct data * data = (struct data *) ntuple_data;
  double x, y, z, E2, scale;
  scale = *(double *) params;

  x = data->x;
  y = data->y;
  z = data->z;

  E2 = x * x + y * y + z * z;

  return E2 > scale;
}

double
val_func (void *ntuple_data, void *params)
{
  (void)(params); /* avoid unused parameter warning */
  struct data * data = (struct data *) ntuple_data;
  double x, y, z;

  x = data->x;
  y = data->y;
  z = data->z;

  return x * x + y * y + z * z;
}
```

[Fig. 14](https://www.gnu.org/software/gsl/doc/html/ntuple.html#fig-ntuples) shows the distribution of the selected events. Note the cut-off at the lower bound.



![_images/ntuple.png](https://www.gnu.org/software/gsl/doc/html/_images/ntuple.png)

Fig. 14 Distribution of selected events

## References and Further Reading

Further information on the use of ntuples can be found in the documentation for the CERN packages PAW and HBOOK (available online).