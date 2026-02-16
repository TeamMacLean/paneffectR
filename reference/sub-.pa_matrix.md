# Subset a pa_matrix object

Subset a pa_matrix object

## Usage

``` r
# S3 method for class 'pa_matrix'
x[i, j, drop = FALSE, ...]
```

## Arguments

- x:

  A `pa_matrix` object.

- i:

  Row indices (orthogroups).

- j:

  Column indices (assemblies).

- drop:

  Ignored (always FALSE to preserve pa_matrix structure).

- ...:

  Additional arguments (ignored).

## Value

A subsetted `pa_matrix` object.
