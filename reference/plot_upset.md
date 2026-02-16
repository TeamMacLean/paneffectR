# Plot UpSet diagram of orthogroup sharing

Creates an UpSet plot showing intersections of orthogroups across
assemblies. Useful for visualizing core, accessory, and unique
orthogroups.

## Usage

``` r
plot_upset(pa, min_size = 1, max_sets = NULL, order_by = "freq", ...)
```

## Arguments

- pa:

  A `pa_matrix` object.

- min_size:

  Integer. Minimum intersection size to show (default: 1).

- max_sets:

  Integer. Maximum number of sets to display. If NULL (default), all
  assemblies are shown.

- order_by:

  Character. How to order intersections: "freq" (by size, default) or
  "degree" (by number of sets in intersection).

- ...:

  Additional arguments passed to UpSetR::upset().

## Value

An UpSetR plot object (can be printed to display).

## Examples

``` r
if (FALSE) { # \dontrun{
plot_upset(pa)
plot_upset(pa, min_size = 2)
plot_upset(pa, order_by = "degree")
} # }
```
