# Get singletons from clustering result

Extracts the singletons tibble from an orthogroup_result object.

## Usage

``` r
get_singletons(orthogroups)
```

## Arguments

- orthogroups:

  An `orthogroup_result` object.

## Value

A tibble with columns `protein_id` and `assembly`.

## Examples

``` r
if (FALSE) { # \dontrun{
singletons <- get_singletons(clusters)
} # }
```
