# Create a pa_matrix object

Presence/absence matrix with metadata.

## Usage

``` r
new_pa_matrix(
  matrix,
  orthogroups,
  assemblies,
  type = "binary",
  threshold = NULL
)
```

## Arguments

- matrix:

  Matrix. Rows = orthogroups, cols = assemblies (0/1 or scores).

- orthogroups:

  A tibble. Orthogroup metadata.

- assemblies:

  A tibble. Assembly metadata.

- type:

  Character. One of "binary", "score", or "count".

- threshold:

  Numeric. Score threshold used (if applicable).

## Value

A `pa_matrix` S3 object.

## Examples

``` r
m <- matrix(c(1L, 1L, 1L, 0L), nrow = 2,
  dimnames = list(c("OG0001", "OG0002"), c("isolate_A", "isolate_B")))
pa <- new_pa_matrix(
  matrix = m,
  orthogroups = tibble::tibble(orthogroup_id = rownames(m), size = c(2L, 1L)),
  assemblies = tibble::tibble(assembly_name = colnames(m), n_orthogroups = c(2L, 1L))
)
pa
#> -- pa_matrix (binary) --
#> 2 orthogroups x 2 assemblies
#> Sparsity: 25%
```
