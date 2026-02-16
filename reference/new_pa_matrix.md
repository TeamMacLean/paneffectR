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
