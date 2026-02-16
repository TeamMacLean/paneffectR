# Build presence/absence matrix from clustering results

Creates a matrix with orthogroups as rows and assemblies as columns. Can
produce binary (0/1), score-based, or count matrices.

## Usage

``` r
build_pa_matrix(
  orthogroups,
  proteins = NULL,
  type = "binary",
  score_column = "custom_score",
  score_threshold = NULL,
  score_aggregation = "max",
  exclude_singletons = FALSE
)
```

## Arguments

- orthogroups:

  An `orthogroup_result` object.

- proteins:

  Optional `protein_collection` with score data.

- type:

  Character. Matrix type: "binary", "count", or "score".

- score_column:

  Character. Column name for scores (default: "custom_score").

- score_threshold:

  Numeric. Only include proteins with score \>= threshold.

- score_aggregation:

  Character. How to aggregate scores: "max", "mean", or "sum".

- exclude_singletons:

  Logical. If FALSE (default), singletons are included as single-member
  orthogroups with IDs like "OG_single_001". If TRUE, singletons are
  excluded from the matrix.

## Value

A `pa_matrix` object.

## Examples

``` r
if (FALSE) { # \dontrun{
pa <- build_pa_matrix(clusters)
pa <- build_pa_matrix(clusters, score_threshold = 5.0)
pa <- build_pa_matrix(clusters, exclude_singletons = TRUE)
} # }
```
