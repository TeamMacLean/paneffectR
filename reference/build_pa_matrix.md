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
# Pre-computed clustering of the bundled example data
visual_dir <- system.file("testdata", "visual", package = "paneffectR")
clusters <- readRDS(file.path(visual_dir, "clusters_visual.rds"))
proteins <- readRDS(file.path(visual_dir, "proteins_visual.rds"))

pa <- build_pa_matrix(clusters)
pa
#> -- pa_matrix (binary) --
#> 62 orthogroups x 6 assemblies
#> Sparsity: 53.8%

# Drop singletons
build_pa_matrix(clusters, exclude_singletons = TRUE)
#> -- pa_matrix (binary) --
#> 50 orthogroups x 6 assemblies
#> Sparsity: 46.7%

# Keep only orthogroups with at least one high-scoring member
# (score filtering needs the protein_collection that carries the scores)
build_pa_matrix(clusters, proteins = proteins, score_threshold = 5.0)
#> -- pa_matrix (binary) --
#> 31 orthogroups x 6 assemblies
#> Sparsity: 32.3%
#> Threshold: 5
```
