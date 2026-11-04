# Plot presence/absence heatmap

Creates a heatmap visualization of the presence/absence matrix using
ComplexHeatmap. Rows are orthogroups, columns are assemblies.

## Usage

``` r
plot_heatmap(
  pa,
  cluster_rows = TRUE,
  cluster_cols = TRUE,
  show_row_names = FALSE,
  show_col_names = TRUE,
  color = NULL,
  row_annotation = NULL,
  col_annotation = NULL,
  distance_method = "binary",
  max_cluster_rows = 5000,
  ...
)
```

## Arguments

- pa:

  A `pa_matrix` object.

- cluster_rows:

  Logical. Whether to cluster rows (default: TRUE).

- cluster_cols:

  Logical. Whether to cluster columns (default: TRUE).

- show_row_names:

  Logical. Whether to show row names (default: FALSE).

- show_col_names:

  Logical. Whether to show column names (default: TRUE).

- color:

  Character vector or color function. Colors for the heatmap. If NULL,
  auto-detected based on pa\$type.

- row_annotation:

  A data.frame with row names matching orthogroup IDs. Columns become
  annotation tracks on the left side.

- col_annotation:

  A data.frame with row names matching assembly names. Columns become
  annotation tracks on the top.

- distance_method:

  Character. Distance metric for clustering: "binary" (default),
  "jaccard", or "bray-curtis".

- max_cluster_rows:

  Integer. Above this many rows, row clustering is done on distinct
  patterns rather than individual rows (default 5000). See "Large
  matrices".

- ...:

  Additional arguments passed to ComplexHeatmap::Heatmap().

## Value

A ComplexHeatmap Heatmap object.

## Large matrices

Clustering every row of a matrix with tens of thousands of orthogroups
is slow and the resulting dendrogram is unreadable. When the matrix has
more than `max_cluster_rows` rows and `cluster_rows = TRUE`, rows are
instead ordered by clustering the distinct presence/absence patterns (at
most `2^n_assemblies`), which gives the same grouping in a fraction of
the time but draws no row dendrogram. For a readable plot, reduce the
matrix first, e.g. with
`build_pa_matrix(..., exclude_singletons = TRUE)` or a
`score_threshold`.

## Examples

``` r
visual_dir <- system.file("testdata", "visual", package = "paneffectR")
pa <- readRDS(file.path(visual_dir, "pa_binary.rds"))
ht <- plot_heatmap(pa)
ComplexHeatmap::draw(ht)


plot_heatmap(pa, cluster_rows = FALSE, distance_method = "jaccard")
```
