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

- ...:

  Additional arguments passed to ComplexHeatmap::Heatmap().

## Value

A ComplexHeatmap Heatmap object.

## Examples

``` r
if (FALSE) { # \dontrun{
plot_heatmap(pa)
plot_heatmap(pa, cluster_rows = FALSE)
} # }
```
