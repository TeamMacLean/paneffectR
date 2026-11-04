# Plot presence/absence heatmap

Creates a heatmap visualization of the presence/absence matrix using
ComplexHeatmap. Rows are orthogroups, columns are assemblies.

## Usage

``` r
plot_heatmap(
  pa,
  style = c("clean", "pan", "mono"),
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

- style:

  Character. One of `"clean"` (default), `"pan"` or `"mono"`. See
  "Styles".

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
  chosen from `style` and `pa$type`.

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

  Additional arguments passed to ComplexHeatmap::Heatmap(). These take
  precedence over the style's settings.

## Value

A ComplexHeatmap Heatmap object.

## Styles

`style` selects a ready-made look. All three use the same layout rules
(light grid between cells on matrices of up to 200 rows, slim
dendrograms, angled column labels, legends labelled "Absent"/"Present")
and differ in colour and in what the rows say:

- `"clean"` (default): one hue on a light grey background. Binary and
  count matrices are blue; score matrices are orange so that they are
  never mistaken for presence/absence.

- `"pan"`: rows are split into **Core** (present in every assembly),
  **Accessory** (present in more than one but not all) and **Unique**
  (present in exactly one) blocks, each with a coloured side bar. On
  binary matrices the cells themselves take the class colour; on count
  and score matrices the cells keep the magnitude ramp and the class is
  shown by the split. Note that
  [`plot_pan_structure()`](https://TeamMacLean.github.io/paneffectR/reference/plot_pan_structure.md)
  further divides Accessory into Accessory and Rare using thresholds;
  the heatmap does not.

- `"mono"`: black on white with grey grid lines, for print and greyscale
  figures.

`color` overrides the cell colours of any style. Count matrices whose
maximum is 1 are drawn with the binary colours.

## Customising further

The return value is a
[`ComplexHeatmap::Heatmap`](https://rdrr.io/pkg/ComplexHeatmap/man/Heatmap.html)
object, so anything that package can do applies. Arguments in `...` are
passed to
[`ComplexHeatmap::Heatmap()`](https://rdrr.io/pkg/ComplexHeatmap/man/Heatmap.html)
and override the style's settings, e.g. `column_title = "Effectors"`,
`rect_gp = grid::gpar(col = NA)` to remove the grid, or
`column_names_rot = 90`. See
[`vignette("styling-heatmaps")`](https://TeamMacLean.github.io/paneffectR/articles/styling-heatmaps.md).

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


# Core / accessory / unique blocks
plot_heatmap(pa, style = "pan")


# Greyscale for print, no row clustering
plot_heatmap(pa, style = "mono", cluster_rows = FALSE, distance_method = "jaccard")
```
