# Plot effector score distributions

Creates density plots or histograms of effector scores from a protein
collection.

## Usage

``` r
plot_scores(
  proteins,
  score_column = "custom_score",
  by_assembly = FALSE,
  threshold = NULL,
  plot_type = "density"
)
```

## Arguments

- proteins:

  A `protein_collection` object with score data.

- score_column:

  Character. Column name for scores (default: "custom_score").

- by_assembly:

  Logical. Whether to facet by assembly (default: FALSE). If TRUE,
  creates separate panels for each assembly.

- threshold:

  Numeric or NULL. If specified, draws a vertical line at this value to
  show where a score cutoff would fall.

- plot_type:

  Character. Plot type: "density" (default) or "histogram".

## Value

A ggplot2 object.

## Examples

``` r
if (FALSE) { # \dontrun{
plot_scores(proteins)
plot_scores(proteins, threshold = 5.0)
plot_scores(proteins, by_assembly = TRUE)
plot_scores(proteins, plot_type = "histogram")
} # }
```
