# Plot assembly clustering dendrogram

Clusters assemblies by shared orthogroup content and displays as a
dendrogram.

## Usage

``` r
plot_dendro(
  pa,
  distance_method = "jaccard",
  cluster_method = "complete",
  labels = TRUE,
  hang = 0.1
)
```

## Arguments

- pa:

  A `pa_matrix` object.

- distance_method:

  Character. Distance metric: "jaccard" (default), "binary", or
  "bray-curtis".

- cluster_method:

  Character. Clustering method for hclust (default: "complete").
  Options: "complete", "single", "average", "ward.D2".

- labels:

  Logical. Whether to show assembly labels (default: TRUE).

- hang:

  Numeric. Fraction of plot height to hang labels below leaves (default:
  0.1).

## Value

A ggplot2 object.

## Examples

``` r
if (FALSE) { # \dontrun{
plot_dendro(pa)
plot_dendro(pa, distance_method = "bray-curtis")
plot_dendro(pa, cluster_method = "ward.D2")
} # }
```
