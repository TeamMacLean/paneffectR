# Plot pan-genome structure

Creates a bar chart showing the composition of the pan-genome by
category: Core (present in all assemblies), Accessory, Rare, Unique
(present in one assembly), and optionally Singletons (proteins not in
any orthogroup).

## Usage

``` r
plot_pan_structure(
  clusters,
  pa = NULL,
  rare_threshold = 0.1,
  accessory_threshold = 0.6,
  exclude_singletons = FALSE,
  colors = NULL,
  show_counts = TRUE,
  title = NULL
)
```

## Arguments

- clusters:

  An `orthogroup_result` object.

- pa:

  Optional pre-built `pa_matrix`. If NULL, built internally.

- rare_threshold:

  Proportion threshold for "Rare" category (default 0.10). Orthogroups
  present in \<= this proportion of assemblies (but \> 1) are "Rare".

- accessory_threshold:

  Proportion threshold for "Accessory" category (default 0.60).
  Orthogroups present in \> rare_threshold and \<= this proportion are
  "Accessory".

- exclude_singletons:

  Logical. If TRUE, singletons are not shown (default FALSE).

- colors:

  Named vector of colors for each category. If NULL, uses default
  palette.

- show_counts:

  Logical. If TRUE, display count labels on bars (default TRUE).

- title:

  Plot title. If NULL, uses default.

## Value

A ggplot object.

## Examples

``` r
if (FALSE) { # \dontrun{
# Basic usage
plot_pan_structure(clusters)

# Custom thresholds
plot_pan_structure(clusters, rare_threshold = 0.15, accessory_threshold = 0.50)

# Without singletons
plot_pan_structure(clusters, exclude_singletons = TRUE)
} # }
```
