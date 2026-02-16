# Plot assembly composition

Creates a stacked bar chart showing the composition of each assembly's
proteins by category: Core, Accessory, Rare, Unique, and optionally
Singletons.

## Usage

``` r
plot_assembly_composition(
  clusters,
  pa = NULL,
  rare_threshold = 0.1,
  accessory_threshold = 0.6,
  exclude_singletons = FALSE,
  colors = NULL,
  position = "stack",
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

- position:

  Bar position: "stack" for counts (default), "fill" for proportions.

- title:

  Plot title. If NULL, uses default.

## Value

A ggplot object.

## Examples

``` r
if (FALSE) { # \dontrun{
# Basic usage
plot_assembly_composition(clusters)

# Show proportions instead of counts
plot_assembly_composition(clusters, position = "fill")

# Without singletons
plot_assembly_composition(clusters, exclude_singletons = TRUE)
} # }
```
