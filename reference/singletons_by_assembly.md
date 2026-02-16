# Summarize singletons by assembly

Returns a summary of singleton counts per assembly, optionally including
the percentage of total proteins that are singletons.

## Usage

``` r
singletons_by_assembly(orthogroups, proteins = NULL)
```

## Arguments

- orthogroups:

  An `orthogroup_result` object.

- proteins:

  Optional `protein_collection` to calculate percentages.

## Value

A tibble with columns:

- `assembly`: Assembly name

- `n_singletons`: Number of singletons in assembly

- `n_total`: Total proteins in assembly (if proteins provided)

- `pct_singleton`: Percentage of proteins that are singletons (if
  proteins provided)

## Examples

``` r
if (FALSE) { # \dontrun{
summary <- singletons_by_assembly(clusters)
summary_with_pct <- singletons_by_assembly(clusters, proteins = proteins)
} # }
```
