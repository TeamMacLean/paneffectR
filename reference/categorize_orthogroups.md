# Categorize orthogroups by presence pattern

Internal helper that assigns each orthogroup to a category based on how
many assemblies contain it.

## Usage

``` r
categorize_orthogroups(
  clusters,
  pa = NULL,
  rare_threshold = 0.1,
  accessory_threshold = 0.6
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
  Orthogroups present in \<= this proportion (but \> rare_threshold) are
  "Accessory".

## Value

A tibble with columns: orthogroup_id, n_present, n_assemblies, category
