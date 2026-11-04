# Order rows by clustering their distinct presence/absence patterns

Rows with an identical pattern across columns are grouped; the distinct
patterns are clustered with
[`calculate_distance()`](https://TeamMacLean.github.io/paneffectR/reference/calculate_distance.md) +
[`hclust()`](https://rdrr.io/r/stats/hclust.html), and every row is
placed according to its pattern's position in that clustering.

## Usage

``` r
pattern_row_order(mat, method)
```

## Arguments

- mat:

  Numeric matrix without NAs (rows are observations).

- method:

  Distance method passed to
  [`calculate_distance()`](https://TeamMacLean.github.io/paneffectR/reference/calculate_distance.md).

## Value

A list with `order` (integer row order) and `pattern` (integer pattern
index per row, in original row order).
