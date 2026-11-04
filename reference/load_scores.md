# Load effector scores from CSV

Parses a scored CSV file, such as those produced by the omnieff
pipeline. The file must have `protein_id` and `custom_score` columns;
see the "Score file format" section of
[`load_proteins()`](https://TeamMacLean.github.io/paneffectR/reference/load_proteins.md)
for details.

## Usage

``` r
load_scores(path)
```

## Arguments

- path:

  Character. Path to scored CSV file.

## Value

A tibble with score columns.

## Examples

``` r
csv <- system.file("testdata", "assembly1_scored.csv", package = "paneffectR")
scores <- load_scores(csv)
scores[, c("protein_id", "custom_score", "score_rank")]
#> # A tibble: 100 × 3
#>    protein_id       custom_score score_rank
#>    <chr>                   <dbl>      <int>
#>  1 assembly1_000267         5.92          1
#>  2 assembly1_003502         5.14          2
#>  3 assembly1_001347         5.06          3
#>  4 assembly1_005840         5.06          4
#>  5 assembly1_006291         5.00          5
#>  6 assembly1_004421         4.99          6
#>  7 assembly1_001504         4.99          7
#>  8 assembly1_002717         4.95          8
#>  9 assembly1_002883         4.93          9
#> 10 assembly1_004234         4.91         10
#> # ℹ 90 more rows
```
