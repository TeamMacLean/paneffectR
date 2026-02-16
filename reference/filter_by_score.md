# Filter presence/absence matrix by score

Removes orthogroups where no member meets the score threshold.

## Usage

``` r
filter_by_score(pa, proteins, threshold, score_column = "custom_score")
```

## Arguments

- pa:

  A `pa_matrix` object.

- proteins:

  A `protein_collection` with score data.

- threshold:

  Numeric. Minimum score to retain an orthogroup.

- score_column:

  Character. Column name for scores.

## Value

A filtered `pa_matrix` object.
